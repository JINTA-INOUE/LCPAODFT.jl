"""
Type-stable CWF eigensystem calculation for the Boltz solver.

The original implementation remains in `Calc_Enk_Cnk_Boltz.jl` and is used
as the fallback for non-CWF models.  The explicit construction of `Cnk`
avoids boxing `MPI_Nkpt` in a comprehension closure on Julia 1.12.
"""
@timeit timer "Calc_Enk_Cnk_Boltz" function Calc_Enk_Cnk_Boltz(material::CWF_model, kpoints::BoltzKPoints)

    Nwann = Int(material.Ngsize)
    NCell = Int(material.NCell)
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    spinsize = material.SpinPol == "on" ? 2 : 1
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    Enk = zeros(Float64, Nwann, spinsize, MPI_Nkpt)

    # Do not use a comprehension that captures MPI_Nkpt here.  On Julia
    # 1.12 that capture boxes MPI_Nkpt and makes all subsequent k-point loop
    # indices type-unstable.
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
    end

    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, local_k = 1:MPI_Nkpt
        fill!(H, 0.0)
        k1, k2, k3 = MPI_kpts[local_k]
        for cell = 1:NCell
            kRn = k1 * cell_list_ijk[cell][1] +
                  k2 * cell_list_ijk[cell][2] +
                  k3 * cell_list_ijk[cell][3]
            phase = cispi(2 * kRn)
            @inbounds for ist = 1:Nwann, jst = 1:Nwann
                H[ist, jst] += HmnR[jst, ist, cell, spin] * phase
            end
        end

        decomposition = eigen(Hermitian(H))
        Cnk[spin][local_k] = decomposition.vectors
        @views Enk[:, spin, local_k] .= decomposition.values
    end

    return Enk, Cnk
end


"""
Calculate the band energies, velocities and optional orbital weights while
holding the eigensystem of only one k point at a time.

`Val(true)` requests the orbital weights used by the decomposed transport
calculation.  Using `Val` keeps the hot k-point loop type-stable in both
modes.  In contrast to `Calc_Enk_Cnk_Boltz`, this routine never builds the
`O(Nwann^2 * Nkpt)` collection of eigenvector matrices.
"""
function Calc_Enk_Vnk_EVec_Boltz(
    material::CWF_model, kpoints::BoltzKPoints, decomp::Bool)

    return _Calc_Enk_Vnk_EVec_Boltz(
        material, kpoints.MPI_kpts, Val(decomp))
end


"""Block-oriented entry point accepting an arbitrary list of k points."""
function Calc_Enk_Vnk_EVec_Boltz(
    material::CWF_model, kpoints::AbstractVector, decomp::Bool)

    return _Calc_Enk_Vnk_EVec_Boltz(material, kpoints, Val(decomp))
end


struct BoltzEigensystemWorkspace{E}
    Enk::Array{Float64,3}
    Vnk::Array{Float64,4}
    EVec::E
    cell_ijk::Matrix{Int32}
    derivative_base::Matrix{ComplexF64}
    H::Vector{Matrix{ComplexF64}}
    Hmnk_x::Vector{Matrix{ComplexF64}}
    Hmnk_y::Vector{Matrix{ComplexF64}}
    Hmnk_z::Vector{Matrix{ComplexF64}}
    Htemp1::Vector{Matrix{ComplexF64}}
    Htemp2::Vector{Matrix{ComplexF64}}
end


function _boltz_eigensystem_workspace(
    material::CWF_model, capacity::Integer,
    ::Val{DECOMP}) where {DECOMP}

    Latvecs = material.Latvecs
    Nwann = Int(material.Ngsize)
    NCell = Int(material.NCell)
    cell_list_ijk = material.cell_list_ijk
    spinsize = material.SpinPol == "on" ? 2 : 1
    ncapacity = Int(capacity)
    ncapacity > 0 || error("eigensystem workspace capacity must be positive")

    Enk = Array{Float64}(undef, Nwann, spinsize, ncapacity)
    Vnk = Array{Float64}(undef, Nwann, 3, spinsize, ncapacity)
    EVec = DECOMP ?
        Array{Float32}(undef, Nwann, Nwann, spinsize, ncapacity) : nothing

    cartesian_cell = Matrix{Float64}(undef, 3, NCell)
    cell_ijk = Matrix{Int32}(undef, 3, NCell)
    derivative_base = Matrix{ComplexF64}(undef, 3, NCell)
    for cell = 1:NCell
        i, j, k = cell_list_ijk[cell]
        cell_ijk[1, cell] = i
        cell_ijk[2, cell] = j
        cell_ijk[3, cell] = k
        cartesian_cell[1, cell] =
            (Latvecs[1, 1] * i + Latvecs[2, 1] * j +
             Latvecs[3, 1] * k) / Ang_to_bohr
        cartesian_cell[2, cell] =
            (Latvecs[1, 2] * i + Latvecs[2, 2] * j +
             Latvecs[3, 2] * k) / Ang_to_bohr
        cartesian_cell[3, cell] =
            (Latvecs[1, 3] * i + Latvecs[2, 3] * j +
             Latvecs[3, 3] * k) / Ang_to_bohr
        derivative_base[1, cell] =
            -im * cartesian_cell[1, cell] * eV2Hartree
        derivative_base[2, cell] =
            -im * cartesian_cell[2, cell] * eV2Hartree
        derivative_base[3, cell] =
            -im * cartesian_cell[3, cell] * eV2Hartree
    end

    # Each Julia thread gets private O(Nwann^2) work arrays.  The much larger
    # k-resolved output arrays remain shared and every thread writes disjoint
    # (spin, k-point) slices.
    thread_slots = Threads.maxthreadid()
    H = [zeros(ComplexF64, Nwann, Nwann) for _ = 1:thread_slots]
    Hmnk_x = [similar(H[1]) for _ = 1:thread_slots]
    Hmnk_y = [similar(H[1]) for _ = 1:thread_slots]
    Hmnk_z = [similar(H[1]) for _ = 1:thread_slots]
    Htemp1 = [similar(H[1]) for _ = 1:thread_slots]
    Htemp2 = [similar(H[1]) for _ = 1:thread_slots]

    return BoltzEigensystemWorkspace(
        Enk, Vnk, EVec, cell_ijk, derivative_base,
        H, Hmnk_x, Hmnk_y, Hmnk_z, Htemp1, Htemp2)
end


function _fill_boltz_eigensystem_impl!(
    workspace::BoltzEigensystemWorkspace, material::CWF_model,
    MPI_kpts::AbstractVector, ::Val{DECOMP}) where {DECOMP}

    Nwann = Int(material.Ngsize)
    NCell = Int(material.NCell)
    HmnR = material.HmnR
    spinsize = material.SpinPol == "on" ? 2 : 1
    ChemP = material.ChemP
    MPI_Nkpt = length(MPI_kpts)
    MPI_Nkpt <= size(workspace.Enk, 3) || error(
        "eigensystem workspace capacity exceeded")

    Enk = workspace.Enk
    Vnk = workspace.Vnk
    EVec = workspace.EVec
    cell_ijk = workspace.cell_ijk
    derivative_base = workspace.derivative_base
    matrix_size = Nwann * Nwann

    Threads.@threads :static for work_index = 1:(spinsize * MPI_Nkpt)
        spin = div(work_index - 1, MPI_Nkpt) + 1
        local_k = mod(work_index - 1, MPI_Nkpt) + 1
        thread = Threads.threadid()
        H = workspace.H[thread]
        Hmnk_x = workspace.Hmnk_x[thread]
        Hmnk_y = workspace.Hmnk_y[thread]
        Hmnk_z = workspace.Hmnk_z[thread]
        Htemp1 = workspace.Htemp1[thread]
        Htemp2 = workspace.Htemp2[thread]
        fill!(H, 0.0)
        fill!(Hmnk_x, 0.0)
        fill!(Hmnk_y, 0.0)
        fill!(Hmnk_z, 0.0)

        k1, k2, k3 = MPI_kpts[local_k]
        for cell = 1:NCell
            i = cell_ijk[1, cell]
            j = cell_ijk[2, cell]
            k = cell_ijk[3, cell]
            kRn = k1 * i + k2 * j + k3 * k
            phase = cispi(2 * kRn)
            prefactor_x = derivative_base[1, cell] * phase
            prefactor_y = derivative_base[2, cell] * phase
            prefactor_z = derivative_base[3, cell] * phase

            hcell_offset = ((spin - 1) * NCell + cell - 1) * matrix_size
            for ist = 1:Nwann
                hindex = hcell_offset + (ist - 1) * Nwann
                matrix_index = ist
                @inbounds for jst = 1:Nwann
                    hvalue = HmnR[hindex + jst]
                    H[matrix_index] += hvalue * phase
                    Hmnk_x[matrix_index] += hvalue * prefactor_x
                    Hmnk_y[matrix_index] += hvalue * prefactor_y
                    Hmnk_z[matrix_index] += hvalue * prefactor_z
                    matrix_index += Nwann
                end
            end
        end

        decomposition = eigen!(Hermitian(H))
        Cn = decomposition.vectors
        @inbounds for band = 1:Nwann
            Enk[band, spin, local_k] =
                decomposition.values[band] * eV2Hartree -
                ChemP * eV2Hartree
        end

        mul!(Htemp1, adjoint(Cn), Hmnk_x)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band, 1, spin, local_k] = real(Htemp2[band, band])
        end

        mul!(Htemp1, adjoint(Cn), Hmnk_y)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band, 2, spin, local_k] = real(Htemp2[band, band])
        end

        mul!(Htemp1, adjoint(Cn), Hmnk_z)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band, 3, spin, local_k] = real(Htemp2[band, band])
        end

        if DECOMP
            @inbounds for band = 1:Nwann, orbital = 1:Nwann
                coefficient = Cn[orbital, band]
                EVec[orbital, band, spin, local_k] =
                    Float32(abs2(coefficient))
            end
        end
    end

    return Enk, Vnk, EVec
end


@timeit timer "Calc_Enk_Vnk_EVec_Boltz" function _fill_boltz_eigensystem!(
    workspace::BoltzEigensystemWorkspace, material::CWF_model,
    MPI_kpts::AbstractVector, decomp::Val{DECOMP}) where {DECOMP}

    return _fill_boltz_eigensystem_impl!(
        workspace, material, MPI_kpts, decomp)
end


function _Calc_Enk_Vnk_EVec_Boltz(
    material::CWF_model, MPI_kpts::AbstractVector,
    decomp::Val{DECOMP}) where {DECOMP}

    workspace = _boltz_eigensystem_workspace(
        material, length(MPI_kpts), decomp)
    return _fill_boltz_eigensystem!(
        workspace, material, MPI_kpts, decomp)
end
