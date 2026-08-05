@timeit timer "Calc_Vnk_CWF" function Calc_Vnk!(
    material::CWF_model, kpoints::BoltzKPoints, Enk::Array{Float64,3},
    Cnk::Vector{Vector{Matrix{ComplexF64}}})

    Latvecs = material.Latvecs
    Nwann = Int(material.Ngsize)
    spinsize = material.SpinPol == "on" ? 2 : 1
    ChemP = material.ChemP
    NCell = Int(material.NCell)
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    cartesian_cell = zeros(Float64, 3, NCell)
    for cell = 1:NCell
        cartesian_cell[1,cell] = (Latvecs[1,1]*cell_list_ijk[cell][1] + Latvecs[2,1]*cell_list_ijk[cell][2] + Latvecs[3,1]*cell_list_ijk[cell][3])/Ang_to_bohr
        cartesian_cell[2,cell] = (Latvecs[1,2]*cell_list_ijk[cell][1] + Latvecs[2,2]*cell_list_ijk[cell][2] + Latvecs[3,2]*cell_list_ijk[cell][3])/Ang_to_bohr
        cartesian_cell[3,cell] = (Latvecs[1,3]*cell_list_ijk[cell][1] + Latvecs[2,3]*cell_list_ijk[cell][2] + Latvecs[3,3]*cell_list_ijk[cell][3])/Ang_to_bohr
    end

    Vnk = zeros(Float64, Nwann, 3, spinsize, MPI_Nkpt)

    # Keep these matrices as explicitly named locals.  Using an `ntuple`
    # closure here boxes Nwann on Julia 1.12 and makes the innermost Fourier
    # loops dynamically dispatched.
    Hmnk_x = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_y = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_z = zeros(ComplexF64, Nwann, Nwann)
    Htemp1 = zeros(ComplexF64, Nwann, Nwann)
    Htemp2 = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, local_k = 1:MPI_Nkpt
        fill!(Hmnk_x, 0.0)
        fill!(Hmnk_y, 0.0)
        fill!(Hmnk_z, 0.0)

        kpoint = MPI_kpts[local_k]
        for cell = 1:NCell
            kRn = dot(kpoint, cell_list_ijk[cell])
            phase = cispi(2 * kRn)
            prefactor_x = -im*cartesian_cell[1,cell]*eV2Hartree*phase
            prefactor_y = -im*cartesian_cell[2,cell]*eV2Hartree*phase
            prefactor_z = -im*cartesian_cell[3,cell]*eV2Hartree*phase

            @inbounds for ist = 1:Nwann, jst = 1:Nwann
                hvalue = HmnR[jst,ist,cell,spin]
                Hmnk_x[ist, jst] += hvalue*prefactor_x
                Hmnk_y[ist, jst] += hvalue*prefactor_y
                Hmnk_z[ist, jst] += hvalue*prefactor_z
            end
        end

        Cn = Cnk[spin][local_k]

        mul!(Htemp1, adjoint(Cn), Hmnk_x)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band,1,spin,local_k] = real(Htemp2[band,band])
        end

        mul!(Htemp1, adjoint(Cn), Hmnk_y)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band,2,spin,local_k] = real(Htemp2[band,band])
        end

        mul!(Htemp1, adjoint(Cn), Hmnk_z)
        mul!(Htemp2, Htemp1, Cn)
        @inbounds for band = 1:Nwann
            Vnk[band,3,spin,local_k] = real(Htemp2[band,band])
        end
    end

    @. Enk = Enk*eV2Hartree - ChemP*eV2Hartree

    
    return Vnk
end
