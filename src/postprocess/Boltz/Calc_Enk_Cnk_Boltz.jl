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
