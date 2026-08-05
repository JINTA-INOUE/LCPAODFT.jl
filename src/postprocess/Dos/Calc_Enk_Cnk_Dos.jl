function Calc_Enk_Cnk_Dos(material::LCPAO_model, kpoints::DosKPoints,
                          iemin, iemax)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol in ("off", "on"), fsize, 2 * fsize)
    spinsize = ifelse(SpinPol == "off", 1, 2)
    cnk_spinsize = ifelse(SpinPol == "on", 2, 1)
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks
    ChemP = material.ChemP

    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    MPI_kpts = kpoints.MPI_kpts
    neg = iemax - iemin + 1

    # The final dimension contains only k points owned by this MPI rank.
    Enk = zeros(Float64, neg, spinsize, MPI_Nkpt)
    Cnk = [Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
           for _ = 1:cnk_spinsize]

    if SpinPol in ("off", "on")
        S = zeros(ComplexF64, Nfsize, Nfsize)
        H = zeros(ComplexF64, Nfsize, Nfsize)

        for spin = 1:spinsize, local_k = 1:MPI_Nkpt
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP,
                       FNAN, natn, ncn, atv_ijk, MPI_kpts[local_k])
            decomposition = eigen(Hermitian(H), Hermitian(S))
            Cnk[spin][local_k] = decomposition.vectors
            @views Enk[:, spin, local_k] .=
                decomposition.values[iemin:iemax] .- ChemP
        end
    elseif SpinPol == "nc"
        S = zeros(ComplexF64, Nfsize, Nfsize)
        H = zeros(ComplexF64, Nfsize, Nfsize)
        tmpH = zeros(ComplexF64, fsize, fsize)

        for local_k = 1:MPI_Nkpt
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP,
                          FNAN, natn, ncn, atv_ijk, MPI_kpts[local_k])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP,
                       FNAN, natn, ncn, atv_ijk, MPI_kpts[local_k])
            fill!(S, 0.0)
            @views S[1:fsize, 1:fsize] .= tmpH
            @views S[fsize + 1:end, fsize + 1:end] .= tmpH

            decomposition = eigen(Hermitian(H), Hermitian(S))
            Cnk[1][local_k] = decomposition.vectors
            @views Enk[:, 1, local_k] .=
                decomposition.values[iemin:iemax] .- ChemP
        end
    else
        error("unsupported spin polarization: $SpinPol")
    end

    if SpinPol == "nc"
        @views Enk[:, 2, :] .= Enk[:, 1, :]
    end

    return Enk, Cnk
end


function Calc_Enk_Cnk_Dos(material::CWF_model, kpoints::DosKPoints,
                          iemin, iemax)

    SpinPol = material.SpinPol
    cnk_spinsize = ifelse(SpinPol == "on", 2, 1)
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Nwann = material.Ngsize
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    ChemP = material.ChemP

    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    MPI_kpts = kpoints.MPI_kpts
    neg = iemax - iemin + 1

    Enk = zeros(Float64, neg, spinsize, MPI_Nkpt)
    Cnk = [Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
           for _ = 1:cnk_spinsize]

    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:cnk_spinsize, local_k = 1:MPI_Nkpt
        fill!(H, 0.0)

        k1, k2, k3 = MPI_kpts[local_k]
        @inbounds for cell = 1:NCell
            cell_ijk = cell_list_ijk[cell]
            kRn = k1 * cell_ijk[1] + k2 * cell_ijk[2] + k3 * cell_ijk[3]
            ex = cispi(2 * kRn)
            for ist = 1:Nwann, jst = 1:Nwann
                H[ist, jst] += HmnR[jst, ist, cell, spin] * ex
            end
        end

        decomposition = eigen(Hermitian(H))
        Cnk[spin][local_k] = decomposition.vectors
        @views Enk[:, spin, local_k] .=
            decomposition.values[iemin:iemax] .- ChemP
    end

    if SpinPol == "nc"
        @views Enk[:, 2, :] .= Enk[:, 1, :]
    end

    return Enk, Cnk
end
