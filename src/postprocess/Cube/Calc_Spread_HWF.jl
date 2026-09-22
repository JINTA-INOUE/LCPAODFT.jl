function Calc_Spread_HWF(material::LCPAO_model, kpoints::KPoints, Wannier_Center, OLPpos2, Cnk, Umjk, sflag)

    SpinPol = material.SpinPol
    if SpinPol ∈ ("off", "on")
        Calc_Spread_HWF_Col(material, kpoints, Wannier_Center, OLPpos2, Cnk, Umjk, sflag)
    else
        Calc_Spread_HWF_NonCol(material, kpoints, Wannier_Center, OLPpos2, Cnk, Umjk, sflag)
    end
end


function Calc_Spread_HWF_Col(material::LCPAO_model, kpoints::KPoints, Wannier_Center, OLPpos2, Cnk, Umjk, sflag)

    Natom = material.Natom
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)

    kmesh = kpoints.kmesh
    AllNkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    NkAB = Int64(div(AllNkpt,kmesh[sflag]))

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    else
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    end

    Wannier_Center2x = Vector{Vector{Float64}}(undef, spinsize)
    Wannier_Center2y = Vector{Vector{Float64}}(undef, spinsize)
    Wannier_Center2z = Vector{Vector{Float64}}(undef, spinsize)
    for spin = 1:spinsize
        Wannier_Center2x[spin] = zeros(Float64, Nocc)
        Wannier_Center2y[spin] = zeros(Float64, Nocc)
        Wannier_Center2z[spin] = zeros(Float64, Nocc)
    end
    

    COLP = zeros(ComplexF64, Nfsize, Nfsize)
    COLPC = zeros(ComplexF64, Nfsize, Nfsize)
    COLPC_occ = zeros(ComplexF64, Nocc, Nocc)
    UCOLPC = zeros(ComplexF64, Nocc, Nocc)
    UCOLPCU = zeros(ComplexF64, Nocc, Nocc)
    Sx = zeros(ComplexF64, fsize, fsize)
    Sy = zeros(ComplexF64, fsize, fsize)
    Sz = zeros(ComplexF64, fsize, fsize)


    for spin = 1:spinsize, ik = 1:NkAB
        HS_matrix!(Sx, OLPpos2[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
    
        mul!(COLP, Cnk[spin][ik]', Sx)
        mul!(COLPC, COLP, Cnk[spin][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[spin][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[spin][ik])

        for μ = 1:Nocc
            # Wannier_Center2x[spin][μ] += UCOLPCU[μ,μ]/NkAB
            Wannier_Center2x[spin][μ] += real(UCOLPCU[μ,μ])/NkAB
        end

        HS_matrix!(Sy, OLPpos2[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])

        mul!(COLP, Cnk[spin][ik]', Sy)
        mul!(COLPC, COLP, Cnk[spin][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[spin][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[spin][ik])

        for μ = 1:Nocc
            Wannier_Center2y[spin][μ] += real(UCOLPCU[μ,μ])/NkAB
        end

        HS_matrix!(Sz, OLPpos2[3], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])

        mul!(COLP, Cnk[spin][ik]', Sz)
        mul!(COLPC, COLP, Cnk[spin][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[spin][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[spin][ik])

        for μ = 1:Nocc
            Wannier_Center2z[spin][μ] += real(UCOLPCU[μ,μ])/NkAB
        end
    end


    for spin = 1:spinsize
        @show spin
        for μ = 1:Nocc
            @show μ, Wannier_Center[spin][μ], Wannier_Center2z[spin][μ]
        end
    end
end


function Calc_Spread_HWF_NonCol(material::LCPAO_model, kpoints::KPoints, Wannier_Center, OLPpos2, Cnk, Umjk, sflag)

    Natom = material.Natom
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize
    Nocc = Int64(Valence_Electrons)

    kmesh = kpoints.kmesh
    AllNkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    NkAB = Int64(div(AllNkpt,kmesh[sflag]))


    Wannier_Center2x = Vector{Vector{Float64}}(undef, spinsize)
    Wannier_Center2y = Vector{Vector{Float64}}(undef, spinsize)
    Wannier_Center2z = Vector{Vector{Float64}}(undef, spinsize)
    for spin = 1:spinsize
        Wannier_Center2x[spin] = zeros(Float64, Nocc)
        Wannier_Center2y[spin] = zeros(Float64, Nocc)
        Wannier_Center2z[spin] = zeros(Float64, Nocc)
    end
    

    COLP = zeros(ComplexF64, Nfsize, Nfsize)
    COLPC = zeros(ComplexF64, Nfsize, Nfsize)
    COLPC_occ = zeros(ComplexF64, Nocc, Nocc)
    UCOLPC = zeros(ComplexF64, Nocc, Nocc)
    UCOLPCU = zeros(ComplexF64, Nocc, Nocc)
    tmpH = zeros(ComplexF64, fsize, fsize)
    Sx = zeros(ComplexF64, Nfsize, Nfsize)
    Sy = zeros(ComplexF64, Nfsize, Nfsize)
    Sz = zeros(ComplexF64, Nfsize, Nfsize)


    for ik = 1:NkAB
        HS_matrix!(tmpH, OLPpos2[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
        @. @views Sx[1:fsize, 1:fsize] = tmpH
        @. @views Sx[fsize+1:end, fsize+1:end] = tmpH

        mul!(COLP, Cnk[1][ik]', Sx)
        mul!(COLPC, COLP, Cnk[1][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[1][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[1][ik])

        for μ = 1:Nocc
            # Wannier_Center2x[1][μ] += UCOLPCU[μ,μ]/NkAB
            Wannier_Center2x[1][μ] += real(UCOLPCU[μ,μ])/NkAB
        end

        HS_matrix!(tmpH, OLPpos2[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
        @. @views Sy[1:fsize, 1:fsize] = tmpH
        @. @views Sy[fsize+1:end, fsize+1:end] = tmpH

        mul!(COLP, Cnk[1][ik]', Sy)
        mul!(COLPC, COLP, Cnk[1][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[1][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[1][ik])

        for μ = 1:Nocc
            Wannier_Center2y[1][μ] += real(UCOLPCU[μ,μ])/NkAB
        end

        HS_matrix!(tmpH, OLPpos2[3], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
        @. @views Sz[1:fsize, 1:fsize] = tmpH
        @. @views Sz[fsize+1:end, fsize+1:end] = tmpH
        
        mul!(COLP, Cnk[1][ik]', Sz)
        mul!(COLPC, COLP, Cnk[1][ik])
        @. @views COLPC_occ = COLPC[1:Nocc,1:Nocc]
        mul!(UCOLPC, Umjk[1][ik]', COLPC_occ)
        mul!(UCOLPCU, UCOLPC, Umjk[1][ik])

        for μ = 1:Nocc
            Wannier_Center2z[1][μ] += real(UCOLPCU[μ,μ])/NkAB
        end
    end


    for μ = 1:Nocc
        z = Wannier_Center[1][μ]/Ang_to_bohr
        z2 = Wannier_Center2z[1][μ]/Ang_to_bohr/Ang_to_bohr
        Ωz = z2 - z^2
        @show μ, z, z2, Ωz
    end
end