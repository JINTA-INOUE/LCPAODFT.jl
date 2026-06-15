function Get_Angle_spin(material::LCPAO_model)

    SpinPol = material.SpinPol
    if SpinPol ≠ "nc"
        error("not support $SpinPol case")
    end

    Nspin = 4
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    Total_NumOrbs = material.Total_NumOrbs
    OLP = material.OLP
    DM = material.DM


    MulP = zeros(Float64, 4)
    Angle_Spin = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Angle_Spin[atom] = zeros(Float64, 2)
    end

    for atom = 1:Natom

        fill!(MulP, 0.0)
        for Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            for spin = 1:Nspin, ist = 1:Total_NumOrbs[atom]
                tmp0 = 0.0
                for jst = 1:Total_NumOrbs[jatom]
                    tmp0 += DM[spin][atom][Rn][ist][jst] * OLP[atom][Rn][ist][jst]
                end

                if spin == 4
                    MulP[spin] -= tmp0
                else
                    MulP[spin] += tmp0
                end
            end
        end


        Nup, Ndown, theta, phi = EulerAngle_Spin(MulP[1], MulP[2], MulP[3], MulP[4])

        MulP[1] = Nup
        MulP[2] = Ndown
        MulP[3] = theta
        MulP[4] = phi

        Angle_Spin[atom][1] = MulP[3]
        Angle_Spin[atom][2] = MulP[4]
    end


    return Angle_Spin
end


function Calc_Band_size(material::LCPAO_model, Dos_Erange)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol ∈ ("off", "on"), fsize, 2*fsize)
    spinsize = ifelse(SpinPol == "off", 1, 2)

    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks
    ChemP = material.ChemP


    # Find iemin, iemax at Γ point
    # iemin : minimal band index
    # iemax : maximum band index
    iemin = 1
    iemax = 1
    n1min = 1
    kpts_zeros = zeros(Float64, 3)
    EΓ = zeros(Float64, Nfsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    tmpH = zeros(ComplexF64, fsize, fsize)

    for spin = 1:spinsize

        if SpinPol ∈ ("off", "on")
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        elseif SpinPol == "nc"
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        end

        iemin0 = 1
        n1min = ifelse(n1min<Nfsize, Nfsize, n1min)
        
        for μ = 1:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[1])
                iemin0 = μ - 1
                break
            end
        end

        iemin0 = ifelse(iemin0<1, 1, iemin0)
        iemax0 = Nfsize

        for μ = iemin0:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[2])
                iemax0 = μ
                break
            end
        end

        iemax0 = ifelse(iemax0>Nfsize, Nfsize, iemax0)
        iemin = ifelse(iemin>iemin0, iemin0, iemin)
        iemax = ifelse(iemax<iemax0, iemax0, iemax)
    end

    if SpinPol ∈ ("off", "on")
        iemin -= max(div(fsize, 20), 10)
        iemax += max(div(fsize, 20), 10)
    elseif SpinPol == "nc"
        iemin -= max(div(fsize, 10), 10)
        iemax += max(div(fsize, 10), 10)
    end

    iemin = ifelse(iemin<1, 1, iemin)
    iemax = ifelse(iemax>Nfsize, Nfsize, iemax)
    iemax = ifelse(iemax>n1min, n1min, iemax)

    return iemin, iemax
end


function Calc_Band_size(material::CWF_model, Dos_Erange)

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    NCell = material.NCell
    HmnR = material.HmnR
    ChemP = material.ChemP
    spinsize = ifelse(SpinPol == "on", 2, 1)    


    # Find iemin, iemax at Γ point
    # iemin : minimal band index
    # iemax : maximum band index
    iemin = 1
    iemax = 1
    n1min = 1
    EΓ = zeros(Float64, Nwann)
    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize
        for cell = 1:NCell, ist = 1:Nwann, jst = 1:Nwann
            H[ist,jst] += HmnR[jst,ist,cell,spin]
        end

        EΓ = eigvals(Hermitian(H))

        iemin0 = 1
        n1min = ifelse(n1min<Nwann, Nwann, n1min)
        
        for μ = 1:Nwann
            if EΓ[μ] > (ChemP + Dos_Erange[1])
                iemin0 = μ - 1
                break
            end
        end

        iemin0 = ifelse(iemin0<1, 1, iemin0)
        iemax0 = Nwann

        for μ = iemin0:Nwann
            if EΓ[μ] > (ChemP + Dos_Erange[2])
                iemax0 = μ
                break
            end
        end

        iemax0 = ifelse(iemax0>Nwann, Nwann, iemax0)
        iemin = ifelse(iemin>iemin0, iemin0, iemin)
        iemax = ifelse(iemax<iemax0, iemax0, iemax)
    end

    if SpinPol ∈ ("off", "on")
        iemin -= max(div(Nwann, 20), 10)
        iemax += max(div(Nwann, 20), 10)
    elseif SpinPol == "nc"
        iemin -= max(div(div(Nwann,2), 10), 10)
        iemax += max(div(div(Nwann,2), 10), 10)
    end

    iemin = ifelse(iemin<1, 1, iemin)
    iemax = ifelse(iemax>Nwann, Nwann, iemax)
    iemax = ifelse(iemax>n1min, n1min, iemax)

    return iemin, iemax
end


function Calc_PDos_Atom_proj(filename::String, material::LCPAO_model, Dos_Erange, DosE, Dos)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    spinsize = ifelse(SpinPol == "off", 1, 2)
    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]
    de = 0.01/eV2Hartree
    Dos_N = floor(Int, (DosEmax-DosEmin)/de)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)


    for atom = 1:Natom
        fill!(DosSum, 0.0)
        for spin = 1:spinsize
            for q = 1:Dos_N
                s1 = 0.0
                s2 = 0.0
                for ie = 2:2:q-1
                    DosSum[ie,spin] = 0.0
                    for ist = 1:Total_NumOrbs[atom]
                        DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                    end
                    s1 += DosSum[ie,spin]
                end

                for ie = 3:2:q-1
                    DosSum[ie,spin] = 0.0
                    for ist = 1:Total_NumOrbs[atom]
                        DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                    end
                    s2 += DosSum[ie,spin]
                end
                ssum[q,spin] = (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
            end
        end

        println("Write $(filename).PDOS.Tetrahedron.atom$(atom)")
        Dos_data = open("$(filename).PDOS.Tetrahedron.atom$(atom)", "w")

        for ie = 1:Dos_N
            fill!(DSum, 0.0)
            for ist = 1:Total_NumOrbs[atom], spin = 1:spinsize
                DSum[spin] += Dos[spin][atom][ist][ie]
            end

            if SpinPol ∈ ("on", "nc")
                @printf(Dos_data, "%lf %lf %lf %lf %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
            else
                @printf(Dos_data, "%lf %lf %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
            end
        end
        close(Dos_data)
    end
end


function Calc_PDos_Orbital_proj(filename::String, material::LCPAO_model, Spe_Num_Relation, Spe_Num_Basis, Dos_Erange, DosE, Dos)
    
    Natom = material.Natom
    SpinPol = material.SpinPol
    atom2spe = material.atom2spe
    Total_NumOrbs = material.Total_NumOrbs
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Lname = ["s", "p", "d", "f"]
    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]
    de = 0.01/eV2Hartree
    Dos_N = floor(Int, (DosEmax-DosEmin)/de)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)

    

    for atom = 1:Natom, L = 0:3
        spe = atom2spe[atom]
        if Spe_Num_Basis[spe][L+1] > 0
            for M = 0:2*L
                LM = 100*L + M + 1

                fill!(DosSum, 0.0)
                for spin = 1:spinsize
                    for q = 1:Dos_N
                        s1 = 0.0
                        s2 = 0.0
                        for ie = 2:2:q-1
                            DosSum[ie,spin] = 0.0
                            for ist = 1:Total_NumOrbs[atom]
                                if LM == Spe_Num_Relation[spe][ist]
                                    DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                                end
                            end
                            s1 += DosSum[ie,spin]
                        end
                        for ie = 3:2:q-1
                            DosSum[ie,spin] = 0.0
                            for ist = 1:Total_NumOrbs[atom]
                                if LM == Spe_Num_Relation[spe][ist]
                                    DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                                end
                            end
                            s2 += DosSum[ie,spin]
                        end
                        ssum[q,spin] = (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
                    end
                end

                println("Write $(filename).PDOS.Tetrahedron.atom$(atom).$(Lname[L+1])$(M+1)")
                PDos_data = open("$(filename).PDOS.Tetrahedron.atom$(atom).$(Lname[L+1])$(M+1)", "w")
                for ie = 1:Dos_N
                    fill!(DSum, 0.0)
                    for ist = 1:Total_NumOrbs[atom]
                        if LM == Spe_Num_Relation[spe][ist]
                            for spin = 1:spinsize
                                DSum[spin] += Dos[spin][atom][ist][ie]
                            end
                        end
                    end

                    if SpinPol ∈ ("on", "nc")
                        @printf(PDos_data, "%lf  %lf  %lf  %lf  %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
                    else
                        @printf(PDos_data, "%lf  %lf  %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
                    end
                end
                close(PDos_data)
            end
        end
    end
end


function Calc_PDos_Orbital_proj(filename::String, material::CWF_model, Dos_Erange, DosE, Dos)
    
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    gsize = material.gsize
    DosEmin, DosEmax = Dos_Erange
    de = 0.01/eV2Hartree
    Dos_N = floor(Int, (DosEmax-DosEmin)/de)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)

    

    for ist = 1:gsize
        for spin = 1:spinsize, q = 1:Dos_N
            s1 = 0.0
            s2 = 0.0
            for ie = 2:2:q-1
                s1 += Dos[spin][ist][ie]
            end
                
            for ie = 3:2:q-1
                s2 += Dos[spin][ist][ie]
            end
            ssum[q,spin] = (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
        end


        println("Write $(filename).PDOS.Tetrahedron.orb$(ist)")
        PDos_data = open("$(filename).PDOS.Tetrahedron.orb$(ist)", "w")
        for ie = 1:Dos_N
            for spin = 1:spinsize
                DSum[spin] = Dos[spin][ist][ie]
            end

            if SpinPol ∈ ("on", "nc")
                @printf(PDos_data, "%lf  %lf  %lf  %lf  %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
            else
                @printf(PDos_data, "%lf  %lf  %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
            end
        end
        close(PDos_data)
    end
end