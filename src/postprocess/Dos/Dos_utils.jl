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


function Get_EVec_Collinear!(material::LCPAO_model, Nkpt, kpts, Cnk, EVec, iemin, iemax)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol == "off", 1, 2)
    OLP = material.OLP
    ChemP = material.ChemP


    Cnk_tmp = zeros(ComplexF64, fsize, fsize)
    SD = zeros(Float32, fsize)
    S = zeros(ComplexF64, fsize, fsize)

    for spin = 1:spinsize, ik = 1:Nkpt
        @. Cnk_tmp = Cnk[spin][ik]
        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik,:])
        for l = iemin:iemax
            
            fill!(SD, 0.0)
            for atom = 1:Natom, jatom = 1:Natom
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                Anum = MP[atom]
                Bnum = MP[jatom]

                for ist = 1:NO0, jst = 1:NO1
                    tmp = conj(Cnk_tmp[Anum+ist,l]) * Cnk_tmp[Bnum+jst,l]
                    SD[Anum+ist] += real(tmp * S[Anum+ist,Bnum+jst])
                end
            end

            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                Anum = MP[atom]
                EVec[spin][ik][atom][ist][l-iemin+1] = SD[Anum+ist]
            end
        end
    end
end


function Get_EVec_NonCollinear!(material::LCPAO_model, Nkpt, kpts, Cnk, EVec, iemin, iemax)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    OLP = material.OLP
    ChemP = material.ChemP
    Angle_spin = Get_Angle_spin(material)


    Cnk_tmp = zeros(ComplexF64, 2*fsize, 2*fsize)
    SD = zeros(Float32, 2*fsize)
    S = zeros(ComplexF64, fsize, fsize)

    for ik = 1:Nkpt
        @. Cnk_tmp = Cnk[1][ik]
        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik,:])
        for l = iemin:iemax
            
            fill!(SD, 0.0)
            for atom = 1:Natom, jatom = 1:Natom
                theta = Angle_spin[atom][1]
                phi = Angle_spin[atom][2]
                sit = sin(theta)
                cot = cos(theta)
                sip = sin(phi)
                cop = cos(phi)
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                Anum = MP[atom]
                Bnum = MP[jatom]

                for ist = 1:NO0, jst = 1:NO1
                    tmp_uu = real(conj(Cnk_tmp[Anum+ist,l]) * Cnk_tmp[Bnum+jst,l] * S[Anum+ist,Bnum+jst])
                    tmp_dd = real(conj(Cnk_tmp[Anum+ist+fsize,l]) * Cnk_tmp[Bnum+jst+fsize,l] * S[Anum+ist,Bnum+jst])
                    tmp_ud_real = real(conj(Cnk_tmp[Anum+ist,l]) * Cnk_tmp[Bnum+jst+fsize,l] * S[Anum+ist,Bnum+jst])
                    tmp_ud_imag = imag(conj(Cnk_tmp[Anum+ist,l]) * Cnk_tmp[Bnum+jst+fsize,l] * S[Anum+ist,Bnum+jst])
                    SD[Anum+ist] += 0.5*(tmp_uu + tmp_dd) + 0.5*cot*(tmp_uu - tmp_dd) + (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
                    SD[Anum+ist+fsize] += 0.5*(tmp_uu + tmp_dd) - 0.5*cot*(tmp_uu - tmp_dd) - (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
                end
            end

            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                Anum = MP[atom]
                EVec[1][ik][atom][ist][l-iemin+1] = SD[Anum+ist]
                EVec[2][ik][atom][ist][l-iemin+1] = SD[Anum+ist+fsize]
            end
        end
    end
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
