function Generate_HWF(material::LCPAO_model, OLPpos, Nk_hwf::Integer, pflag::Vector{Int64})    

    Natom = material.Natom
    Nspin = material.Nspin
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(fsize=="nc", 2*fsize, fsize)


    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    elseif SpinPol == "nc"
        Nocc = Int64(Valence_Electrons)
    end


    kpts_hwf = zeros(Float64, Nk_hwf+1)
    for ik = 0:Nk_hwf
        kpts_hwf[ik+1] = ik/Nk_hwf
    end


    tmpH = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    UVT = zeros(ComplexF64, Nocc, Nocc)
    UVT2 = zeros(ComplexF64, Nocc, Nocc)

    Sop = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Sop[spin] = zeros(ComplexF64, Nocc, Nocc)
    end
    

    Enk1 = Vector{Vector{ComplexF64}}(undef, spinsize)
    Enk2 = Vector{Vector{ComplexF64}}(undef, spinsize)
    Enk3 = Vector{Vector{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Enk1[spin] = zeros(Float64, Nfsize)
        Enk2[spin] = zeros(Float64, Nfsize)
        Enk3[spin] = zeros(Float64, Nfsize)
    end

    Cnk1 = Vector{Matrix{ComplexF64}}(undef, spinsize)
    Cnk2 = Vector{Matrix{ComplexF64}}(undef, spinsize)
    Cnk3 = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk1[spin] = zeros(ComplexF64, Nfsize, Nfsize)
        Cnk2[spin] = zeros(ComplexF64, Nfsize, Nfsize)
        Cnk3[spin] = zeros(ComplexF64, Nfsize, Nfsize)
    end




    zjk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        zjk[spin] = Vector{Vector{Float64}}(undef, Nk_hwf)
        for ik = 1:Nk_hwf
            zjk[spin][ik] = zeros(Float64, Nocc)
        end
    end


    Umjk = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Umjk[spin] = zeros(ComplexF64, Nocc, Nocc)
    end




    kpoints = Vector{Vector{Float64}}(undef, 2)
    for i = 1:2
        kpoints[i] = zeros(Float64, 3)
    end
    kloop = [[1,2,3], [2,3,1], [3,1,2]]
    k_axis = ["a-axis", "b-axis", "c-axis"]


    for k = 1:3
        if pflag[k] == 1
            for i1 = 0:Nk_hwf-1

                if i1 == 0
                    diag_flag = 0
                elseif i1 == Nk_hwf-1
                    diag_flag = 3
                else
                    diag_flag = 2
                end

                if mod(i1,2) == 0
                    if i1 ≠ Nk_hwf-1
                        Diag_k1k2!(
                            diag_flag, kpoints,
                            material,
                            tmpH, S, H, Enk1, Enk2, Cnk1, Cnk2)
                        Overlap_k1k2!(kpoints, material, Nocc, OLPpos, Cnk1, Cnk2, Sop)
                    else
                        Diag_k1k2!(
                            diag_flag, kpoints,
                            material,
                            tmpH, S, H, Enk1, Enk3, Cnk1, Cnk3)
                        Overlap_k1k2!(kpoints, material, Nocc, OLPpos, Cnk1, Cnk3, Sop)
                    end
                else
                    if i1 ≠ Nk_hwf-1
                        Diag_k1k2!(
                            diag_flag, kpoints,
                            material,
                            tmpH, S, H, Enk2, Enk1, Cnk2, Cnk1)
                        Overlap_k1k2!(kpoints, material, Nocc, OLPpos, Cnk2, Cnk1, Sop)
                    else
                        Diag_k1k2!(
                            diag_flag, kpoints,
                            material,
                            tmpH, S, H, Enk2, Enk3, Cnk2, Cnk3)
                        Overlap_k1k2!(kpoints, material, Nocc, OLPpos, Cnk2, Cnk3, Sop)
                    end
                end


                if i1 == 0
                    Enk3 = deepcopy(Enk1)
                    Cnk3 = deepcopy(Cnk1)
                end



                for spin = 1:spinsize
                    U, Σ, VT = svd(Sop[spin])
                    mul!(UVT, U, VT')
                    if i1 == 0
                        for ist = 1:Nocc, jst = 1:Nocc
                            Umjk[spin][ist,jst] = UVT[ist,jst]
                        end
                    else
                        mul!(UVT2, Umjk[spin], UVT)
                        for ist = 1:Nocc, jst = 1:Nocc
                            Umjk[spin][ist,jst] = UVT2[ist,jst]
                        end
                    end
                end
            end
        end
    end



    return zjk, Umjk
end


function Diag_k1k2!(
    diag_flag, kpoints::Vector{Vector{Float64}},
    material::LCPAO_model,
    tmpH, S, H, Enk1, Enk2, Cnk1, Cnk2)

    Natom = material.Natom
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks

    spinsize = ifelse(SpinPol=="on", 2, 1)


    if diag_flag == 0
        recalc = [true, true]
    elseif diag_flag == 1
        recalc = [true, false]
    elseif diag_flag == 2
        recalc = [false, true]
    elseif diag_flag == 3
        recalc = [false, false]
    end



    for ik = 1:2
        if recalc[ik]
            if SpinPol ∈ ("off", "on")
                for spin = 1:spinsize
                    HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpoints[ik])
                    if ik == 1
                        Enk1[spin], Cnk1[spin] = eigen(Hermitian(H), Hermitian(S))
                    elseif ik == 2
                        Enk2[spin], Cnk2[spin] = eigen(Hermitian(H), Hermitian(S))
                    end
                end   
            else
                HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpoints[ik])
                HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpoints[ik])
                @. S[1:fsize, 1:fsize] = tmpH
                @. S[fsize+1:end, fsize+1:end] = tmpH
                if ik == 1
                    Enk1[1], Cnk1[1] = eigen(Hermitian(H), Hermitian(S))
                elseif ik == 2
                    Enk2[1], Cnk2[1] = eigen(Hermitian(H), Hermitian(S))
                end    
            end
        end
    end
end


function Overlap_k1k2!(
    kpoints::Vector{Vector{Float64}},
    material::LCPAO_model,
    Nocc, OLPpos, Cnk1, Cnk2, Sop)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    OLP = material.OLP
    spinsize = ifelse(SpinPol=="on", 2, 1)


    if SpinPol ∈ ("off", "on")
        Overlap_k1k2_Col!(
            spinsize, Natom, FNAN, natn, ncn, atv_ijk, Gxyz, Total_NumOrbs, MP,
            Recvecs, Nocc, Cnk1, Cnk2, kpoints, OLP, OLPpos, Sop)
    else
        Overlap_k1k2_NonCol!(
            spinsize, Natom, FNAN, natn, ncn, atv_ijk, Gxyz, Total_NumOrbs, MP,
            Recvecs, Nocc, Cnk1, Cnk2, kpoints, OLP, OLPpos, Sop)
    end
end


function Overlap_k1k2_Col!(
    spinsize, Natom, FNAN, natn, ncn, atv_ijk, Gxyz, Total_NumOrbs, MP,
    Recvecs, Nocc, Cnk1, Cnk2, kpoints::Vector{Vector{Float64}}, OLP, OLPpos, Sop)

    dka = kpoints[2][1] - kpoints[1][1]
    dkb = kpoints[2][2] - kpoints[1][2]
    dkc = kpoints[2][3] - kpoints[1][3]

    dkx = dka*Recvecs[1,1] + dkb*Recvecs[2,1] + dkc*Recvecs[3,1]
    dky = dka*Recvecs[1,2] + dkb*Recvecs[2,2] + dkc*Recvecs[3,2]
    dkz = dka*Recvecs[1,3] + dkb*Recvecs[2,3] + dkc*Recvecs[3,3]

    for spin = 1:spinsize, μ = 1:Nocc, ν = 1:Nocc
        Sum_re = 0.0
        Sum_im = 0.0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            cell = ncn[atom][Rn]+1
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]
            Anum = MP[atom]
            Bnum = MP[jatom]
            l1, l2, l3 = atv_ijk[cell]

            kRn = 2*pi*(kpoints[2][1]*l1 + kpoints[2][2]*l2 + kpoints[2][3]*l3)
            kRn = kRn - dkx*Gxyz[atom][1] - dky*Gxyz[atom][2] - dkz*Gxyz[atom][3]
            si = sin(kRn)
            co = cos(kRn)

            for ist = 1:NO0, jst = 1:NO1
                tmp1r = real(Cnk1[spin][Anum+ist,μ])*real(Cnk2[spin][Bnum+jst,ν]) + imag(Cnk1[spin][Anum+ist,μ])*imag(Cnk1[spin][Bnum+jst,ν])
                tmp1i = real(Cnk1[spin][Anum+ist,μ])*imag(Cnk2[spin][Bnum+jst,ν]) - imag(Cnk1[spin][Anum+ist,μ])*real(Cnk1[spin][Bnum+jst,ν])

                tmp2r = co*tmp1r - si*tmp1i
                tmp2i = co*tmp1i + si*tmp1r

                tmp3r = OLP[atom][Rn][ist][jst]
                tmp3i = -dkx*OLPpos[1][atom][Rn][ist][jst] - dky*OLPpos[2][atom][Rn][ist][jst] - dkz*OLPpos[3][atom][Rn][ist][jst]

                Sum_re += tmp2r*tmp3r - tmp2i*tmp3i
                Sum_im += tmp2r*tmp3i + tmp2i*tmp3r
            end
        end

        Sop[spin][μ,ν] = Sum_re + im*Sum_im
    end
end


function Overlap_k1k2_NonCol!(
    spinsize, Natom, FNAN, natn, ncn, atv_ijk, Gxyz, Total_NumOrbs, MP,
    Recvecs, Nocc, Cnk1, Cnk2, kpoints::Vector{Vector{Float64}}, OLP, OLPpos, Sop)

    dka = kpoints[2][1] - kpoints[1][1]
    dkb = kpoints[2][2] - kpoints[1][2]
    dkc = kpoints[2][3] - kpoints[1][3]

    dkx = dka*Recvecs[1,1] + dkb*Recvecs[2,1] + dkc*Recvecs[3,1]
    dky = dka*Recvecs[1,2] + dkb*Recvecs[2,2] + dkc*Recvecs[3,2]
    dkz = dka*Recvecs[1,3] + dkb*Recvecs[2,3] + dkc*Recvecs[3,3]

    for spin = 1:spinsize, μ = 1:Nocc, ν = 1:Nocc
        Sum_re = 0.0
        Sum_im = 0.0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            cell = ncn[atom][Rn]+1
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]
            Anum = MP[atom]
            Bnum = MP[jatom]
            l1, l2, l3 = atv_ijk[cell]

            kRn = 2*pi*(kpoints[2][1]*l1 + kpoints[2][2]*l2 + kpoints[2][3]*l3)
            kRn = kRn - dkx*Gxyz[atom][1] - dky*Gxyz[atom][2] - dkz*Gxyz[atom][3]
            si = sin(kRn)
            co = cos(kRn)

            for ist = 1:NO0, jst = 1:NO1
                tmp1r  = real(Cnk1[spin][Anum+ist,μ])*real(Cnk2[spin][Bnum+jst,ν]) + imag(Cnk1[spin][Anum+ist,μ])*imag(Cnk1[spin][Bnum+jst,ν])
                tmp1r += real(Cnk1[spin][fsize+Anum+ist,μ])*real(Cnk2[spin][fsize+Bnum+jst,ν]) + imag(Cnk1[spin][fsize+Anum+ist,μ])*imag(Cnk1[spin][fsize+Bnum+jst,ν])
                tmp1r  = real(Cnk1[spin][Anum+ist,μ])*imag(Cnk2[spin][Bnum+jst,ν]) - imag(Cnk1[spin][Anum+ist,μ])*real(Cnk1[spin][Bnum+jst,ν])
                tmp1r += real(Cnk1[spin][fsize+Anum+ist,μ])*imag(Cnk2[spin][fsize+Bnum+jst,ν]) - imag(Cnk1[spin][fsize+Anum+ist,μ])*real(Cnk1[spin][fsize+Bnum+jst,ν])

                tmp2r = co*tmp1r - si*tmp1i
                tmp2i = co*tmp1i + si*tmp1r

                tmp3r = OLP[atom][Rn][ist][jst]
                tmp3i = -dkx*OLPpos[1][atom][Rn][ist][jst] - dky*OLPpos[2][atom][Rn][ist][jst] - dkz*OLPpos[3][atom][Rn][ist][jst]

                Sum_re += tmp2r*tmp3r - tmp2i*tmp3i
                Sum_im += tmp2r*tmp3i + tmp2i*tmp3r
            end
        end

        Sop[spin][μ,ν] = Sum_re + im*Sum_im
    end
end