function Calc_EVec_Dos(material::LCPAO_model, kpoints::KPoints, Cnk, iemin, iemax)

    Natom = material.Natom
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Total_NumOrbs = material.Total_NumOrbs
    Nkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    neg = iemax-iemin+1

    EVec = Vector{Vector{Vector{Vector{Vector{Float32}}}}}(undef, spinsize)
    for spin = 1:spinsize
        EVec[spin] = Vector{Vector{Vector{Vector{Float32}}}}(undef, Nkpt)
        for ik = 1:Nkpt
            EVec[spin][ik] = Vector{Vector{Vector{Float32}}}(undef, Natom)
            for atom = 1:Natom
                EVec[spin][ik][atom] = Vector{Vector{Float32}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    EVec[spin][ik][atom][ist] = zeros(Float32, neg)
                end
            end
        end
    end

    if SpinPol ∈ ("off", "on")
        Get_EVec_Collinear!(material, Nkpt, kpts, Cnk, EVec, iemin, iemax)
    else
        Get_EVec_NonCollinear!(material, Nkpt, kpts, Cnk, EVec, iemin, iemax)
    end


    return EVec
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
        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
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
        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
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
                    tmp_uu = real(conj(Cnk_tmp[Anum+ist,l])*Cnk_tmp[Bnum+jst,l]*S[Anum+ist,Bnum+jst])
                    tmp_dd = real(conj(Cnk_tmp[Anum+ist+fsize,l])*Cnk_tmp[Bnum+jst+fsize,l]*S[Anum+ist,Bnum+jst])
                    tmp_ud_real = real(conj(Cnk_tmp[Anum+ist,l])*Cnk_tmp[Bnum+jst+fsize,l]*S[Anum+ist,Bnum+jst])
                    tmp_ud_imag = imag(conj(Cnk_tmp[Anum+ist,l])*Cnk_tmp[Bnum+jst+fsize,l]*S[Anum+ist,Bnum+jst])
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


function Calc_EVec_Dos(material::CWF_model, kpoints::KPoints, Cnk, iemin, iemax)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    gsize = material.gsize
    Nwann = material.Ngsize
    Nkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    neg = iemax-iemin+1

    EVec = Vector{Vector{Vector{Vector{Float32}}}}(undef, spinsize)
    for spin = 1:spinsize
        EVec[spin] = Vector{Vector{Vector{Float32}}}(undef, Nkpt)
        for ik = 1:Nkpt
            EVec[spin][ik] = Vector{Vector{Float32}}(undef, gsize)
            for ist = 1:gsize
                EVec[spin][ik][ist] = zeros(Float32, neg)
            end
        end
    end

    if SpinPol ∈ ("off", "on")
        @inbounds for spin = 1:spinsize, ik = 1:Nkpt, ist = 1:gsize, μ = iemin:iemax
            tmp = conj(Cnk[spin][ik][ist,μ]) * Cnk[spin][ik][ist,μ]
            EVec[spin][ik][ist][μ-iemin+1] = real(tmp)
        end
    else
        @inbounds for ik = 1:Nkpt, ist = 1:gsize, μ = iemin:iemax
            tmp1 = conj(Cnk[1][ik][ist,μ]) * Cnk[1][ik][ist,μ]
            tmp2 = conj(Cnk[1][ik][ist+gsize,μ]) * Cnk[1][ik][ist+gsize,μ]
            EVec[1][ik][ist][μ-iemin+1] = real(tmp1)
            EVec[2][ik][ist][μ-iemin+1] = real(tmp2)
        end
    end


    return EVec
end
