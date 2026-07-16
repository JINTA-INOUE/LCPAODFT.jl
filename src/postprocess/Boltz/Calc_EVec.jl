@timeit timer "Calc_EVec" function Calc_EVec(boltz_setup::Boltz_Setup, Cnk)

    material = boltz_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    kmesh = boltz_setup.kmesh
    Nkpt = prod(kmesh)
    Nstate = Set_Nstate(material)
    mat_type = boltz_setup.mat_type
    plane_type = boltz_setup.plane_type

    if plane_type
        if kmesh[3] ≠ 1
            error("please check plane_type, kmesh.")
        end
    end


    EVec = Vector{Vector{Vector{Vector{Float32}}}}(undef, spinsize)
    for spin = 1:spinsize
        EVec[spin] = Vector{Vector{Vector{Float32}}}(undef, Nkpt)
        for ik = 1:Nkpt
            EVec[spin][ik] = Vector{Vector{Float32}}(undef, Nstate)
            for ist = 1:Nstate
                EVec[spin][ik][ist] = zeros(Float32, Nstate)
            end
        end
    end


    if mat_type == "LCPAO"
        # return Calc_EVec_LCPAO!(boltz_setup, Cnk, EVec)
    elseif mat_type == "CWF"
        return Calc_EVec_Wannier!(boltz_setup, Cnk, EVec)
    elseif mat_type == "HWF"
        # return Calc_EVec_HWFs!(boltz_setup)
    else
        error("please check WF_mode")
    end
end


function Calc_EVec_LCPAO!(SpinPol, material::LCPAO_model, Cnk, EVec)

    if SpinPol ∈ ("off", "on")
        Calc_EVec_LCPAO_Col!(material, Nkpt, kpts, Cnk, EVec)
    else
        Calc_EVec_LCPAO_NonCol!(material, Nkpt, kpts, Cnk, EVec)
    end
end


function Calc_EVec_LCPAO_Col!(material::LCPAO_model, Nkpt, kpts, Cnk, EVec)

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


function Calc_EVec_LCPAO_NonCol!(material::LCPAO_model, Nkpt, kpts, Cnk, EVec)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize
    OLP = material.OLP
    ChemP = material.ChemP
    Angle_spin = Get_Angle_spin(material)


    Cnk_tmp = zeros(ComplexF64, 2*fsize, 2*fsize)
    SD = zeros(Float32, 2*fsize)
    S = zeros(ComplexF64, fsize, fsize)

    for ik = 1:Nkpt
        @. Cnk_tmp = Cnk[1][ik]
        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
        for μ = 1:Nfsize
            fill!(SD, 0.0)
            for atom = 1:Natom, jatom = 1:Natom
                theta, phi = Angle_spin[atom]
                sit = sin(theta)
                cot = cos(theta)
                sip = sin(phi)
                cop = cos(phi)
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                Anum = MP[atom]
                Bnum = MP[jatom]

                @inbounds for ist = 1:NO0, jst = 1:NO1
                    tmp_uu = real(conj(Cnk_tmp[Anum+ist,μ])*Cnk_tmp[Bnum+jst,μ]*S[Anum+ist,Bnum+jst])
                    tmp_dd = real(conj(Cnk_tmp[Anum+ist+fsize,μ])*Cnk_tmp[Bnum+jst+fsize,μ]*S[Anum+ist,Bnum+jst])
                    tmp_ud_real = real(conj(Cnk_tmp[Anum+ist,μ])*Cnk_tmp[Bnum+jst+fsize,μ]*S[Anum+ist,Bnum+jst])
                    tmp_ud_imag = imag(conj(Cnk_tmp[Anum+ist,μ])*Cnk_tmp[Bnum+jst+fsize,μ]*S[Anum+ist,Bnum+jst])
                    SD[Anum+ist] += 0.5*(tmp_uu + tmp_dd) + 0.5*cot*(tmp_uu - tmp_dd) + (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
                    SD[Anum+ist+fsize] += 0.5*(tmp_uu + tmp_dd) - 0.5*cot*(tmp_uu - tmp_dd) - (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
                end
            end


            jst = 0
            @inbounds for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                jst += 1
                Anum = MP[atom]
                EVec[1][ik][jst][μ] = SD[Anum+ist]
                EVec[1][ik][jst+fsize][μ] = SD[Anum+ist+fsize]
            end
        end
    end
end


function Calc_EVec_Wannier!(boltz_setup::Boltz_Setup, Cnk, EVec)

    material = boltz_setup.material
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    kmesh = boltz_setup.kmesh
    Nkpt = prod(kmesh)

    for spin = 1:spinsize, ik = 1:Nkpt
        _Cnk = Cnk[spin][ik]
        @inbounds for ist = 1:Nwann, μ = 1:Nwann
            tmp = conj(_Cnk[ist,μ]) * _Cnk[ist,μ]
            EVec[spin][ik][ist][μ] = real(tmp)
        end
    end
end


function Calc_EVec_HWFs!(boltz_setup::Boltz_Setup, Cnk, EVec)

    material = boltz_setup.material
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    kmesh = boltz_setup.kmesh
    Nkpt = prod(kmesh)


    for spin = 1:spinsize, ik = 1:Nkpt
        _Cnk = Cnk[spin][ik]
        @inbounds for ist = 1:Nwann, μ = 1:Nwann
            tmp = conj(_Cnk[ist,μ]) * _Cnk[ist,μ]
            EVec[spin][ik][ist][μ] = real(tmp)
        end
    end
    
    return EVec
end