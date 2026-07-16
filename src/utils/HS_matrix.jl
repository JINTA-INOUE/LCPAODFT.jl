@timeit timer "HS_matrix!" function HS_matrix!(
    S, 
    OLP::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fill!(S, 0.0)
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        _OLP = OLP[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            S[Anum+ist,Bnum+jst] += _OLP[ist][jst]*ex
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix!(
    S, H, 
    OLP::Vector{Vector{Vector{Vector{Float64}}}}, 
    Hks::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fill!(S, 0.0)
    fill!(H, 0.0)
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        _OLP = OLP[atom][Rn]
        _Hks = Hks[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            S[Anum+ist,Bnum+jst] += _OLP[ist][jst]*ex
            H[Anum+ist,Bnum+jst] += _Hks[ist][jst]*ex
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix!(
    S, 
    OLP::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn)
    
    fill!(S, 0.0)
    hst = 0
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        NO1 = Total_NumOrbs[jatom]
        Bnum = MP[jatom]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            S[Anum+ist,Bnum+jst] += OLP[hst]
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix!(
    S, 
    OLP::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fill!(S, 0.0)
    hst = 0
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            S[Anum+ist,Bnum+jst] += OLP[hst]*ex
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix!(
    S, H, 
    OLP::Vector{Float64}, 
    Hks::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn)
    
    fill!(S, 0.0)
    fill!(H, 0.0)
    hst = 0
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        NO1 = Total_NumOrbs[jatom]
        Bnum = MP[jatom]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            S[Anum+ist,Bnum+jst] += OLP[hst]
            H[Anum+ist,Bnum+jst] += Hks[hst]
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix!(
    S, H, 
    OLP::Vector{Float64}, 
    Hks::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fill!(S, 0.0)
    fill!(H, 0.0)
    hst = 0
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            S[Anum+ist,Bnum+jst] += OLP[hst]*ex
            H[Anum+ist,Bnum+jst] += Hks[hst]*ex
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix_NC!(
    H, 
    Hks::Vector{Vector{Float64}},
    iHks::Vector{Vector{Vector{Vector{Vector{Float64}}}}},
    Natom, Total_NumOrbs, MP, FNAN, natn)
    
    fsize = sum(Total_NumOrbs)

    _Hks_uu = Hks[1]
    _Hks_dd = Hks[2]
    _Hks_ud = Hks[3]
    _iHks_ud1 = Hks[4]
    iHks1 = iHks[1]
    iHks2 = iHks[2]
    iHks3 = iHks[3]

    hst = 0
    fill!(H, 0.0)
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        NO1 = Total_NumOrbs[jatom]
        Bnum = MP[jatom]
        _iHks_uu = iHks1[atom][Rn]
        _iHks_dd = iHks2[atom][Rn]
        _iHks_ud2 = iHks3[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            Hks_uu = _Hks_uu[hst]
            Hks_dd = _Hks_dd[hst]
            Hks_ud = _Hks_ud[hst]
            iHks_ud1 = _iHks_ud1[ist][jst]
            iHks_uu = _iHks_uu[ist][jst]
            iHks_dd = _iHks_dd[ist][jst]
            iHks_ud2 = _iHks_ud2[ist][jst]
            H[Anum+ist,Bnum+jst] += Hks_uu + im*iHks_uu
            H[fsize+Anum+ist,fsize+Bnum+jst] += Hks_dd + im*iHks_dd
            H[Anum+ist,fsize+Bnum+jst] += Hks_ud + im*(iHks_ud1 + iHks_ud2)
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix_NC!(
    H, 
    Hks::Vector{Vector{Float64}},
    iHks::Vector{Vector{Vector{Vector{Vector{Float64}}}}},
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fsize = sum(Total_NumOrbs)

    _Hks_uu = Hks[1]
    _Hks_dd = Hks[2]
    _Hks_ud = Hks[3]
    _iHks_ud1 = Hks[4]
    iHks1 = iHks[1]
    iHks2 = iHks[2]
    iHks3 = iHks[3]

    hst = 0
    fill!(H, 0.0)
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        _iHks_uu = iHks1[atom][Rn]
        _iHks_dd = iHks2[atom][Rn]
        _iHks_ud2 = iHks3[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            Hks_uu = _Hks_uu[hst]
            Hks_dd = _Hks_dd[hst]
            Hks_ud = _Hks_ud[hst]
            iHks_ud1 = _iHks_ud1[hst]
            iHks_uu = _iHks_uu[ist][jst]
            iHks_dd = _iHks_dd[ist][jst]
            iHks_ud2 = _iHks_ud2[ist][jst]
            H[Anum+ist,Bnum+jst] += (Hks_uu + im*iHks_uu)*ex
            H[fsize+Anum+ist,fsize+Bnum+jst] += (Hks_dd + im*iHks_dd)*ex
            H[Anum+ist,fsize+Bnum+jst] += (Hks_ud + im*(iHks_ud1 + iHks_ud2))*ex
        end
    end
end


@timeit timer "HS_matrix!" function HS_matrix_NC!(
    H, 
    Hks::Vector{Vector{Vector{Vector{Vector{Float64}}}}},
    iHks::Vector{Vector{Vector{Vector{Vector{Float64}}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    ka, kb, kc = kpts
    fsize = sum(Total_NumOrbs)

    Hks1 = Hks[1]
    Hks2 = Hks[2]
    Hks3 = Hks[3]
    Hks4 = Hks[4]
    iHks1 = iHks[1]
    iHks2 = iHks[2]
    iHks3 = iHks[3]

    fill!(H, 0.0)
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        _Hks_uu = Hks1[atom][Rn]
        _Hks_dd = Hks2[atom][Rn]
        _Hks_ud = Hks3[atom][Rn]
        _iHks_ud1 = Hks4[atom][Rn]
        _iHks_uu = iHks1[atom][Rn]
        _iHks_dd = iHks2[atom][Rn]
        _iHks_ud2 = iHks3[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            Hks_uu = _Hks_uu[ist][jst]
            Hks_dd = _Hks_dd[ist][jst]
            Hks_ud = _Hks_ud[ist][jst]
            iHks_ud1 = _iHks_ud1[ist][jst]
            iHks_uu = _iHks_uu[ist][jst]
            iHks_dd = _iHks_dd[ist][jst]
            iHks_ud2 = _iHks_ud2[ist][jst]
            H[Anum+ist,Bnum+jst]             += (Hks_uu + im*iHks_uu)*ex
            H[fsize+Anum+ist,fsize+Bnum+jst] += (Hks_dd + im*iHks_dd)*ex
            H[Anum+ist,fsize+Bnum+jst]       += (Hks_ud + im*(iHks_ud1 + iHks_ud2))*ex
        end
    end
end