#=
function MPI_Generate_HWF_slab(material::LCPAO_model, OLPpos, kpoints::KPoints, sflag::Integer)    

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Recvecs = material.Recvecs
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)

    kmesh = kpoints.kmesh
    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    NkAB = Int64(div(AllNkpt,kmesh[sflag]))

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    else
        Nocc = Int64(Valence_Electrons)
    end

    slab_c_len = sqrt(Recvecs[sflag,1]^2 + Recvecs[sflag,2]^2 + Recvecs[sflag,3]^2)

    if myrank == 0
        @show sflag
        @show kmesh
        @show AllNkpt
        @show NkAB
        @show slab_c_len
        @show Total_SpinS
        @show SpinPol
        @show Valence_Electrons
        @show Nocc
    end
    MPI.Barrier(comm)




    OLPexp = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Natom)
    for atom = 1:Natom
        OLPexp[atom] = Vector{Vector{Vector{ComplexF64}}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            OLPexp[atom][Rn] = Vector{Vector{ComplexF64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                OLPexp[atom][Rn][ist] = zeros(ComplexF64, Total_NumOrbs[natn[atom][Rn]])
            end
        end
    end

    if sflag == 1
        Set_OLPexp_HWF!(material, 1.0, 0.0, 0.0, OLPpos, OLPexp)
    elseif sflag == 2
        Set_OLPexp_HWF!(material, 0.0, 1.0, 0.0, OLPpos, OLPexp)
    else
        Set_OLPexp_HWF!(material, 0.0, 0.0, 1.0, OLPpos, OLPexp)
    end



    Mmnkb_HWF = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Mmnkb_HWF[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Mmnkb_HWF[spin][ik] = zeros(ComplexF64, BANDNUM, BANDNUM)
        end
    end



    tmpH = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    U = zeros(ComplexF64, Nocc, Nocc)
    VT = zeros(ComplexF64, Nocc, Nocc)
    UVT = zeros(ComplexF64, Nocc, Nocc)


    kpts1 = zeros(Float64, 3)
    kpts2 = zeros(Float64, 3)
    Lambda = zeros(ComplexF64, Nocc)
    Mmnkb = zeros(ComplexF64, Nocc, Nocc)
    Cnk2 = zeros(ComplexF64, Nfsize, Nfsize)


    # Enk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        # Enk[spin] = Vector{Vector{Float64}}(undef, AllNkpt)
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            # Enk[spin][ik] = zeros(Float64, Nfsize)
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end



    Wannier_Center = Vector{Vector{Float64}}(undef, spinsize)
    for spin = 1:spinsize
        Wannier_Center[spin] = zeros(Float64, Nocc)
    end


    Umjk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umjk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Umjk[spin][ik] = zeros(ComplexF64, Nocc, Nocc)
        end
    end


    for spin = 1:spinsize, ik = 1:MPI_Nkpt

        @. kpts1 = MPI_kpts[ik]
        @. kpts2 = MPI_kpts[ik]
        kpts2[sflag] = 1.0
        @show myrank, spin, ik, kpts1, kpts2
        Diag_k1k2!(material, spin, tmpH, S, H, Cnk[spin][ik], Cnk2, kpts1, kpts2)
        Overlap_k1k2!(material, Nocc, T, OLPexp, Cnk[spin][ik], Cnk2, Mmnkb, kpts1, kpts2)


        U, _, VT = svd(Mmnkb)
        mul!(UVT, U, VT')
        Lambda, Umjk[spin][ik] = eigen(UVT)

        
        for μ = 1:Nocc
            phase = -atan(imag(Lambda[μ])/real(Lambda[μ]))
            if phase < 0.0
                phase += 2*pi
            end
            Wannier_Center[spin][μ] += phase/slab_c_len/Ang_to_bohr/NkAB
        end
    end

    for spin = 1:spinsize
        MPI.Allreduce!(Wannier_Center[spin], MPI.SUM, comm)
    end
    MPI.Barrier(comm)

    if myrank == 0
        println("\n\nWannier_Center")
        for spin = 1:spinsize
            @show Wannier_Center[spin]
        end
    end
    MPI.Barrier(comm)


    return Cnk, Wannier_Center, Umjk
end
=#


@timeit timer "Generate_HWF_slab" function Generate_HWF_slab(material::LCPAO_model, OLPpos, kpoints::KPoints, sflag::Integer)    

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Recvecs = material.Recvecs
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)

    kmesh = kpoints.kmesh
    AllNkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    NkAB = Int64(div(AllNkpt,kmesh[sflag]))

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    elseif SpinPol == "nc"
        Nocc = Int64(Valence_Electrons)
    end

    slab_c_len = sqrt(Recvecs[sflag,1]^2 + Recvecs[sflag,2]^2 + Recvecs[sflag,3]^2)

    @show sflag
    @show kmesh
    @show AllNkpt
    @show NkAB
    @show slab_c_len
    @show Total_SpinS
    @show SpinPol
    @show Valence_Electrons
    @show Nocc



    OLPexp = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Natom)
    for atom = 1:Natom
        OLPexp[atom] = Vector{Vector{Vector{ComplexF64}}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            OLPexp[atom][Rn] = Vector{Vector{ComplexF64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                OLPexp[atom][Rn][ist] = zeros(ComplexF64, Total_NumOrbs[natn[atom][Rn]])
            end
        end
    end

    if sflag == 1
        Set_OLPexp_HWF!(material, 1.0/kmesh[1], 0.0, 0.0, OLPpos, OLPexp)
    elseif sflag == 2
        Set_OLPexp_HWF!(material, 0.0, 1.0/kmesh[2], 0.0, OLPpos, OLPexp)
    else
        Set_OLPexp_HWF!(material, 0.0, 0.0, 1.0/kmesh[3], OLPpos, OLPexp)
    end


    tmpH = zeros(ComplexF64, fsize, fsize)
    T = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    U = zeros(ComplexF64, Nocc, Nocc)
    V = zeros(ComplexF64, Nocc, Nocc)
    UVT = zeros(ComplexF64, Nocc, Nocc)


    kpts1 = zeros(Float64, 3)
    kpts2 = zeros(Float64, 3)
    Lambda = zeros(ComplexF64, Nocc)
    UVT_temp = zeros(ComplexF64, Nocc, Nocc)
    Mmnkb = zeros(ComplexF64, Nocc, Nocc)


    # Enk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        # Enk[spin] = Vector{Vector{Float64}}(undef, AllNkpt)
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, AllNkpt)
        for ik = 1:AllNkpt
            # Enk[spin][ik] = zeros(Float64, Nfsize)
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end


    WC_temp = zeros(Float64, Nocc)
    Wannier_Center = Vector{Vector{Float64}}(undef, spinsize)
    for spin = 1:spinsize
        Wannier_Center[spin] = zeros(Float64, Nocc)
    end


    Umjk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umjk[spin] = Vector{Matrix{ComplexF64}}(undef, AllNkpt)
        for ik = 1:AllNkpt
            Umjk[spin][ik] = zeros(ComplexF64, Nocc, Nocc)
        end
    end


    for spin = 1:spinsize, ik = 1:AllNkpt

        @. kpts1 = kpts[ik]
        @. kpts2 = kpts[ik]
        kpts2[sflag] += 1.0
        @show spin, ik, kpts1, kpts2
        Diaginalize!(material, spin, tmpH, S, H, Cnk[spin][ik], kpts1)
        Overlap_k1k2!(material, Nocc, T, OLPexp, Cnk[spin][ik], Mmnkb, kpts1, kpts2)

        U, _, V = svd(Mmnkb)
        mul!(UVT, U, V')
        Lambda, UVT_temp = eigen(UVT)

        for μ = 1:Nocc
            phase = -atan(imag(Lambda[μ]), real(Lambda[μ]))
            if phase < 0.0
                phase += 2*pi
            end
            WC_temp[μ] = phase/slab_c_len/Ang_to_bohr
        end

        WC_idx = sortperm(WC_temp)
        @. WC_temp = WC_temp[WC_idx]
        for μ = 1:Nocc
            Wannier_Center[spin][μ] += WC_temp[μ]
        end

        for ist = 1:Nocc, μ = 1:Nocc
            Umjk[spin][ik][ist,μ] = UVT_temp[ist,WC_idx[μ]]
        end
    end

    for spin = 1:spinsize, μ = 1:Nocc
        Wannier_Center[spin][μ] = Wannier_Center[spin][μ]/AllNkpt
    end

    println("\n\nWannier_Center")
    for spin = 1:spinsize
        # sort!(Wannier_Center[spin])
        @show Wannier_Center[spin]
    end


    return Cnk, Wannier_Center, Umjk
end


@timeit timer "Set_OLPexp_HWF!" function Set_OLPexp_HWF!(material::LCPAO_model, dka, dkb, dkc, OLPpos, OLPexp)
    
    Recvecs = material.Recvecs
    Natom = material.Natom
    Gxyz = material.Gxyz
    FNAN = material.FNAN
    natn = material.natn
    Total_NumOrbs = material.Total_NumOrbs
    OLP = material.OLP

    dkx = dka*Recvecs[1,1] + dkb*Recvecs[2,1] + dkc*Recvecs[3,1]
    dky = dka*Recvecs[1,2] + dkb*Recvecs[2,2] + dkc*Recvecs[3,2]
    dkz = dka*Recvecs[1,3] + dkb*Recvecs[2,3] + dkc*Recvecs[3,3]

    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        kRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        ex = cis(-kRn)
        _OLP = OLP[atom][Rn]
        _OLPpos1 = OLPpos[1][atom][Rn]
        _OLPpos2 = OLPpos[2][atom][Rn]
        _OLPpos3 = OLPpos[3][atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            tmp = _OLP[ist][jst] - im*(dkx*_OLPpos1[ist][jst] + dky*_OLPpos2[ist][jst] + dkz*_OLPpos3[ist][jst])
            OLPexp[atom][Rn][ist][jst] = tmp*ex
        end
    end
end


@timeit timer "Diaginalize!" function Diaginalize!(material::LCPAO_model, spin, tmpH, S, H, Cnk, kpts::Vector{Float64})

    SpinPol = material.SpinPol
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks

    if SpinPol ∈ ("off", "on")
        HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
        Cnk[:,:] = eigvecs(Hermitian(H), Hermitian(S))
    else
        HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
        @. S[1:fsize, 1:fsize] = tmpH
        @. S[fsize+1:end, fsize+1:end] = tmpH
        Cnk[:,:] = eigvecs(Hermitian(H), Hermitian(S))
    end
end


function Set_Tmatrix!(
    T, 
    OLPexp::Vector{Vector{Vector{Vector{ComplexF64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts1::Vector{Float64}, kpts2::Vector{Float64})
    
    dka = kpts2[1] - kpts1[1]
    dkb = kpts2[2] - kpts1[2]
    dkc = kpts2[3] - kpts1[3]

    ka, kb, kc = kpts2
    fill!(T, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        l1, l2, l3 = atv_ijk[cell]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        kRn = ka*l1 + kb*l2 + kc*l3
        kRn2 = kRn - dka*l1 - dkb*l2 - dkc*l3
        Anum = MP[atom]
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        ex2 = cispi(-2*kRn2)
        # ex2 = cispi(2*kRn2)
        _OLPexp = OLPexp[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            tmp = _OLPexp[ist][jst]
            T[Anum+ist,Bnum+jst] += 0.5*tmp*ex
            T[Bnum+jst,Anum+ist] += 0.5*tmp*ex2
        end
    end
end


@timeit timer "Overlap_k1k2" function Overlap_k1k2!(material::LCPAO_model, Nocc, T, OLPexp, Cnk, Mmnkb, kpts1, kpts2)

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

    if SpinPol ∈ ("off", "on")
        Overlap_k1k2_Col!(
            Natom, FNAN, natn, ncn, atv_ijk, Gxyz, Total_NumOrbs, MP,
            Recvecs, Nocc, Cnk, OLP, OLPpos, Mmnkb, kpts1, kpts2)
    else
        Overlap_k1k2_NonCol!(
            Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP,
            Nocc, Cnk, T, OLPexp, Mmnkb, kpts1, kpts2)
    end
end


function Overlap_k1k2_NonCol!(
    Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP,
    Nocc, Cnk, T, OLPexp, Mmnkb, kpts1, kpts2)

    fsize = sum(Total_NumOrbs)
    fill!(Mmnkb, 0.0)
    Set_Tmatrix!(T, OLPexp, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts1, kpts2)

    for μ = 1:Nocc, ν = 1:Nocc
        Sum = 0.0 + im*0.0
        for atom = 1:Natom, jatom = 1:Natom
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]
            Anum = MP[atom]
            Bnum = MP[jatom]
            @inbounds for ist = 1:NO0, jst = 1:NO1
                Sum += conj(Cnk[Anum+ist,μ])*Cnk[Bnum+jst,ν]*T[Anum+ist,Bnum+jst]
                Sum += conj(Cnk[fsize+Anum+ist,μ])*Cnk[fsize+Bnum+jst,ν]*T[Anum+ist,Bnum+jst]
            end
        end

        Mmnkb[μ,ν] = Sum
    end
end
