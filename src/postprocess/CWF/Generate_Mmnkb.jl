@timeit timer "Generate_Mmnkb" function Generate_Mmnkb(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, MinN, MaxN)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    spinsize = cwf_setup.spinsize
    filename = cwf_setup.filename
    material = cwf_setup.material
    SpinPol = material.SpinPol
    mlwf_kpoints = cwf_setup.mlwf_kpoints
    MPI_Nkpt = mlwf_kpoints.MPI_Nkpt
    tot_bvector = mlwf_kpoints.tot_bvector
    bvector = mlwf_kpoints.bvector
    BANDNUM = MaxN - MinN + 1

    myrank == 0 && println("<Set_OLPexp>")
    OLPexp = Set_OLPexp(material, tot_bvector, bvector)



    Mmnkb = Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Mmnkb[spin] = Vector{Vector{Matrix{ComplexF64}}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Mmnkb[spin][ik] = Vector{Matrix{ComplexF64}}(undef, tot_bvector)
            for ib = 1:tot_bvector
                Mmnkb[spin][ik][ib] = zeros(ComplexF64, BANDNUM, BANDNUM)
            end
        end
    end


    if SpinPol ∈ ("off", "on")
        Generate_Mmnkb_Col!(filename, MinN, MaxN, OLPexp, Mmnkb, material, mlwf_kpoints)
    elseif SpinPol == "nc"
        Generate_Mmnkb_NonCol!(filename, MinN, MaxN, OLPexp, Mmnkb, material, mlwf_kpoints)
    else
        error("please check SpinPol.")
    end


    Write_Variable_work_file(filename, myrank, Mmnkb, "Mmnkb")
    if myrank == 0
        if SpinPol ∈ ("off", "nc")
            Generate_Mmnkb_Spindeg1(filename, BANDNUM, mlwf_kpoints)
        elseif SpinPol == "on"
            Generate_Mmnkb_Spindeg2(filename, BANDNUM, mlwf_kpoints)
        else
            error("please check SpinPol.")
        end
    end
    MPI.Barrier(comm)
end


function Generate_Mmnkb_Spindeg1(filename::String, BANDNUM, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    kmesh = mlwf_kpoints.kmesh
    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    frac_bv_int = mlwf_kpoints.frac_bv_int
    work_dirname = pwd()*"/"*filename*"_work_cwf"
    
    data_Mmnkb = open(filename*".mmn", "w")
    @printf(data_Mmnkb, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
    @printf(data_Mmnkb, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Mmnkb$rank.jld2"
        data = jldopen(work_file, "r")
        Mmnkb = data["Mmnkb"]
        _Write_Mmnkb(data_Mmnkb, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[1])
        close(data)
    end
    close(data_Mmnkb)
end


function Generate_Mmnkb_Spindeg2(filename::String, BANDNUM, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    kmesh = mlwf_kpoints.kmesh
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    frac_bv_int = mlwf_kpoints.frac_bv_int
    work_dirname = pwd()*"/"*filename*"_work_cwf"
    
    data_Mmnkb1 = open(filename*"_1.mmn", "w")
    data_Mmnkb2 = open(filename*"_2.mmn", "w")
    @printf(data_Mmnkb1, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
    @printf(data_Mmnkb1, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)
    @printf(data_Mmnkb2, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
    @printf(data_Mmnkb2, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Mmnkb$rank.jld2"
        data = jldopen(work_file, "r")
        Mmnkb = data["Mmnkb"]
        _Write_Mmnkb(data_Mmnkb1, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[1])
        _Write_Mmnkb(data_Mmnkb2, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[2])
        close(data)
    end
    close(data_Mmnkb1)
    close(data_Mmnkb2)
end


function _Write_Mmnkb(data, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb)
    
    kmesh1, kmesh2, kmesh3 = kmesh

    for ik = 1:MPI_Nkpt, ib = 1:tot_bvector
        ik2 = ik + knum
        i = div(ik2-1, kmesh2*kmesh3)
        j = div(ik2-1 - i*kmesh2*kmesh3, kmesh3)
        k = ik2 - i*kmesh2*kmesh3 - j*kmesh3 - 1
        ii = ((i+frac_bv_int[ib][1])%kmesh1 + kmesh1)%kmesh1
        jj = ((j+frac_bv_int[ib][2])%kmesh2 + kmesh2)%kmesh2
        kk = ((k+frac_bv_int[ib][3])%kmesh3 + kmesh3)%kmesh3

        ikb = ii*kmesh2*kmesh3 + jj*kmesh3 + kk

        g1 = div((i+frac_bv_int[ib][1])-ii, kmesh1)
        g2 = div((j+frac_bv_int[ib][2])-jj, kmesh2)
        g3 = div((k+frac_bv_int[ib][3])-kk, kmesh3)
                
        @printf(data, "%d %d  %d %d %d\n", ik2, ikb+1, g1, g2, g3)

        @inbounds for ν = 1:BANDNUM, μ = 1:BANDNUM
            @printf(data, " %22.14f %22.14f\n", real(Mmnkb[ik][ib][μ,ν]), imag(Mmnkb[ik][ib][μ,ν]))
        end
    end
end


function Generate_Mmnkb_Col!(filename::String, MinN, MaxN, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Recvecs = material.Recvecs
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    BANDNUM = MaxN - MinN + 1

    AllNkpt = mlwf_kpoints.AllNkpt
    MPI_Nkpt = mlwf_kpoints.MPI_Nkpt
    MPI_kpts = mlwf_kpoints.MPI_kpts
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    kplusb = mlwf_kpoints.kplusb
    frac_bv = mlwf_kpoints.frac_bv
    knum = MPkpts[myrank+1]

    myrank == 0 && println("Calculating Mmnk...")


    k1 = zeros(Float64, 3)
    k2 = zeros(Float64, 3)
    dk = zeros(Float64, 3)
    Cnk1 = zeros(ComplexF64, fsize, fsize)
    Cnk2 = zeros(ComplexF64, fsize, fsize)


    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        ik2 = ik + knum
        @printf("\tmyrank = %3d %d/%d\n", myrank, ik2, AllNkpt)
        for ib = 1:tot_bvector
            kk = kplusb[ik2][ib]+1
            @. k1 = MPI_kpts[ik]
            @. k2 = MPI_kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
            
            Read_Cnk_work!(filename, ik2, Cnk1)
            Read_Cnk_work!(filename, kk, Cnk2)
            Calc_Mmnkb_Col!(Mmnkb[spin][ik][ib], MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], Cnk1, Cnk2)
        end
    end
    MPI.Barrier(comm)
end


function Generate_Mmnkb_NonCol!(filename::String, MinN, MaxN, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Recvecs = material.Recvecs
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    atv_ijk = material.atv_ijk
    Nfsize = 2*sum(Total_NumOrbs)
    BANDNUM = MaxN - MinN + 1

    AllNkpt = mlwf_kpoints.AllNkpt
    MPI_Nkpt = mlwf_kpoints.MPI_Nkpt
    MPI_kpts = mlwf_kpoints.MPI_kpts
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    kplusb = mlwf_kpoints.kplusb
    frac_bv = mlwf_kpoints.frac_bv
    knum = MPkpts[myrank+1]


    myrank == 0 && println("Calculating Mmnk...")

    k1 = zeros(Float64, 3)
    k2 = zeros(Float64, 3)
    dk = zeros(Float64, 3)
    Cnk1 = zeros(ComplexF64, Nfsize, Nfsize)
    Cnk2 = zeros(ComplexF64, Nfsize, Nfsize)

    for ik = 1:MPI_Nkpt
        ik2 = ik + knum
        @printf("\tmyrank = %3d %d/%d\n", myrank, ik2, AllNkpt)
        for ib = 1:tot_bvector
            kk = kplusb[ik2][ib]+1
            @. k1 = MPI_kpts[ik]
            @. k2 = MPI_kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
            
            Read_Cnk_work!(filename, ik2, Cnk1)
            Read_Cnk_work!(filename, kk, Cnk2)
            Calc_Mmnkb_NonCol!(Mmnkb[1][ik][ib], MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], Cnk1, Cnk2)
        end
    end
    MPI.Barrier(comm)
end


@timeit timer "Calc_Mmnkb" function Calc_Mmnkb_Col!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, Cnk1, Cnk2)
    
    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]

    band_offset = MinN - 1
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        cell  = ncn[atom][Rn]+1
        l1, l2, l3 = atv_ijk[cell]
        NO0  = Total_NumOrbs[atom]
        NO1  = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        kRn1 = -2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn
        kRn2 =  2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn
        ex1 = cis(kRn2)
        ex2 = cis(kRn1)

        _OLPexp = OLPexp[atom][Rn]

        for jst = 1:NO1, ist = 1:NO0
            aidx = Anum + ist
            bidx = Bnum + jst
            sij = _OLPexp[ist][jst]
            α1 = 0.5*ex1*sij
            α2 = 0.5*ex2*sij
            
            for ν = 1:BANDNUM
                νidx = ν + band_offset
                c2b = Cnk2[bidx,νidx]
                c2a = Cnk2[aidx,νidx]

                @inbounds for μ = 1:BANDNUM
                    μidx = μ + band_offset
                    c1a = Cnk1[aidx,μidx]
                    c1b = Cnk1[bidx,μidx]

                    Mmnkb[μ,ν] += α1*conj(c1a)*c2b + α2*conj(c1b)*c2a
                end
            end
        end
    end
end


@timeit timer "Calc_Mmnkb" function Calc_Mmnkb_NonCol!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, Cnk1, Cnk2)
    
    fsize = sum(Total_NumOrbs)
    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]

    band_offset = MinN - 1    
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        cell  = ncn[atom][Rn]+1
        l1, l2, l3 = atv_ijk[cell]
        NO0  = Total_NumOrbs[atom]
        NO1  = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        kRn1 = -2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn
        kRn2 =  2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn
        ex1 = cis(kRn2)
        ex2 = cis(kRn1)

        _OLPexp = OLPexp[atom][Rn]

        for jst = 1:NO1, ist = 1:NO0
            aidx1 = Anum + ist
            aidx2 = Anum + ist + fsize
            bidx1 = Bnum + jst
            bidx2 = Bnum + jst + fsize

            sij = _OLPexp[ist][jst]
            α1 = 0.5*ex1*sij
            α2 = 0.5*ex2*sij
            
            for ν = 1:BANDNUM
                νidx = ν + band_offset
                c2b1 = Cnk2[bidx1,νidx]
                c2b2 = Cnk2[bidx2,νidx]
                c2a1 = Cnk2[aidx1,νidx]
                c2a2 = Cnk2[aidx2,νidx]
                
                @inbounds for μ = 1:BANDNUM
                    μidx = μ + band_offset
                    c1a1 = Cnk1[aidx1,μidx]
                    c1a2 = Cnk1[aidx2,μidx]
                    c1b1 = Cnk1[bidx1,μidx]
                    c1b2 = Cnk1[bidx2,μidx]

                    Mmnkb[μ,ν] += α1*conj(c1a1)*c2b1 + α2*conj(c1b1)*c2a1
                    Mmnkb[μ,ν] += α1*conj(c1a2)*c2b2 + α2*conj(c1b2)*c2a2
                end
            end
        end
    end
end