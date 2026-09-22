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
    MPI.Barrier(comm)
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


function Generate_Mmnkb_Spindeg1(filename::AbstractString, BANDNUM, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    kmesh = mlwf_kpoints.kmesh
    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    frac_bv_int = mlwf_kpoints.frac_bv_int
    open(_cwf_seed_path(filename)*".mmn", "w") do data_Mmnkb
        @printf(data_Mmnkb, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
        @printf(data_Mmnkb, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

        for rank = 0:nprocs-1
            MPI_Nkpt = length(MPI_krange[rank+1])
            knum = MPkpts[rank+1]
            jldopen(_cwf_work_file(filename, "Mmnkb$rank"), "r") do file
                Mmnkb = file["Mmnkb"]
                _Write_Mmnkb(data_Mmnkb, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[1])
            end
        end
    end
end


function Generate_Mmnkb_Spindeg2(filename::AbstractString, BANDNUM, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    kmesh = mlwf_kpoints.kmesh
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    tot_bvector = mlwf_kpoints.tot_bvector
    frac_bv_int = mlwf_kpoints.frac_bv_int
    seed = _cwf_seed_path(filename)
    open(seed*"_1.mmn", "w") do data_Mmnkb1
        open(seed*"_2.mmn", "w") do data_Mmnkb2
            @printf(data_Mmnkb1, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
            @printf(data_Mmnkb1, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)
            @printf(data_Mmnkb2, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
            @printf(data_Mmnkb2, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

            for rank = 0:nprocs-1
                MPI_Nkpt = length(MPI_krange[rank+1])
                knum = MPkpts[rank+1]
                jldopen(_cwf_work_file(filename, "Mmnkb$rank"), "r") do file
                    Mmnkb = file["Mmnkb"]
                    _Write_Mmnkb(data_Mmnkb1, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[1])
                    _Write_Mmnkb(data_Mmnkb2, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb[2])
                end
            end
        end
    end
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


function Generate_Mmnkb_Col!(filename::AbstractString, MinN, MaxN, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

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
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    atv_ijk = material.atv_ijk
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
    work = zeros(ComplexF64, maximum(Total_NumOrbs), BANDNUM)
    Cnk1 = zeros(ComplexF64, fsize, fsize)
    Cnk2 = zeros(ComplexF64, fsize, fsize)

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        ik2 = ik + knum
        Read_Cnk_work!(filename, spin, ik2, Cnk1)
        @printf("\tmyrank = %3d %d %d/%d\n", myrank, spin, ik2, AllNkpt)
        for ib = 1:tot_bvector
            kk = kplusb[ik2][ib]+1
            @. k1 = MPI_kpts[ik]
            @. k2 = MPI_kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
            
            Read_Cnk_work!(filename, spin, kk, Cnk2)
            Calc_Mmnkb_Col!(Mmnkb[spin][ik][ib], MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], work, Cnk1, Cnk2)
        end
    end
    MPI.Barrier(comm)
end


function Generate_Mmnkb_NonCol!(filename::AbstractString, MinN, MaxN, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

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
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    atv_ijk = material.atv_ijk
    Nfsize = 2*fsize
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
    work = zeros(ComplexF64, maximum(Total_NumOrbs), BANDNUM)
    Cnk1 = zeros(ComplexF64, Nfsize, Nfsize)
    Cnk2 = zeros(ComplexF64, Nfsize, Nfsize)

    for ik = 1:MPI_Nkpt
        ik2 = ik + knum
        Read_Cnk_work!(filename, 1, ik2, Cnk1)
        @printf("\tmyrank = %3d %d %d/%d\n", myrank, 1, ik2, AllNkpt)
        for ib = 1:tot_bvector
            kk = kplusb[ik2][ib]+1
            @. k1 = MPI_kpts[ik]
            @. k2 = MPI_kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
            
            Read_Cnk_work!(filename, 1, kk, Cnk2)
            Calc_Mmnkb_NonCol!(Mmnkb[1][ik][ib], MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, fsize, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], work, Cnk1, Cnk2)
        end
    end
    MPI.Barrier(comm)
end


@timeit timer "Calc_Mmnkb" function Calc_Mmnkb_Col!(
    Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, work, Cnk1, Cnk2)

    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]
    offset = MinN - 1
    bands = offset+1:offset+BANDNUM

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn] + 1
        l1, l2, l3 = atv_ijk[cell]
        Anum = MP[atom]
        Bnum = MP[jatom]
        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        ex1 = cis(2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn)
        ex2 = cis(-2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn)
        block = OLPexp[atom][Rn]
        NO0, NO1 = size(block)
        arange = Anum+1:Anum+NO0
        brange = Bnum+1:Bnum+NO1

        C1a = view(Cnk1, arange, bands)
        C1b = view(Cnk1, brange, bands)
        C2a = view(Cnk2, arange, bands)
        C2b = view(Cnk2, brange, bands)

        work_a = view(work, 1:NO0, :)
        mul!(work_a, block, C2b)
        mul!(Mmnkb, adjoint(C1a), work_a, 0.5*ex1, 1.0)

        work_b = view(work, 1:NO1, :)
        mul!(work_b, transpose(block), C2a)
        mul!(Mmnkb, adjoint(C1b), work_b, 0.5*ex2, 1.0)
    end
end


@timeit timer "Calc_Mmnkb" function Calc_Mmnkb_NonCol!(
    Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, fsize, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, work, Cnk1, Cnk2)

    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]
    offset = MinN - 1
    bands = offset+1:offset+BANDNUM

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn] + 1
        l1, l2, l3 = atv_ijk[cell]
        Anum = MP[atom]
        Bnum = MP[jatom]
        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        ex1 = cis(2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn)
        ex2 = cis(-2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn)
        block = OLPexp[atom][Rn]
        NO0, NO1 = size(block)
        a1range = Anum+1:Anum+NO0
        a2range = Anum+fsize+1:Anum+fsize+NO0
        b1range = Bnum+1:Bnum+NO1
        b2range = Bnum+fsize+1:Bnum+fsize+NO1

        C1a1 = view(Cnk1, a1range, bands)
        C1a2 = view(Cnk1, a2range, bands)
        C1b1 = view(Cnk1, b1range, bands)
        C1b2 = view(Cnk1, b2range, bands)
        C2a1 = view(Cnk2, a1range, bands)
        C2a2 = view(Cnk2, a2range, bands)
        C2b1 = view(Cnk2, b1range, bands)
        C2b2 = view(Cnk2, b2range, bands)

        work_a = view(work, 1:NO0, :)
        mul!(work_a, block, C2b1)
        mul!(Mmnkb, adjoint(C1a1), work_a, 0.5*ex1, 1.0)
        mul!(work_a, block, C2b2)
        mul!(Mmnkb, adjoint(C1a2), work_a, 0.5*ex1, 1.0)

        work_b = view(work, 1:NO1, :)
        mul!(work_b, transpose(block), C2a1)
        mul!(Mmnkb, adjoint(C1b1), work_b, 0.5*ex2, 1.0)
        mul!(work_b, transpose(block), C2a2)
        mul!(Mmnkb, adjoint(C1b2), work_b, 0.5*ex2, 1.0)
    end
end
