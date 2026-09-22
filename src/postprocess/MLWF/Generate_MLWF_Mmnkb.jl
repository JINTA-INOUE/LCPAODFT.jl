@timeit timer "Generate_MLWF_Mmnkb" function Generate_MLWF_Mmnkb(mlwf_setup::MLWF_Setup, mlwf_kpoints::MLWF_KPoints, Nk, BANDNUM)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    filename = mlwf_setup.filename
    material = mlwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    MPI_Nkpt = mlwf_kpoints.MPI_Nkpt
    tot_bvector = mlwf_kpoints.tot_bvector
    bvector = mlwf_kpoints.bvector

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
        Generate_MLWF_Mmnkb_Col!(filename, BANDNUM, Nk, OLPexp, Mmnkb, material, mlwf_kpoints)
    elseif SpinPol == "nc"
        Generate_MLWF_Mmnkb_NonCol!(filename, BANDNUM, Nk, OLPexp, Mmnkb, material, mlwf_kpoints)
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


function Generate_MLWF_Mmnkb_Spindeg1(filename::String, BANDNUM, mlwf_kpoints::MLWF_KPoints)

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


function Generate_MLWF_Mmnkb_Spindeg2(filename::String, BANDNUM, mlwf_kpoints::MLWF_KPoints)

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


function _Write_MLWF_Mmnkb(data, BANDNUM, MPI_Nkpt, knum, kmesh, tot_bvector, frac_bv_int, Mmnkb)
    
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


function Generate_MLWF_Mmnkb_Col!(filename::String, BANDNUM, Nk, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

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
        MinN = Nk[spin,ik,1] + 1
        Read_Cnk_work!(filename, spin, ik2, Cnk1)
        @printf("\tmyrank = %3d %d %d/%d\n", myrank, spin, ik2, AllNkpt)
        for ib = 1:tot_bvector
            kk = kplusb[ik2][ib]+1
            @. k1 = MPI_kpts[ik]
            @. k2 = MPI_kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
        
            Read_Cnk_work!(filename, spin, kk, Cnk2)
            Calc_Mmnkb_Col!(Mmnkb[spin][ik][ib], MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, 
                            MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], work, Cnk1, Cnk2)
        end
    end
    MPI.Barrier(comm)
end


function Generate_MLWF_Mmnkb_NonCol!(filename::String, BANDNUM, Nk, OLPexp, Mmnkb, material::LCPAO_model, mlwf_kpoints::MLWF_KPoints)

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
        MinN = Nk[1,ik,1] + 1
        Read_Cnk_work!(filename, 1, ik2, Cnk1)
        @printf("\tmyrank = %3d %d/%d\n", myrank, ik2, AllNkpt)
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