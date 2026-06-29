@timeit timer "Generate_Amnk" function Generate_Amnk(filename::String, SpinPol::String, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        if SpinPol ∈ ("off", "nc")
            Generate_Amnk_Spindeg1(filename, BANDNUM, Nwann, mlwf_kpoints)
        elseif SpinPol == "on"
            Generate_Amnk_Spindeg2(filename, BANDNUM, Nwann, mlwf_kpoints)
        else
            error("please check SpinPol.")
        end
    end
    MPI.Barrier(comm)
end


function Generate_Amnk_Spindeg1(filename::String, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    work_dirname = pwd()*"/"*filename*"_work_cwf"

    data_Amnk = open(filename*".amn", "w")
    @printf(data_Amnk, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
    @printf(data_Amnk, "%d %d %d\n", BANDNUM, Nkpt, Nwann)

    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Amnk$rank.jld2"
        data = jldopen(work_file, "r")
        Amnk = data["Amnk"]
        _Write_Amnk(data_Amnk, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[1])
        close(data)
    end
    close(data_Amnk)
end


function Generate_Amnk_Spindeg2(filename::String, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    work_dirname = pwd()*"/"*filename*"_work_cwf"

    data_Amnk1 = open(filename*"_1.amn", "w")
    data_Amnk2 = open(filename*"_2.amn", "w")
    @printf(data_Amnk1, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
    @printf(data_Amnk1, "%d %d %d\n", BANDNUM, Nkpt, Nwann)
    @printf(data_Amnk2, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
    @printf(data_Amnk2, "%d %d %d\n", BANDNUM, Nkpt, Nwann)

    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Amnk$rank.jld2"
        data = jldopen(work_file, "r")
        Amnk = data["Amnk"]
        _Write_Amnk(data_Amnk1, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[1])
        _Write_Amnk(data_Amnk2, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[2])
        close(data)
    end
    close(data_Amnk1)
    close(data_Amnk2)
end


function _Write_Amnk(data, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk)
    @inbounds for ik = 1:MPI_Nkpt, m = 1:Nwann, μ = 1:BANDNUM
        @printf(data, "%5d %5d %5d %22.14f %22.14f\n", μ, m, ik+knum, real(Amnk[ik][μ][m]), imag(Amnk[ik][μ][m]))
    end
end