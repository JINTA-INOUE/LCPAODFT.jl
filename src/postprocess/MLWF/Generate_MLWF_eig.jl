@timeit timer "Generate_MLWF_eig" function Generate_MLWF_eig(filename::String, SpinPol::String, ChemP, BANDNUM, Nk, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        if SpinPol ∈ ("off", "nc")
            Generate_MLWF_Enk_Spindeg1(filename, BANDNUM, Nk, ChemP, mlwf_kpoints)
        elseif SpinPol == "on"
            Generate_MLWF_Enk_Spindeg2(filename, BANDNUM, Nk, ChemP, mlwf_kpoints)
        else
            error("please check SpinPol.")
        end
    end
    MPI.Barrier(comm)
end


function Generate_MLWF_Enk_Spindeg1(filename::String, BANDNUM, Nk, ChemP, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    work_dirname = pwd()*"/"*filename*"_work_cwf"

    Nk_spin = Nk[1,:,1]

    data_Enk = open(filename*".eig", "w")

    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Enk$rank.jld2"
        data = jldopen(work_file, "r")
        Enk = data["Enk"]
        _Write_MLWF_Enk(data_Enk, MPI_Nkpt, knum, BANDNUM, Nk_spin, ChemP, Enk[1])
        close(data)
    end
    close(data_Enk)
end


function Generate_MLWF_Enk_Spindeg2(filename::String, BANDNUM, Nk, ChemP, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    work_dirname = pwd()*"/"*filename*"_work_cwf"

    Nk_up = Nk[1,:,1]
    Nk_dn = Nk[2,:,1]

    data_Enk1 = open(filename*"_1.eig", "w")
    data_Enk2 = open(filename*"_2.eig", "w")
    for rank = 0:nprocs-1
        MPI_Nkpt = length(MPI_krange[rank+1])
        knum = MPkpts[rank+1]
        work_file = work_dirname*"/"*filename*"_Enk$rank.jld2"
        data = jldopen(work_file, "r")
        Enk = data["Enk"]
        _Write_MLWF_Enk(data_Enk1, MPI_Nkpt, knum, BANDNUM, Nk_up, ChemP, Enk[1])
        _Write_MLWF_Enk(data_Enk2, MPI_Nkpt, knum, BANDNUM, Nk_dn, ChemP, Enk[2])
        close(data)
    end
    close(data_Enk1)
    close(data_Enk2)
end


function _Write_MLWF_Enk(data, MPI_Nkpt, knum, BANDNUM, Nk, ChemP, Enk)
    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        @printf(data, "%5d %5d %22.14f\n", μ, ik+knum, (Enk[ik][μ+Nk[ik]]-ChemP)*27.2113845)
    end
end