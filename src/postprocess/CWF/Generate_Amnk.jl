@timeit timer "Generate_Amnk" function Generate_Amnk(filename::AbstractString, SpinPol::String, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

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


function Generate_Amnk_Spindeg1(filename::AbstractString, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    open(_cwf_seed_path(filename)*".amn", "w") do data_Amnk
        @printf(data_Amnk, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
        @printf(data_Amnk, "%d %d %d\n", BANDNUM, Nkpt, Nwann)

        for rank = 0:nprocs-1
            MPI_Nkpt = length(MPI_krange[rank+1])
            knum = MPkpts[rank+1]
            jldopen(_cwf_work_file(filename, "Amnk$rank"), "r") do file
                Amnk = file["Amnk"]
                _Write_Amnk(data_Amnk, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[1])
            end
        end
    end
end


function Generate_Amnk_Spindeg2(filename::AbstractString, BANDNUM, Nwann, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = mlwf_kpoints.AllNkpt
    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    seed = _cwf_seed_path(filename)
    open(seed*"_1.amn", "w") do data_Amnk1
        open(seed*"_2.amn", "w") do data_Amnk2
            @printf(data_Amnk1, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
            @printf(data_Amnk1, "%d %d %d\n", BANDNUM, Nkpt, Nwann)
            @printf(data_Amnk2, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
            @printf(data_Amnk2, "%d %d %d\n", BANDNUM, Nkpt, Nwann)

            for rank = 0:nprocs-1
                MPI_Nkpt = length(MPI_krange[rank+1])
                knum = MPkpts[rank+1]
                jldopen(_cwf_work_file(filename, "Amnk$rank"), "r") do file
                    Amnk = file["Amnk"]
                    _Write_Amnk(data_Amnk1, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[1])
                    _Write_Amnk(data_Amnk2, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk[2])
                end
            end
        end
    end
end


function _Write_Amnk(data, MPI_Nkpt, knum, Nwann, BANDNUM, Amnk)
    @inbounds for ik = 1:MPI_Nkpt, m = 1:Nwann, μ = 1:BANDNUM
        @printf(data, "%5d %5d %5d %22.14f %22.14f\n", μ, m, ik+knum, real(Amnk[ik][μ][m]), imag(Amnk[ik][μ][m]))
    end
end
