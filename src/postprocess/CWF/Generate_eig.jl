@timeit timer "Generate_eig" function Generate_eig(filename::AbstractString, SpinPol::String, ChemP, MinN, MaxN, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        if SpinPol ∈ ("off", "nc")
            Generate_Enk_Spindeg1(filename, MinN, MaxN, ChemP, mlwf_kpoints)
        elseif SpinPol == "on"
            Generate_Enk_Spindeg2(filename, MinN, MaxN, ChemP, mlwf_kpoints)
        else
            error("please check SpinPol.")
        end
    end
    MPI.Barrier(comm)
end


function Generate_Enk_Spindeg1(filename::AbstractString, MinN, MaxN, ChemP, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    open(_cwf_seed_path(filename)*".eig", "w") do data_Enk
        for rank = 0:nprocs-1
            MPI_Nkpt = length(MPI_krange[rank+1])
            knum = MPkpts[rank+1]
            jldopen(_cwf_work_file(filename, "Enk$rank"), "r") do file
                Enk = file["Enk"]
                _Write_Enk(data_Enk, MPI_Nkpt, knum, MinN, MaxN, ChemP, Enk[1])
            end
        end
    end
end


function Generate_Enk_Spindeg2(filename::AbstractString, MinN, MaxN, ChemP, mlwf_kpoints::MLWF_KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_krange = mlwf_kpoints.MPI_krange
    MPkpts = mlwf_kpoints.MPkpts
    seed = _cwf_seed_path(filename)
    open(seed*"_1.eig", "w") do data_Enk1
        open(seed*"_2.eig", "w") do data_Enk2
            for rank = 0:nprocs-1
                MPI_Nkpt = length(MPI_krange[rank+1])
                knum = MPkpts[rank+1]
                jldopen(_cwf_work_file(filename, "Enk$rank"), "r") do file
                    Enk = file["Enk"]
                    _Write_Enk(data_Enk1, MPI_Nkpt, knum, MinN, MaxN, ChemP, Enk[1])
                    _Write_Enk(data_Enk2, MPI_Nkpt, knum, MinN, MaxN, ChemP, Enk[2])
                end
            end
        end
    end
end


function _Write_Enk(data, MPI_Nkpt, knum, MinN, MaxN, ChemP, Enk)
    BANDNUM = MaxN - MinN + 1
    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        @printf(data, "%5d %5d %22.14f\n", μ, ik+knum, (Enk[ik][μ+MinN-1]-ChemP)*27.2113845)
    end
end
