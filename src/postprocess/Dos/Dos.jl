function DosMain(filepath::String, kmesh, Erange::Vector{Float64};
                 mode::String="Dos", de_Dos=0.01)

    Threads.nthreads() == 1 || error(
        "MPI-flat DosMain requires exactly one Julia thread per MPI process; " *
        "start Julia with --threads=1",
    )
    provided_thread_level = MPI.Init(; threadlevel=:single)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    BLAS.set_num_threads(1)
    start_time = time()

    length(kmesh) == 3 || error("kmesh must contain three dimensions")
    length(Erange) == 2 || error("Erange must contain its lower and upper limits")
    de_Dos > 0 || error("de_Dos must be positive")

    mode = lowercase(mode)
    mode in ("all", "dos", "pdos") ||
        error("mode must be one of: all, dos, pdos")

    filename, _ = splitext(basename(filepath))
    if myrank == 0
        println("filename = $filename")
        println("<DOS MPI-flat configuration>")
        println("\t$nprocs MPI processes × 1 Julia thread")
        println("\t$(BLAS.get_num_threads()) BLAS thread per process")
        println("\tMPI thread level: $provided_thread_level")
    end

    # The model is read independently by each process. This avoids serializing
    # a large, deeply nested Julia object and keeps all later matrix data local.
    model = select_model(filepath)
    if model == 1
        material = Load_LCPAODFT_model(filepath)
        if myrank == 0
            Print_LCPAO_model(filepath, material)
        end
    elseif model == 2
        material = Load_CWF_model(filepath)
        if myrank == 0
            Print_CWF_model(filepath, material)
        end
    else
        error("unsupported model file: $filepath")
    end

    Dos_Erange = Erange ./ eV2Hartree
    if myrank == 0
        println("Erange = $Erange")
        println("de_Dos = $de_Dos")
        println("kmesh = $kmesh")
    end

    kpoints = _dos_kpoints(kmesh)
    halo_plan = _dos_halo_plan(kpoints, kmesh)

    iemin, iemax = Calc_Band_size(material, Dos_Erange)
    neg = iemax - iemin + 1

    Enk, Cnk = Calc_Enk_Cnk_Dos(material, kpoints, iemin, iemax)
    if mode in ("all", "pdos")
        EVec = Calc_EVec_Dos(material, kpoints, Cnk, iemin, iemax)
    end
    Cnk = nothing

    Enk = _dos_exchange_halo(Enk, halo_plan)
    if mode in ("all", "pdos")
        EVec = _dos_exchange_halo(EVec, halo_plan)
        # _dos_report_distributed_memory(halo_plan, kpoints.AllNkpt, Enk, EVec)
    else
        # _dos_report_distributed_memory(halo_plan, kpoints.AllNkpt, Enk)
    end

    if mode == "all"
        Calc_DosMain(filename, material, Enk, halo_plan, neg, kmesh,
                     Dos_Erange; de_Dos)
        Calc_PDosMain(filename, material, Enk, EVec, halo_plan, neg, kmesh,
                      Dos_Erange; de_Dos)
    elseif mode == "dos"
        Calc_DosMain(filename, material, Enk, halo_plan, neg, kmesh,
                     Dos_Erange; de_Dos)
    else
        Calc_PDosMain(filename, material, Enk, EVec, halo_plan, neg, kmesh,
                      Dos_Erange; de_Dos)
    end

    MPI.Barrier(comm)
    if myrank == 0
        # println("DOS calculation completed in $(round(time() - start_time; digits=3)) s")
    end
end
