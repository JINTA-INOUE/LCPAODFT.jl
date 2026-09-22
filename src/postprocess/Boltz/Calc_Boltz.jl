function Calc_Boltz(boltz_setup::Boltz_Setup)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    kpoints = _boltz_kpoints(boltz_setup.kmesh)
    myrank == 0 && println(
        "<Calc_TDF_streaming>  Generate TDF using bounded k-point blocks")
    TDF_Energy, TDF, TDF_decomp =
        Calc_TDF_streaming(boltz_setup, kpoints)

    if myrank == 0
        println("<Calc_Sigma>")
        Calc_Sigma(boltz_setup, TDF_Energy, TDF)

        println("<Calc_SigmaS>")
        Calc_SigmaS(boltz_setup, TDF_Energy, TDF)

        println("<Calc_Seebeck>")
        Calc_Seebeck(boltz_setup, TDF_Energy, TDF)
    end

    if boltz_setup.decomp && myrank == 0
        println("<Calc_Sigma_decomp>")
        Calc_Sigma_decomp(boltz_setup, TDF_Energy, TDF_decomp)

        println("<Calc_Seebeck_decomp>")
        Calc_Seebeck_decomp(boltz_setup, TDF_Energy, TDF_decomp)
    end

    # Spectra are no longer needed after all output routines have returned.
    TDF_decomp = nothing
    TDF = nothing
    GC.gc()

    MPI.Barrier(comm)
    if myrank == 0
        println()
        println("Boltz completed with $nprocs MPI process(es) × " *
                "$(Threads.nthreads()) Julia thread(s)")
        @show LCPAODFT.timer
    end
end
