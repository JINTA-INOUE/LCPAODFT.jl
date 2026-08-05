function Calc_Boltz(boltz_setup::Boltz_Setup)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    material = boltz_setup.material

    kpoints = _boltz_kpoints(boltz_setup.kmesh)
    halo_plan = _boltz_halo_plan(kpoints, boltz_setup.kmesh)

    myrank == 0 && println("<Calc_Enk_Cnk>")
    Enk, Cnk = Calc_Enk_Cnk_Boltz(material, kpoints)

    myrank == 0 && println("<Calc_Vnk!>")
    Vnk = Calc_Vnk!(material, kpoints, Enk, Cnk)

    if boltz_setup.decomp
        myrank == 0 && println("<Calc_EVec>")
        EVec = Calc_EVec(boltz_setup, kpoints, Cnk)
    end
    Cnk = nothing

    Enk_halo = _boltz_exchange_halo(Enk, halo_plan)
    Vnk_halo = _boltz_exchange_halo(Vnk, halo_plan)
    if boltz_setup.decomp
        EVec_halo = _boltz_exchange_halo(EVec, halo_plan)
        # _boltz_report_distributed_memory(halo_plan, kpoints.AllNkpt, Enk_halo, Vnk_halo, EVec_halo)
    else
        # _boltz_report_distributed_memory(halo_plan, kpoints.AllNkpt, Enk_halo, Vnk_halo)
    end

    myrank == 0 && println("<Calc_TDF>  Generate TDF using Tetrahedron method")
    TDF_Energy, TDF = Calc_TDF(boltz_setup, Enk_halo, Vnk_halo, halo_plan)

    if myrank == 0
        println("<Calc_Sigma>")
        Calc_Sigma(boltz_setup, TDF_Energy, TDF)

        println("<Calc_SigmaS>")
        Calc_SigmaS(boltz_setup, TDF_Energy, TDF)

        println("<Calc_Seebeck>")
        Calc_Seebeck(boltz_setup, TDF_Energy, TDF)
    end

    if boltz_setup.decomp
        myrank == 0 && println("<Calc_TDF_decomp>  Generate TDF using Tetrahedron method")
        TDF_Energy, TDF_decomp = Calc_TDF_decomp(boltz_setup, Enk_halo, EVec_halo, Vnk_halo, halo_plan)

        if myrank == 0
            println("<Calc_Sigma_decomp>")
            Calc_Sigma_decomp(boltz_setup, TDF_Energy, TDF_decomp)

            println("<Calc_Seebeck_decomp>")
            Calc_Seebeck_decomp(boltz_setup, TDF_Energy, TDF_decomp)
        end
    end

    MPI.Barrier(comm)
    if myrank == 0
        println()
        println("Boltz MPI-flat completed with $nprocs MPI processes")
        @show LCPAODFT.timer
    end
end
