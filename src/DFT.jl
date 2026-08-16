function DFT(dft_setup::DFT_Setup)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    LCPAODFT.reset_timer!(LCPAODFT.timer)

    
    Latvecs = dft_setup.Latvecs
    Nspin = dft_setup.Nspin
    Natom = dft_setup.Natom
    Nspecies = dft_setup.Nspecies
    atom2spe = dft_setup.atom2spe
    Gxyz = dft_setup.Gxyz
    Atoms_pao = dft_setup.Atoms_pao
    Atoms_symbol = dft_setup.Atoms_symbol
    system = dft_setup.system
    xc_type = dft_setup.xc_type
    SpinPol = dft_setup.SpinPol
    SO_switch = dft_setup.SO_switch
    xc_type = dft_setup.xc_type
    Ngrid = dft_setup.Ngrid
    Grid_Origin = dft_setup.Grid_Origin
    kmesh = dft_setup.kmesh
    E_Temp = dft_setup.E_Temp
    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    Init_Mixing_weight = dft_setup.Init_Mixing_weight
    Min_Mixing_weight = dft_setup.Min_Mixing_weight
    Max_Mixing_weight = dft_setup.Max_Mixing_weight
    Num_Mixing_Pulay = dft_setup.Num_Mixing_Pulay
    Start_Pulay_SCF = dft_setup.Start_Pulay_SCF
    time_rev = dft_setup.time_rev
    cal_mode = dft_setup.cal_mode
    verbosity = dft_setup.verbosity
    fileout = dft_setup.fileout



    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)

    pao = Vector{PAO}(undef, Nspecies)
    pspot = Vector{Pspot}(undef, Nspecies)
    for spe = 1:Nspecies
        pspot[spe] = Read_VPS(Spe_symbol[spe], Spe_extra[spe], xc_type, SO_switch; verbosity)
        pao[spe] = Read_PAO(pspot[spe].Spe_Core_Charge, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe]; verbosity)
    end

    

    Atoms_Cut1 = zeros(Float64, Natom)
    Atoms_Core_Charge = zeros(Float64, Natom)
    Total_NumOrbs = zeros(Int32, Natom)
    fsize = 0
	for atom = 1:Natom
		spe = atom2spe[atom]
		Atoms_Cut1[atom] = pao[spe].Spe_Atom_Cut1
        Atoms_Core_Charge[atom] = pspot[spe].Spe_Core_Charge
        Total_NumOrbs[atom] = pao[spe].Spe_Total_NumOrbs
        fsize += Total_NumOrbs[atom]
	end


    Shift_K_Point = 1.0e-12
    kpoints = KPoints(kmesh, time_rev, Shift_K_Point)

    dft_options = DFT_Options(Mixing_method, SCF_criterion, SCF_max, 
                              Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
                              Num_Mixing_Pulay, -1, Start_Pulay_SCF, 3, Gxyz, false, time_rev)
    
    energy = Init_Energy()
    force = Init_Force(Natom, Gxyz, Atoms_symbol)

    ucell = UCell(Nspin, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
    electron = Electron(kpoints, cal_mode, SpinPol, E_Temp, Atoms_Core_Charge, fsize, system)

    
    KSsolve_SCF!(1, dft_setup, pao, pspot, ucell, electron, kpoints, dft_options, energy, force, false, fileout)



    myrank == 0 && println("")
    myrank == 0 && println("The calculation was normally finished.")


    MPI.Finalized()
end


function DFT(geoopt_setup::GeoOpt_Setup)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.get_num_threads()

    LCPAODFT.reset_timer!(LCPAODFT.timer)

    
    dft_setup = geoopt_setup.dft_setup
    Latvecs = dft_setup.Latvecs
    Recvecs = dft_setup.Recvecs
    Nspin = dft_setup.Nspin
    Natom = dft_setup.Natom
    Nspecies = dft_setup.Nspecies
    atom2spe = dft_setup.atom2spe
    Gxyz = dft_setup.Gxyz
    Atoms_pao = dft_setup.Atoms_pao
    Atoms_symbol = dft_setup.Atoms_symbol
    system = dft_setup.system
    xc_type = dft_setup.xc_type
    SpinPol = dft_setup.SpinPol
    SO_switch = dft_setup.SO_switch
    xc_type = dft_setup.xc_type
    Ngrid = dft_setup.Ngrid
    Grid_Origin = dft_setup.Grid_Origin
    kmesh = dft_setup.kmesh
    E_Temp = dft_setup.E_Temp
    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    Init_Mixing_weight = dft_setup.Init_Mixing_weight
    Min_Mixing_weight = dft_setup.Min_Mixing_weight
    Max_Mixing_weight = dft_setup.Max_Mixing_weight
    Num_Mixing_Pulay = dft_setup.Num_Mixing_Pulay
    Start_Pulay_SCF = dft_setup.Start_Pulay_SCF
    Extra_CHistory = geoopt_setup.Extra_CHistory
    time_rev = dft_setup.time_rev
    filename = dft_setup.filename
    verbosity = dft_setup.verbosity
    fileout = dft_setup.fileout
    send_email = dft_setup.send_email



    geo_optim = Init_Geo_Optim(Atoms_symbol, geoopt_setup)
    Geo_Opt_Max = geo_optim.Geo_Opt_Max
    


    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)

    pao = Vector{PAO}(undef, Nspecies)
    pspot = Vector{Pspot}(undef, Nspecies)
    for spe = 1:Nspecies
        pspot[spe] = Read_VPS(Spe_symbol[spe], Spe_extra[spe], xc_type, SO_switch; verbosity)
        pao[spe] = Read_PAO(pspot[spe].Spe_Core_Charge, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe]; verbosity)
    end

    

    Atoms_Cut1 = zeros(Float64, Natom)
    Atoms_Core_Charge = zeros(Float64, Natom)
    Total_NumOrbs = zeros(Int32, Natom)
    fsize = 0
	for atom = 1:Natom
		spe = atom2spe[atom]
		Atoms_Cut1[atom] = pao[spe].Spe_Atom_Cut1
        Atoms_Core_Charge[atom] = pspot[spe].Spe_Core_Charge
        Total_NumOrbs[atom] = pao[spe].Spe_Total_NumOrbs
        fsize += Total_NumOrbs[atom]
	end


    Gxyz_Optim = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Gxyz_Optim[atom] = zeros(Float64, 3)
        for xyz = 1:3
            Gxyz_Optim[atom][xyz] = Gxyz[atom][xyz]
        end
    end

    His_Gxyz = Vector{Vector{Float64}}(undef, Extra_CHistory)
    for i = 1:Extra_CHistory
        His_Gxyz[i] = zeros(Float64, 3*Natom)
    end
    k = 0
    for atom = 1:Natom
        His_Gxyz[1][k+1] = Gxyz[atom][1]
        His_Gxyz[1][k+2] = Gxyz[atom][2]
        His_Gxyz[1][k+3] = Gxyz[atom][3]
        k += 3
    end




    Shift_K_Point = 1.0e-12
    kpoints = KPoints(kmesh, time_rev, Shift_K_Point)

    dft_options = DFT_Options(Mixing_method, SCF_criterion, SCF_max, 
                              Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
                              Num_Mixing_Pulay, -1, Start_Pulay_SCF, Extra_CHistory, His_Gxyz, false, time_rev)
    
    electron = Electron(kpoints, SpinPol, E_Temp, Atoms_Core_Charge, fsize, system)

    energy = Init_Energy()
    force = Init_Force(Natom, Gxyz_Optim, Atoms_symbol)



    GeoOpt_po = false



    for GeoOpt_iter = 1:Geo_Opt_Max

        if myrank == 0
            println("\n\n")
            println("*******************************************************")
            println("\tGeoOpt_iter = $GeoOpt_iter")
            println("*******************************************************")
        end


        ucell = UCell(Nspin, Latvecs, Natom, atom2spe, Gxyz_Optim, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
        Set_Geo_Optim!(geo_optim, ucell.system_grid)

        # Be sure to save each step of the structural optimization process.
        KSsolve_SCF!(GeoOpt_iter, dft_setup, pao, pspot, ucell, electron, kpoints, dft_options, energy, force, true, false)



        Geo_Optim!(GeoOpt_iter, geo_optim, energy, force)

        Geo_Opt_convergence = geo_optim.Geo_Opt_convergence
        if Geo_Opt_convergence
            GeoOpt_po = true
            geo_optim.Geo_Opt_Conv_iter = GeoOpt_iter
            break
        else
            @. Gxyz_Optim = deepcopy(geo_optim.Gxyz)
            myrank == 0 && print("\n\n")
        end
    end


    


    if myrank == 0
        println("")
        if GeoOpt_po
            println("Geometric Optimizations convergence")
        else
            println("Geometric Optimizations not convergence")
        end

        ForceAll = force.ForceAll
        println("Final Atomic Forces (a.u.)")
        for atom = 1:Natom
            @printf("  atom = %d  %4s  %15.12f  %15.12f  %15.12f\n", atom, Atoms_symbol[atom], ForceAll[atom,1], ForceAll[atom,2], ForceAll[atom,3])
        end

        Gxyz_frac_Optim = Calc_Gxyz_frac(Natom, Gxyz_Optim, Recvecs)
        println("Final Atomic positions (Ang.)")
        for atom = 1:Natom
            @printf("  atom = %d  %4s  %15.12f  %15.12f  %15.12f\n",  atom, Atoms_symbol[atom], Gxyz_Optim[atom][1]/Ang_to_bohr, Gxyz_Optim[atom][2]/Ang_to_bohr, Gxyz_Optim[atom][3]/Ang_to_bohr)
        end

        println("Final Atomic positions (Frac)")
        for atom = 1:Natom
            @printf("  atom = %d  %4s  %15.12f  %15.12f  %15.12f\n",  atom, Atoms_symbol[atom], Gxyz_frac_Optim[atom][1], Gxyz_frac_Optim[atom][2], Gxyz_frac_Optim[atom][3])
        end
    end
    MPI.Barrier(comm)



    if myrank == 0
        Write_xyzfile(filename, Natom, Gxyz_Optim, force.ForceAll, Atoms_symbol)
    end
    MPI.Barrier(comm)


    if verbosity>=1
        myrank == 0 && println("")
        Print_TimerOutput(LCPAODFT.timer, comm)
    end
    MPI.Barrier(comm)


    myrank == 0 && println("")
    myrank == 0 && println("The calculation was normally finished.")


    if send_email && myrank == 0
        sendmail(filename, GeoOpt_po, geo_optim)
    end


    MPI.Finalized()
end
