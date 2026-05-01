function KSsolve_SCF(dft_setup::DFT_Setup)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.get_num_threads()

    LCPAODFT.reset_timer!(LCPAODFT.timer)


    Natom = dft_setup.Natom
    Nspecies = dft_setup.Nspecies
    Nspin = dft_setup.Nspin
    atom2spe = dft_setup.atom2spe
    Gxyz = dft_setup.Gxyz
    Latvecs = dft_setup.Latvecs
    gLatvecs = dft_setup.gLatvecs
    Recvecs = dft_setup.Recvecs
    Atoms_pao = dft_setup.Atoms_pao
    Atoms_Symbol = dft_setup.Atoms_symbol
    Init_Atoms_Nspin = dft_setup.Init_Atoms_Nspin
    Init_Atoms_Angle = dft_setup.Init_Atoms_Angle
    system = dft_setup.system
    xc_type = dft_setup.xc_type
    SpinPol = dft_setup.SpinPol
    SO_switch = dft_setup.SO_switch
    xc_type = dft_setup.xc_type
    Ngrid = dft_setup.Ngrid
    kmesh = dft_setup.kmesh
    E_Temp = dft_setup.E_Temp
    pao_file = dft_setup.pao_file
    pspot_file = dft_setup.pspot_file
    Grid_Origin = dft_setup.Grid_Origin
    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    Init_Mixing_weight = dft_setup.Init_Mixing_weight
    Min_Mixing_weight = dft_setup.Min_Mixing_weight
    Max_Mixing_weight = dft_setup.Max_Mixing_weight
    Num_Mixing_Pulay = dft_setup.Num_Mixing_Pulay
    Start_Pulay_SCF = dft_setup.Start_Pulay_SCF
    time_rev = dft_setup.time_rev
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
	for atom = 1:Natom
		spe = atom2spe[atom]
		Atoms_Cut1[atom] = pao[spe].Spe_Atom_Cut1
        Atoms_Core_Charge[atom] = pspot[spe].Spe_Core_Charge
        Total_NumOrbs[atom] = pao[spe].Spe_Total_NumOrbs
	end


    
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
    system_grid = ucell.system_grid
    Total_Hsize = system_grid.Total_Hsize
    
    Ham = Hamiltonian(SpinPol, pao, pspot, system_grid)
    
    xc_func = XC_Func(xc_type, SpinPol, Nspin, Ngrid, gLatvecs)

    

    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)
    ADensity_Grid, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, pao, pspot, ucell)



    dVHart_Grid = zeros(Float64, prod(Ngrid))
	Vpot_Grid = Vector{Vector{Float64}}(undef, Nspin)
	for spin = 1:Nspin
		Vpot_Grid[spin] = zeros(Float64, prod(Ngrid))
	end



    Hks = Vector{Vector{Float64}}(undef, Nspin)
    DM = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        Hks[spin] = zeros(Float64, Total_Hsize)
        DM[spin] = zeros(Float64, Total_Hsize)
    end



    kpoints = KPoints(system, SpinPol, kmesh, time_rev)

    dft_options = DFT_Options(Mixing_method, SCF_criterion, SCF_max, 
                              Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
                              Num_Mixing_Pulay, -1, Start_Pulay_SCF, false, time_rev)
    

    dft_mixing = DFT_Mixing(Nspin, dft_options, system_grid)
    electron = Electron(kpoints, SpinPol, E_Temp, Atoms_Core_Charge, sum(Total_NumOrbs), system)


    TotalZ = electron.TotalZ
    mulliken_charge = Mulliken_Charge(SpinPol, system_grid, TotalZ)



    scf_po = false
    dEele = 1.0
    pEele = 100.0

    


    myrank == 0 && println("\n")
    myrank == 0 && println("Start Self-Consistent Loop\n")
    
    
    for SCF_iter = 1:SCF_max

        myrank == 0 && @show SCF_iter

        if SCF_iter ≠ 1
            if Mixing_method == "RMM-DIISH"
                if SpinPol == "off"
                    Solve_Poisson!(dft_mixing, 2*Density_Grid[1]-2*ADensity_Grid, dVHart_Grid)
                elseif SpinPol ∈ ("on", "nc")
                    Solve_Poisson!(dft_mixing, Density_Grid[1]+Density_Grid[2]-2*ADensity_Grid, dVHart_Grid)
                end
            elseif Mixing_method ∈ ("Simple", "Kerker", "RMM-DIISK")
                error("not support")
            end
        end

        
        Set_XC_Grid!(xc_func, PCCDensity_Grid, Density_Grid)
        Set_Vpot_Grid!(xc_func, dVHart_Grid, Vpot_Grid)
        Set_Hamiltonian!(Ham, ucell, Orbs_Grid, Vpot_Grid, Hks)

        if Mixing_method == "RMM-DIISH"
            dft_mixing.ChemP = electron.ChemP
            Mixing_H!(SCF_iter, Hks, dft_options, dft_mixing)
        end



        # solve Hc = ϵSc
        if system == "Cluster"
            error("not support")
        elseif system == "Crystal"
            Crystal_DFT!(Ham, system_grid, electron, Hks, DM)
        end



        Eele = electron.Eele
        if SCF_iter ≠ 1
            dEele = abs(Eele - pEele)
        end
        pEele = Eele



        # check SCF-convergence
        if dEele < SCF_criterion && dft_mixing.NormRD[1] < 10*SCF_criterion
            scf_po = true
            dft_mixing.is_convergence = true
            dft_mixing.Conv_iter = SCF_iter

            # calculate Mulliken_Charge and Density_Grid using convergence eigen vector
            Mulliken_Charge!(mulliken_charge, DM, Ham.OLP)

            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "nc"
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
            end

            break
        end



        # mixing
        if Mixing_method ∈ ("Simple", "Kerker", "RMM-DIISK")
            error("not support $Mixing_method.")        
        elseif Mixing_method == "RMM-DIISH"

            update_Eele!(Eele, dft_mixing)

            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "nc"
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
                diagonalize_nc_density!( Density_Grid )
            end
        else
            error("Please check Mixing_method")
        end



        Mulliken_Charge!(mulliken_charge, DM, Ham.OLP)
        

        if myrank == 0

            Total_Mul_up = sum(mulliken_charge.InitN_USpin)
            Total_Mul_dn = sum(mulliken_charge.InitN_DSpin)
            Total_SpinS = mulliken_charge.Total_SpinS

            println(" Sum of MulP: up    = $Total_Mul_up down           = $Total_Mul_dn")
            println("              total = $(Total_Mul_up+Total_Mul_dn) ideal(neutral) = $TotalZ")
            if SpinPol == "nc"
                Total_SpinAngle0 = mulliken_charge.Total_SpinAngle0
                Total_SpinAngle1 = mulliken_charge.Total_SpinAngle1

                println("       Total Spin Moment    (muB) $(2*Total_SpinS)   Angles $(Total_SpinAngle0/pi*180) $(Total_SpinAngle1/pi*180)")
            else
                println("       Total Spin Moment (muB) = $(2*Total_SpinS)")
            end

            println("Mixing_weight = $(dft_options.Mixing_weight)")
            println("Eele = $Eele  dEele = $dEele")
            println("ChemP = $(electron.ChemP)")
            println("NormRD = $(sqrt(abs(dft_mixing.NormRD[1])))  Criterion = $SCF_criterion")
            println("")
        end
        MPI.Barrier(comm)
    end




    if myrank == 0
        if scf_po
            println("SCF convergence")

            Total_Mul_up = sum(mulliken_charge.InitN_USpin)
            Total_Mul_dn = sum(mulliken_charge.InitN_DSpin)
            Total_SpinS = mulliken_charge.Total_SpinS

            println(" Sum of MulP: up    = $Total_Mul_up down           = $Total_Mul_dn")
            println("              total = $(Total_Mul_up+Total_Mul_dn) ideal(neutral) = $TotalZ")
            if SpinPol == "nc"
                Total_SpinAngle0 = mulliken_charge.Total_SpinAngle0
                Total_SpinAngle1 = mulliken_charge.Total_SpinAngle1

                println("       Total Spin Moment    (muB) $(2*Total_SpinS)   Angles $(Total_SpinAngle0/pi*180) $(Total_SpinAngle1/pi*180)")
            else
                println("       Total Spin Moment (muB) = $(2*Total_SpinS)")
            end
            println("Eele = $(electron.Eele)  dEele = $dEele")
            println("NormRD = $(sqrt(abs(dft_mixing.NormRD[1])))  Criterion = $SCF_criterion")
        else
            println("SCF not convergence")
            println("Eele = $(electron.Eele)  dEele = $dEele")
            println("NormRD = $(sqrt(abs(dft_mixing.NormRD[1])))  Criterion = $SCF_criterion")
        end


        system_charge = 0
        dipole_moment = Calc_dipole_moment(SpinPol, system_grid, Atoms_Core_Charge, system_charge, Density_Grid)
        if verbosity >= 1
            println("\n")
            println("Dipole moment")
            println("Total       $(dipole_moment[1,1])  $(dipole_moment[1,2])  $(dipole_moment[1,3])")
            println("Core        $(dipole_moment[2,1])  $(dipole_moment[2,2])  $(dipole_moment[2,3])")
            println("Electron    $(dipole_moment[3,1])  $(dipole_moment[3,2])  $(dipole_moment[3,3])")
            println("Back ground $(dipole_moment[4,1])  $(dipole_moment[4,2])  $(dipole_moment[4,3])")
        end
    end
    MPI.Barrier(comm)



    
    myrank == 0 && println("\n")
    myrank == 0 && println("<Energy> Energy calculation ...")
    if SpinPol == "nc"
        iDM = Vector{Vector{Float64}}(undef, 2)
        for spin = 1:2
            iDM[spin] = zeros(Float64, Total_Hsize)
        end
        Calc_iDM_Crystal_NonCollinear!(electron, system_grid, iDM)
    else
        iDM = nothing
    end
    energy = Init_Energy()
    force = Init_Force(Natom, Gxyz, Atoms_Symbol)
    Total_Energy!(energy, force, DM, iDM, 
                  ADensity_Grid, PCCDensity_Grid, Density_Grid, 
                  dVHart_Grid, Ham, system_grid, pao, pspot)
    energy.Eele = electron.Eele
    energy.ChemP = electron.ChemP
    if myrank == 0
        println("<Energy>")
        Print_Energy(energy)
    end

    
    myrank == 0 && println("\n")
    myrank == 0 && println("<Force> Force calculation ...")



    DM_Vec = Set_DM2DM_Vec(DM, system_grid)
    if !isnothing(iDM)
        iDM_Vec = Set_DM2DM_Vec(iDM, system_grid)
    else
        iDM_Vec = nothing
    end

    OLP_Vec = Set_HVNA2HVNA_Vec(Ham.OLP, system_grid)
    Hks_Vec = Set_DM2DM_Vec(Hks, system_grid)
    if !isnothing(Ham.iHNL)
        iHks_Vec = Set_DM2DM_Vec(Ham.iHNL, system_grid)
    else
        iHks_Vec = nothing
    end


    Force!(force, electron,
           DM_Vec, iDM_Vec, Orbs_Grid,
           ADensity_Grid, PCCDensity_Grid, 
           dVHart_Grid, xc_func.Vxc_Grid, Vpot_Grid,
           Ham, ucell, pao, pspot)


    # write .jld2 file
    if fileout && myrank == 0
        WriteFile(Atoms_Core_Charge, mulliken_charge, dft_setup, system_grid, energy, force, DM_Vec, iDM_Vec, OLP_Vec, Hks_Vec, iHks_Vec)
    end


    if verbosity>=1 && myrank==0
        println("")
        @show LCPAODFT.timer
    end
    
    #=
    CpyCell = system_grid.CpyCell
    TCpyCell = (2*CpyCell + 1)^3 - 1
    material = LCPAO_model(
        Natom, Nspecies, Nspin, atom2spe,
        Atoms_Symbol, Atoms_Cut1, Atoms_pao, Atoms_Core_Charge,
        Init_Atoms_Nspin, Init_Atoms_Angle, mulliken_charge.Angle_Spin,
        Latvecs, Recvecs, Gxyz, TCpyCell, system_grid.atv, system_grid.atv_ijk, 
        system_grid.FNAN, system_grid.natn, system_grid.ncn, Total_NumOrbs, system_grid.MP, Grid_Origin,
        Ngrid, SO_switch, SpinPol, xc_type, time_rev, 
        E_Temp, kmesh, SCF_criterion, pao_file, pspot_file,
        OLP_Vec, Hks_Vec, iHks_Vec, DM_Vec, iDM_Vec, energy.ChemP, energy.Eele, energy.Etot, force.ForceAll
    )=#


    myrank == 0 && println("")
    myrank == 0 && println("The calculation was normally finished.")

    MPI.Finalized()
end
