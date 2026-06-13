function KSsolve_SCF!(
    GeoOpt_iter,
    dft_setup::DFT_Setup,
    pao::Vector{PAO}, pspot::Vector{Pspot}, 
    ucell::UCell, electron::CrystalBloch, 
    kpoints::KPoints,
    dft_options::DFT_Options,
    energy::Energy, force::Force, 
    work_fileout::Bool, fileout::Bool)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nspin = dft_setup.Nspin
    Natom = dft_setup.Natom
    gLatvecs = dft_setup.gLatvecs
    Init_Atoms_Nspin = dft_setup.Init_Atoms_Nspin
    Init_Atoms_Angle = dft_setup.Init_Atoms_Angle
    system = dft_setup.system
    xc_type = dft_setup.xc_type
    SpinPol = dft_setup.SpinPol
    Ngrid = dft_setup.Ngrid
    NN = prod(Ngrid)
    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    filename = dft_setup.filename
    restart = dft_setup.restart
    filepath = dft_setup.filepath
    verbosity = dft_setup.verbosity
    send_email = dft_setup.send_email
    system_grid = ucell.system_grid
    Total_Hsize = system_grid.Total_Hsize
    Atoms_Core_Charge = electron.Core_Charge
    Extra_CHistory = dft_options.Extra_CHistory
    His_Gxyz = dft_options.His_Gxyz


    work_dirname = pwd()*"/"*filename*"_work"
    filename_work_Hks = work_dirname*"/"*filename*"_work_Hks.jld2"

    filename_work_rho = Vector{String}(undef, Extra_CHistory)
    for i = 1:Extra_CHistory
        file = work_dirname*"/"*filename*"_work_rho$i.jld2"
        filename_work_rho[i] = file
    end




    
    Ham = Hamiltonian(SpinPol, pao, pspot, system_grid)
    
    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)
    ADensity_Grid, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, pao, pspot, ucell)



    dVHart_Grid = zeros(Float64, NN)
	Vpot_Grid = Vector{Vector{Float64}}(undef, Nspin)
	for spin = 1:Nspin
		Vpot_Grid[spin] = zeros(Float64, NN)
	end



    Hks = Vector{Vector{Float64}}(undef, Nspin)
    DM = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        Hks[spin] = zeros(Float64, Total_Hsize)
        DM[spin] = zeros(Float64, Total_Hsize)
    end

    if SpinPol == "nc"
        iDM = Vector{Vector{Float64}}(undef, 2)
        for spin = 1:2
            iDM[spin] = zeros(Float64, Total_Hsize)
        end
    else
        iDM = nothing
    end


    # Read restart file for Geometric Optimization
    SucceedReadingHksfile = 0
    SucceedReadingrhofile = 0
    if GeoOpt_iter > 1
        SucceedReadingHksfile = Read_restartFile_Hks!(filename_work_Hks, system_grid, Hks)
        extpln_coes = Extp_Charge(GeoOpt_iter, Extra_CHistory, Natom, His_Gxyz, system_grid.Gxyz)
        SucceedReadingrhofile = Read_restartFile_rho!(filename_work_rho, SpinPol, extpln_coes, Extra_CHistory, Nspin, Ngrid, ADensity_Grid, Density_Grid)
        if SucceedReadingHksfile == 1 || SucceedReadingrhofile == 1
            myrank == 0 && println("<Restart>  Found restart files")
        else
            myrank == 0 && println("<Restart>  Could not find restart files ($SucceedReadingHksfile, $SucceedReadingrhofile)")
        end
    end


    # Read restart file for SCF calc
    if restart
        Read_restartFile_DM!(filepath, system_grid, DM)

        myrank == 0 && println("<Set_Density_Grid>  Calculation Electron Density")
        
        if SpinPol == "off"
            Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
        elseif SpinPol == "on"
            Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
        elseif SpinPol == "nc"
            Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
            diagonalize_nc_density!(Density_Grid)
        end

        SucceedReadingrhofile = 1
    end



    xc_func = XC_Func(xc_type, SpinPol, Nspin, Ngrid, gLatvecs)
    dft_mixing = DFT_Mixing(Nspin, dft_options, system_grid)

    TotalZ = electron.TotalZ
    mulliken_charge = Mulliken_Charge(SpinPol, system_grid, TotalZ)



    scf_po = false
    dEele = 1.0
    pEele = 100.0



    myrank == 0 && println("\n")
    myrank == 0 && println("Start Self-Consistent Loop\n")
    
    
    for SCF_iter = 1:SCF_max

        myrank == 0 && @show SCF_iter

        if SCF_iter ≠ 1 || SucceedReadingrhofile == 1
            myrank == 0 && println("<Poisson>  Poisson's equation using FFT")
            Solve_Poisson!(SpinPol, dft_mixing, Density_Grid, ADensity_Grid, dVHart_Grid)
        end


        Set_XC_Grid!(xc_func, PCCDensity_Grid, Density_Grid)
        Set_Vpot_Grid!(xc_func, dVHart_Grid, Vpot_Grid)

        if SucceedReadingHksfile == 0 || SucceedReadingrhofile == 1
            myrank == 0 && println("<Set_Hamiltonian>  Hamiltonian matrix for dVH+Vxc ...")
            Set_Hamiltonian!(Ham, ucell, Orbs_Grid, Vpot_Grid, Hks)
        end
        SucceedReadingHksfile = 0
        SucceedReadingrhofile = 1
        

        if Mixing_method == "RMM-DIISH"
            dft_mixing.ChemP = electron.ChemP
            Mixing_H!(SCF_iter, Hks, dft_options, dft_mixing)
        end



        # solve Hc = ϵSc
        if system == "Cluster"
            error("not support")
        elseif system == "Crystal"
            myrank == 0 && println("<Crystal_DFT>  Solving the eigenvalue problem ...")
            Crystal_DFT!(Ham, system_grid, electron, kpoints, Hks, DM)
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

            myrank == 0 && println("<Set_Density_Grid>  Calculation Electron Density")

            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "nc"
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
            end

            # Solve_Poisson!(SpinPol, dft_mixing, Density_Grid, ADensity_Grid, dVHart_Grid)
            # Set_XC_Grid!(xc_func, PCCDensity_Grid, Density_Grid)
            # Set_Vpot_Grid!(xc_func, dVHart_Grid, Vpot_Grid)

            break
        end

        


        # mixing
        if Mixing_method ∈ ("Simple", "Kerker", "RMM-DIISK")
            error("not support $Mixing_method.")        
        elseif Mixing_method == "RMM-DIISH"

            update_Eele!(Eele, dft_mixing)

            myrank == 0 && println("<Set_Density_Grid>  Calculation Electron Density")

            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
            elseif SpinPol == "nc"
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
                diagonalize_nc_density!(Density_Grid)
            end
        else
            error("Please check Mixing_method")
        end



        Mulliken_Charge!(mulliken_charge, DM, Ham.OLP)
        

        if myrank == 0

            Total_Mul_up = sum(mulliken_charge.InitN_USpin)
            Total_Mul_dn = sum(mulliken_charge.InitN_DSpin)
            Total_SpinS = mulliken_charge.Total_SpinS

            @printf(" Sum of MulP: up    = %5.12f down           = %5.12f\n", Total_Mul_up, Total_Mul_dn)
            @printf("              total = %5.12f ideal(neutral) = %5.12f\n", Total_Mul_up+Total_Mul_dn, TotalZ)
            if SpinPol == "nc"
                Total_SpinAngle0 = mulliken_charge.Total_SpinAngle0
                Total_SpinAngle1 = mulliken_charge.Total_SpinAngle1

                @printf("       Total Spin Moment    (muB) %5.12f   Angles %5.12f %5.12f\n", 2*Total_SpinS, Total_SpinAngle0/pi*180, Total_SpinAngle1/pi*180)
            else
                @printf("       Total Spin Moment (muB) = %5.12f\n", 2*Total_SpinS)
            end

            @printf("Mixing_weight = %5.12f\n", dft_options.Mixing_weight)
            @printf("Eele = %5.12f  dEele = %5.12f\n", electron.Eele, dEele)
            @printf("ChemP = %5.12f\n", electron.ChemP)
            @printf("NormRD = %5.12f  Criterion = %5.12f\n", sqrt(abs(dft_mixing.NormRD[1])), SCF_criterion)
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

            @printf(" Sum of MulP: up    = %5.12f down           = %5.12f\n", Total_Mul_up, Total_Mul_dn)
            @printf("              total = %5.12f ideal(neutral) = %5.12f\n", Total_Mul_up+Total_Mul_dn, TotalZ)
            if SpinPol == "nc"
                Total_SpinAngle0 = mulliken_charge.Total_SpinAngle0
                Total_SpinAngle1 = mulliken_charge.Total_SpinAngle1

                @printf("       Total Spin Moment    (muB) %5.12f   Angles %5.12f %5.12f\n", 2*Total_SpinS, Total_SpinAngle0/pi*180, Total_SpinAngle1/pi*180)
            else
                @printf("       Total Spin Moment (muB) = %5.12f\n", 2*Total_SpinS)
            end
            @printf("Eele = %5.12f  dEele = %5.12f\n", electron.Eele, dEele)
            @printf("NormRD = %5.12f  Criterion = %5.12f\n", sqrt(abs(dft_mixing.NormRD[1])), SCF_criterion)
        else
            @printf("Eele = %5.12f  dEele = %5.12f\n", electron.Eele, dEele)
            @printf("NormRD = %5.12f  Criterion = %5.12f\n", sqrt(abs(dft_mixing.NormRD[1])), SCF_criterion)
            println("SCF not convergence")
        end


        system_charge = 0
        dipole_moment = Calc_dipole_moment(SpinPol, system_grid, Atoms_Core_Charge, system_charge, Density_Grid)
        if verbosity >= 1
            println("\n")
            println("Dipole moment")
            @printf("Total       %5.12f  %5.12f  %5.12f\n", dipole_moment[1,1], dipole_moment[1,2], dipole_moment[1,3])
            @printf("Core        %5.12f  %5.12f  %5.12f\n", dipole_moment[2,1], dipole_moment[2,2], dipole_moment[2,3])
            @printf("Electron    %5.12f  %5.12f  %5.12f\n", dipole_moment[3,1], dipole_moment[3,2], dipole_moment[3,3])
            @printf("Back ground %5.12f  %5.12f  %5.12f\n", dipole_moment[4,1], dipole_moment[4,2], dipole_moment[4,3])
        end
    end
    MPI.Barrier(comm)





    myrank == 0 && println("\n")
    myrank == 0 && println("<Energy> Energy calculation ...")
    
    Calc_iDM_Crystal_NonCollinear!(electron, kpoints, system_grid, iDM)

    Total_Energy!(energy, force, DM, iDM, 
                  ADensity_Grid, PCCDensity_Grid, Density_Grid, 
                  dVHart_Grid, Ham, system_grid, pao, pspot)
    energy.Eele = electron.Eele
    energy.ChemP = electron.ChemP
    if myrank == 0
        println("<Energy>")
        Print_Energy(energy)
    end
    MPI.Barrier(comm)




    myrank == 0 && println("\n")
    myrank == 0 && println("<Force> Force calculation ...")

    DM_Vec = Set_DM2DM_Vec(DM, system_grid)
    iDM_Vec = Set_DM2DM_Vec(iDM, system_grid)

    Force!(force, electron, kpoints,
           DM_Vec, iDM_Vec, Orbs_Grid,
           ADensity_Grid, PCCDensity_Grid, 
           dVHart_Grid, xc_func.Vxc_Grid, Vpot_Grid,
           Ham, ucell, pao, pspot)
    myrank == 0 && println("\n")
    


    # write $filename_work/$filename_work.jld2 file
    if work_fileout && myrank == 0
        WriteFile(filename_work_Hks, system_grid, energy, force, DM, Hks, Ham.iHNL)
        File_Shift(filename_work_rho, Extra_CHistory)
        WriteFile!(filename_work_rho[1], SpinPol, Ngrid, ADensity_Grid, Density_Grid)
    end
    MPI.Barrier(comm)


    # write .jld2 file
    if fileout && myrank == 0
        OLP_Vec = Set_HVNA2HVNA_Vec(Ham.OLP, system_grid)
        Hks_Vec = Set_DM2DM_Vec(Hks, system_grid)
        iHks_Vec = Set_DM2DM_Vec(Ham.iHNL, system_grid)
        WriteFile(Atoms_Core_Charge, mulliken_charge, dft_setup, system_grid, energy, force, DM_Vec, iDM_Vec, OLP_Vec, Hks_Vec, iHks_Vec)
    end
    MPI.Barrier(comm)



    if verbosity>=1 && myrank==0
        println("")
        println("")
        @show LCPAODFT.timer
    end



    if send_email && iszero((GeoOpt_iter-1)%5) && myrank == 0
        sendmail(filename, scf_po, dft_mixing)
    end



    MPI.Barrier(comm)
end

