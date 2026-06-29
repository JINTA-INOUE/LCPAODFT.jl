function Generate_CWF(cwf_setup::CWF_Setup)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    
    material = cwf_setup.material
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    
    spinsize = cwf_setup.spinsize
    kmesh = cwf_setup.kmesh
    Ngsize = cwf_setup.Ngsize
    weight_type = cwf_setup.weight_type
    MLWF_kpts = cwf_setup.MLWF_kpts
    Dis_Energy = cwf_setup.Dis_Energy
    CWF_HmnR = cwf_setup.CWF_HmnR
    CWF_Wannier = cwf_setup.CWF_Wannier
    CWF_SOC = cwf_setup.CWF_SOC
    CWF2MLWF = cwf_setup.CWF2MLWF
    write_coef = cwf_setup.write_coef
    filename = cwf_setup.filename
    verbose = cwf_setup.verbose
    
    
    Shift_K_Point = 0.0
    KP_flag = "Gcenter"
    kpoints = KPoints(kmesh, false, Shift_K_Point; KP_flag)
    

    myrank == 0 && println("<Calc_Enk_Cnk>")
    Enk, Cnk = Calc_Enk_Cnk(material, kpoints, 2)
    

    
    MinN, MaxN, BANDNUM = Calc_BANDNUM_KS_state(weight_type, MLWF_kpts, Dis_Energy, material)
    myrank == 0 && Check_CWF_Band(weight_type, MinN, MaxN, BANDNUM, Ngsize)
    MPI.Barrier(comm)


    
    myrank == 0 && println("<Calc_Amnk>")
    Amnk = Calc_Amnk(MinN, BANDNUM, Enk, Cnk, cwf_setup, kpoints)


    myrank == 0 && println("<Calc_Smk_Umnk>")
    Σmk, Umnk = Calc_Smk_Umnk(BANDNUM, Amnk, Ngsize, material, kpoints)


    myrank == 0 && println("<Calc_DMfunc>")
    DMfunc = Calc_DMfunc(spinsize, Ngsize, Σmk, kpoints)
    

    if CWF_HmnR

        myrank == 0 && println("<Calc_HmnR>")
        NCell, cell_list, cell_list_ijk = Get_cell_list(kmesh)
        HmnR = zeros(ComplexF64, Ngsize, Ngsize, NCell, spinsize)
        Calc_HmnR!(HmnR, spinsize, MinN, MaxN, NCell, Ngsize, cell_list_ijk, Umnk, Enk, kpoints)

        myrank == 0 && println("<Write_CWF_HmnR>")
        myrank == 0 && Write_CWF_HmnR(cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk, HmnR)
        MPI.Barrier(comm)
    end



    if CWF_Wannier || write_coef

        Latvecs = material.Latvecs
        Natom = material.Natom
        Nspecies = material.Nspecies
        atom2spe = material.atom2spe
        Gxyz = material.Gxyz
        Atoms_pao = material.Atoms_pao
        Total_NumOrbs = material.Total_NumOrbs
        Ecut = cwf_setup.Ecut
        Ngrid = Calc_Ngrid(Ecut, Latvecs)
        Atoms_Cut1 = material.Atoms_Cut1
        Grid_Origin = material.Grid_Origin

        Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)


        pao = Vector{PAO}(undef, Nspecies)
        for spe = 1:Nspecies
            pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
        end
        
        ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs)
        Orbs_Grid = Set_Orbitals_Grid(pao, ucell)

        CWF_ExpnCoef = Set_CWF_ExpnCoef(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk)
            
        if CWF_Wannier
            Set_CWF_Grid(CWF_ExpnCoef, Orbs_Grid, ucell, cwf_setup)
        end
    end


    
    if CWF2MLWF
        myrank == 0 && println("<CWF2Wannier90>")
        work_dirname = pwd()*"/"*filename*"_work_cwf"
        myrank == 0 && mkpath(work_dirname)
        MPI.Barrier(comm)

        Write_Cnk_work(filename, SpinPol, Cnk, kpoints)
        Write_Variable_work_file(filename, myrank, Enk, "Enk")
        Write_Variable_work_file(filename, myrank, Amnk, "Amnk")
        CWF2Wannier90(MinN, MaxN, cwf_setup)
        # rm(work_dirname, force=true)
    end


    if CWF_SOC
        myrank == 0 && println("<Calc_HmnR>")
        NCell, cell_list, cell_list_ijk = Get_cell_list(kmesh)
        HmnR = zeros(ComplexF64, Ngsize, Ngsize, NCell, spinsize)
        Calc_HmnR!(HmnR, spinsize, MinN, MaxN, NCell, Ngsize, cell_list_ijk, Umnk, Enk, kpoints)

        myrank == 0 && println("<Calc_CWF_SOC_Strength>")
        Wannier_SOC = Calc_CWF_SOC_Strength(cwf_setup, kpoints, MinN, MaxN, Umnk, Cnk)

        myrank == 0 && println("<Write_CWF_HmnR>")
        myrank == 0 && Write_CWF_HmnR(cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk, HmnR, Wannier_SOC)
        MPI.Barrier(comm)
    end
    



    if verbose>=1 && myrank==0
        println("")
        @show LCPAODFT.timer
    end
    MPI.Barrier(comm)
end


function Generate_CWF(cwf_setup::CWF_Setup_MO)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    
    material = cwf_setup.material
    SpinPol = material.SpinPol
    xc_type = material.xc_type
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    
    spinsize = cwf_setup.spinsize
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type
    MLWF_kpts = cwf_setup.MLWF_kpts
    kmesh = cwf_setup.kmesh
    Ngsize = cwf_setup.Ngsize
    CWF_HmnR = cwf_setup.CWF_HmnR
    CWF_Wannier = cwf_setup.CWF_Wannier
    CWF_SOC = cwf_setup.CWF_SOC
    CWF2MLWF = cwf_setup.CWF2MLWF
    write_coef = cwf_setup.write_coef
    guide_out = cwf_setup.guide_out
    filename = cwf_setup.filename
    verbose = cwf_setup.verbose
    
    

    Shift_K_Point = 1.0e-12
    KP_flag = "Gcenter"
    kpoints = KPoints(kmesh, false, Shift_K_Point; KP_flag)

    
    
    myrank == 0 && println("<Calc_Enk_Cnk> for <Set_CWF_Guiding_MOs>")
    Enk, Cnk = Calc_Enk_Cnk(material, kpoints, 2)
    

    MinN, MaxN, BANDNUM = Calc_BANDNUM_KS_state(weight_type, MLWF_kpts, Dis_Energy, material)
    myrank == 0 && Check_CWF_Band(weight_type, MinN, MaxN, BANDNUM, Ngsize)
    MPI.Barrier(comm)


    myrank == 0 && println("<Set_CWF_Guiding_MOs>")
    CWF_Guiding_MOs = Set_CWF_Guiding_MOs(cwf_setup, kpoints, Enk, Cnk)



    Shift_K_Point = 0.0
    KP_flag = "Gcenter"
    kpoints = KPoints(kmesh, false, Shift_K_Point; KP_flag)

    myrank == 0 && println("<Calc_Enk_Cnk!>")
    Calc_Enk_Cnk!(material, kpoints, Enk, Cnk)

    
    
    myrank == 0 && println("<Calc_Amnk>")
    Amnk = Calc_Amnk(MinN, BANDNUM, Enk, Cnk, cwf_setup, kpoints, CWF_Guiding_MOs)

    myrank == 0 && println("<Calc_Smk_Umnk>")
    Σmk, Umnk = Calc_Smk_Umnk(BANDNUM, Amnk, Ngsize, material, kpoints)

    myrank == 0 && println("<Calc_DMfunc>")
    DMfunc = Calc_DMfunc(spinsize, Ngsize, Σmk, kpoints)
    


    if CWF_HmnR

        myrank == 0 && println("<Calc_HmnR>")
        NCell, cell_list, cell_list_ijk = Get_cell_list(kmesh)
        HmnR = zeros(ComplexF64, Ngsize, Ngsize, NCell, spinsize)
        Calc_HmnR!(HmnR, spinsize, MinN, MaxN, NCell, Ngsize, cell_list_ijk, Umnk, Enk, kpoints)

        myrank == 0 && println("<Write_CWF_HmnR>")
        myrank == 0 && Write_CWF_HmnR(cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk, HmnR)
        MPI.Barrier(comm)
    end
    

    #=
    if guide_out
        Latvecs = material.Latvecs
        Natom = material.Natom
        Nspecies = material.Nspecies
        atom2spe = material.atom2spe
        Gxyz = material.Gxyz
        Atoms_pao = material.Atoms_pao
        Total_NumOrbs = material.Total_NumOrbs
        Ecut = cwf_setup.Ecut
        Ngrid = Calc_Ngrid(Ecut, Latvecs)
        Atoms_Cut1 = material.Atoms_Cut1
        Grid_Origin = material.Grid_Origin

        Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)


        pao = Vector{PAO}(undef, Nspecies)
        for spe = 1:Nspecies
            pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
        end
        
        ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs)
        Orbs_Grid = Set_Orbitals_Grid(pao, ucell)

        Set_CWF_GuideGrid(CWF_Guiding_MOs, Orbs_Grid, ucell, cwf_setup)        
    end
    =#


    if CWF_Wannier || write_coef

        Latvecs = material.Latvecs
        Natom = material.Natom
        Nspecies = material.Nspecies
        atom2spe = material.atom2spe
        Gxyz = material.Gxyz
        Atoms_pao = material.Atoms_pao
        Total_NumOrbs = material.Total_NumOrbs
        Ecut = cwf_setup.Ecut
        Ngrid = Calc_Ngrid(Ecut, Latvecs)
        Atoms_Cut1 = material.Atoms_Cut1
        Grid_Origin = material.Grid_Origin

        Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)


        pao = Vector{PAO}(undef, Nspecies)
        for spe = 1:Nspecies
            pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
        end
        
    
        ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs)
        Orbs_Grid = Set_Orbitals_Grid(pao, ucell)

        CWF_ExpnCoef = Set_CWF_ExpnCoef(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk)
         
        if CWF_Wannier
            Set_CWF_Grid(CWF_ExpnCoef, Orbs_Grid, ucell, cwf_setup)
        end
    end 



    if CWF2MLWF
        myrank == 0 && println("<CWF2Wannier90>")
        work_dirname = pwd()*"/"*filename*"_work_cwf"
        myrank == 0 && mkpath(work_dirname)
        MPI.Barrier(comm)

        Write_Cnk_work(filename, SpinPol, Cnk, kpoints)
        Write_Variable_work_file(filename, myrank, Enk, "Enk")
        Write_Variable_work_file(filename, myrank, Amnk, "Amnk")
        CWF2Wannier90(MinN, MaxN, cwf_setup)
        # rm(work_dirname, force=true)
    end



    if verbose>=1 && myrank==0
        println("")
        @show LCPAODFT.timer
    end
    MPI.Barrier(comm)
end
