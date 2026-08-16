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

        myrank == 0 && println("<Set_CWF_ExpnCoef>")
        CWF_ExpnCoef = Set_CWF_ExpnCoef(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk)
            
        if CWF_Wannier
            myrank == 0 && println("<Set_CWF_Grid>")
            Set_CWF_Grid(cwf_setup, CWF_ExpnCoef)
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
    MPI.Finalized()
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
    

    if CWF_Wannier || write_coef

        myrank == 0 && println("<Set_CWF_ExpnCoef>")
        CWF_ExpnCoef = Set_CWF_ExpnCoef(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk)
            
        if CWF_Wannier
            myrank == 0 && println("<Set_CWF_Grid>")
            Set_CWF_Grid(cwf_setup, CWF_ExpnCoef)
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
    MPI.Finalized()
end
