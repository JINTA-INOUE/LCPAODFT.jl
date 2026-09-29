function w90_tools(mlwf_setup::MLWF_Setup, modes::Vector{String}, filepath::Vector{String})

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    if nprocs >= 2
        error("please run serial.")
    end

    @. modes = lowercase(modes)
    for mode in modes
        if mode ∉ ("wf", "soc")
            error("please check modes")
        end

        if mode == "wf"
            if isnothing(mlwf_setup.MLWF_Plot_Cube) || isnothing(mlwf_setup.MLWF_Plot_SuperCells)
                println("not define MLWF_Plot_Cube and MLWF_Plot_SuperCells")
                error("please check input")
            end
        end
    end


    
    for file in filepath
        _, ext = splitext(basename(file))
        if ext ≠ ".chk"
            error("please check filepath")
        end
    end

    Nfile = length(filepath)
    chk = Vector{w90_Chk}(undef, Nfile)
    for (spin,file) in enumerate(filepath)
        chk[spin] = Read_chk(file)
    end


    material = mlwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    MLWF_kmesh = mlwf_setup.MLWF_kmesh
    WANNUM = mlwf_setup.WANNUM
    filename = mlwf_setup.filename

    for spin = 1:Nfile
        if chk[spin].checkpoint ≠ "postwann"
            error("please check input file")
        end

        if chk[spin].Nkpt ≠ prod(MLWF_kmesh)
            error("please check input")
        end

        if chk[spin].WANNUM ≠ WANNUM
            error("please check input")
        end
    end

    Shift_K_Point = 0.0
    KP_flag = "Gcenter"
    kpoints = KPoints(MLWF_kmesh, false, Shift_K_Point; KP_flag)

    myrank == 0 && println("<Calc_Enk_Cnk>")
    Enk, Cnk = Calc_Enk_Cnk(material, kpoints, 2)

    Read_Cnk!(filename*"_AO_Poly", material, kpoints, Cnk)
    
    myrank == 0 && println("<Set_Outer_InnerNk>")
    BANDNUM, Nk = Set_Outer_InnerNk(Enk, kpoints, mlwf_setup)

    for spin = 1:Nfile
        if chk[spin].BANDNUM ≠ BANDNUM
            error("please check input")
        end
    end

    for mode in modes
        if mode == "wf"
            MLWF_ExpnCoef = Set_MLWF_ExpnCoef(mlwf_setup, kpoints, chk, BANDNUM, Nk, Cnk)
            Set_MLWF_Grid(mlwf_setup, MLWF_ExpnCoef)
        end
    end
    MPI.Barrier(comm)
end


function Read_Cnk!(filename, material, kpoints, Cnk)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nkpt = kpoints.Nkpt

    for spin = 1:spinsize, ik = 1:Nkpt
        file = joinpath(pwd(), filename*"_work_cwf", filename*"_Cnk$(spin)_$ik.jld2")
        data = jldopen(file, "r")
        Cnk[spin][ik] = data["Cnk"]
        close(data)
    end
end