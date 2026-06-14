function DosMain(filepath::String, kmesh, Erange::Vector{Float64}; mode::String="Dos")

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    if nprocs > 1
        error("please run serial.")
    end

    if length(kmesh) ≠ 3
        error("pleas check kmesh")
    end

    if length(Erange) ≠ 2
        error("pleas check Erange")
    end


    mode = lowercase(mode)
    if mode ∉ ("all", "dos", "pdos")
        error("please check mode.")
    end


    filename, _ = splitext(basename(filepath))
    println("filename = $filename")
    
    model = select_model(filepath)
    if model == 1
        material = Load_LCPAODFT_model(filepath)
        Print_LCPAO_model(filepath, material)
    elseif model == 2
        material = Load_CWF_model(filepath)
        Print_CWF_model(filepath, material)
    else
        error("please check filepath")
    end
    

    

    Dos_Erange = zeros(Float64, 2)
    @. Dos_Erange = Erange/eV2Hartree


    KP_flag = "Gcenter"
    kpoints = KPoints(kmesh, false, 0.0; KP_flag)


    iemin, iemax = Calc_Band_size(material, Dos_Erange)
    neg = iemax-iemin+1
    @show iemin, iemax, neg

    Enk, Cnk = Calc_Enk_Cnk_Dos(material, kpoints, iemin, iemax)


    if mode ∈ ("all", "pdos")
        EVec = Calc_EVec_Dos(material, kpoints, Cnk, iemin, iemax)
    end
    

    
    if mode == "all"
        Calc_DosMain(filename, material, Enk, neg, kmesh, Dos_Erange)
        Calc_PDosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
    elseif mode == "dos"
        Calc_DosMain(filename, material, Enk, neg, kmesh, Dos_Erange)
    elseif mode == "pdos"
        Calc_PDosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
    end
end
