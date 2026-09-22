function select_model(filepath::String)

    base_filepath = basename(filepath)
    file = split(base_filepath, ".")
    file_end = file[end]

    if file_end ∉ ("jld2", "scfout")
        error("please check filepath.")
    end

    if length(file) == 2 
        if file_end == "scfout"
            model = 0
        else file_end == "jld2"     # Read LCPAO_model
            model = 1
        end
    else
        if file[2] ∈ ("CWF", "CWF_SOC")
            model = 2   # Read CWF_model
        elseif file[2] ∈ ("MLWF", "MLWF_SOC")
            model = 3   # Read MLWF_model
        else
            error("please check filename")
        end
    end

    return model
end