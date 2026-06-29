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
        # Read CWF_model
        model = 2
    end

    return model
end