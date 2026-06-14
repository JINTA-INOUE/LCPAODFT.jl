function select_model(filepath::String)

    base_filepath = basename(filepath)
    file = split(base_filepath, ".")

    if file[end] ≠ "jld2"
        error("please check filepath.")
    end
    if length(file) == 2
        # Read LCPAO_model
        model = 1
    else
        # Read CWF_model
        model = 2
    end

    return model
end
