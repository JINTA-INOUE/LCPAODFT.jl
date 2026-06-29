function Write_Variable_work_file(filename::String, myrank, variable, name::String)
    work_dirname = pwd()*"/"*filename*"_work_cwf"
    work_file = work_dirname*"/"*filename*"_$(name)$(myrank).jld2"
    jldopen(work_file, "w") do file
        file["$name"] = variable
    end
end