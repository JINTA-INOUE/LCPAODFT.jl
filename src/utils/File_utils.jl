# shift the index of stored data
function File_Shift(filename::Vector{String}, Extra_CHistory)
    if length(filename) ≠ Extra_CHistory
        error("please check filename, Extra_CHistory")
    end

    for i = Extra_CHistory-1:-1:1
        file1 = filename[i]
        file2 = filename[i+1]
        mv(file1, file2, force=true)
    end
end


function Generate_rhoFile(filename::String, Ngrid, Density_Grid)
    jldopen("$filename", "w") do file
        file["Ngrid"] = Ngrid
        file["Density_Grid"] = Density_Grid
        file["Dates"] = now()
    end
end