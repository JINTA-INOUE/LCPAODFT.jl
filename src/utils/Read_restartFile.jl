@timeit timer "Read_restartFile" function Read_restartFile_DM!(filepath::String, system_grid::System_Grid, DM)
    data = jldopen(filepath, "r")
    DM_Vec = data["DM"]
    DM .= Set_DM_Vec2DM(DM_Vec, system_grid)
end


@timeit timer "Read_restartFile" function Read_restartFile_Hks!(filepath::String, system_grid::System_Grid, Hks)
    
    data = jldopen(filepath, "r")
    old_Natom = data["Natom"]
    old_TCpyCell = data["TCpyCell"]
    old_FNAN = data["FNAN"]
    old_natn = data["natn"]
    old_ncn = data["ncn"]
    old_Total_NumOrbs = data["Total_NumOrbs"]
    old_atv_ijk = data["atv_ijk"]
    old_Total_Hsize = data["Total_Hsize"]

    SucceedReadingHksfile = 0
    Natom = system_grid.Natom
    CpyCell = system_grid.CpyCell
    TCpyCell = (2*CpyCell + 1)^3 - 1
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    Total_NumOrbs = system_grid.Total_NumOrbs
    atv_ijk = system_grid.atv_ijk
    Total_Hsize = system_grid.Total_Hsize

    # check system_grid
    if Natom ≠ old_Natom || TCpyCell ≠ old_TCpyCell || Total_Hsize ≠ old_Total_Hsize
        close(data)
        return SucceedReadingHksfile
    end

    if FNAN ≠ old_FNAN || natn ≠ old_natn || ncn ≠ old_ncn
        close(data)
        return SucceedReadingHksfile
    end

    if Total_NumOrbs ≠ old_Total_NumOrbs || atv_ijk ≠ old_atv_ijk
        close(data)
        return SucceedReadingHksfile
    end


    Hks .= data["Hks"]
    close(data)
    SucceedReadingHksfile = 1
    
    
    return SucceedReadingHksfile
end


@timeit timer "Read_restartFile" function Read_restartFile_rho!(
    filepath::Vector{String}, 
    SpinPol::String, 
    extpln_coes, Extra_CHistory, 
    Nspin, Ngrid, 
    ADensity_Grid, Density_Grid)
    
    SucceedReadingrhofile = 0    
    @show extpln_coes

    data = jldopen(filepath[begin], "r")
    old_Ngrid = data["Ngrid"]
    if old_Ngrid ≠ Ngrid
        close(data)
        return SucceedReadingrhofile
    end
    SucceedReadingrhofile = 1


    NN = prod(Ngrid)
    Density_Grid_tmp = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        Density_Grid_tmp[spin] = zeros(Float64, NN)
    end


    # Read Density_Grid
    Density_Grid_tmp .= data["Density_Grid"]
    close(data)


    # extrapolate Density_Grid
    if SpinPol ∈ ("off", "on")
        for spin = 1:Nspin
            @. Density_Grid[spin] = ADensity_Grid + extpln_coes[1]*Density_Grid_tmp[spin]
        end
    elseif SpinPol == "nc"
        for spin = 1:2
            @. Density_Grid[spin] = ADensity_Grid + extpln_coes[1]*Density_Grid_tmp[spin]
        end
        @. Density_Grid[3] = extpln_coes[1]*Density_Grid_tmp[3]
        @. Density_Grid[4] = extpln_coes[1]*Density_Grid_tmp[4]
    else
        error("please check SpinPol")
    end
    

    for i = 2:Extra_CHistory
        data = jldopen(filepath[i], "r")
        Density_Grid_tmp .= data["Density_Grid"]
        for spin = 1:Nspin, j = 1:NN
            Density_Grid[spin][j] += extpln_coes[i]*Density_Grid_tmp[spin][j]
        end
        close(data)
    end

    
    return SucceedReadingrhofile
end