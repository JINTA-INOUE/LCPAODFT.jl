function Set_DM_Vec2DM(DM, system_grid::System_Grid)

    Nspin = length(DM)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    Total_Hsize = system_grid.Total_Hsize

    DM_1D = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
		DM_1D[spin] = zeros(Float64, Total_Hsize)
	end
    
    for spin = 1:Nspin
        counts = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            counts += 1
            DM_1D[spin][counts] = DM[spin][atom][Rn][ist][jst]
        end
	end

    return DM_1D
end


function Set_DM2DM_Vec(DM, system_grid::System_Grid)

    if isnothing(DM)
        return nothing
    end

    Nspin = length(DM)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    DM_Vec = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspin)
	for spin = 1:Nspin
		DM_Vec[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			DM_Vec[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				DM_Vec[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					DM_Vec[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    for spin = 1:Nspin
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            hst += 1
            DM_Vec[spin][atom][Rn][ist][jst] = DM[spin][hst]
        end
	end

    return DM_Vec
end


function Set_HNL2HNL_Vec(HNL, system_grid::System_Grid)

    if isnothing(HNL)
        return nothing
    end
    
    Nspin = length(HNL)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HNL_Vec = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspin)
    for spin = 1:Nspin
        HNL_Vec[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
        for atom = 1:Natom
            HNL_Vec[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                HNL_Vec[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    HNL_Vec[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end

    
    for spin = 1:Nspin
        counts = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            counts += 1
            HNL_Vec[spin][atom][Rn][ist][jst] = HNL[spin][counts]
        end
    end


    return HNL_Vec
end


function Set_HVNA2HVNA_Vec(HVNA, system_grid::System_Grid)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HVNA_Vec = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
    for atom = 1:Natom
        HVNA_Vec[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            HVNA_Vec[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                HVNA_Vec[atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
            end
        end
    end

    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        counts += 1
        HVNA_Vec[atom][Rn][ist][jst] = HVNA[counts]
    end


    return HVNA_Vec
end