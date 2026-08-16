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


function Set_MPI_iDM2iDM(MPI_iDM, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    Total_NumOrbs = system_grid.Total_NumOrbs
    

    iDM = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 2)
	for spin = 1:2
		iDM[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			iDM[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				iDM[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					iDM[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    MPI_iDM1 = MPI_iDM[1]
    MPI_iDM2 = MPI_iDM[2]
    hst = 0
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        iDM1 = iDM[1][atom][Rn]
        iDM2 = iDM[2][atom][Rn]
        for ist = 1:NO0, jst = 1:NO1
            hst += 1
            iDM1[ist][jst] = MPI_iDM1[hst]
            iDM2[ist][jst] = MPI_iDM2[hst]
        end
    end

    iDM1 = iDM[1]
    iDM2 = iDM[2]
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(iDM1[atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iDM2[atom][Rn][ist], MPI.SUM, comm)
    end
end


function Set_MPI_iHks2iHks(MPI_iHks, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    Total_NumOrbs = system_grid.Total_NumOrbs
    
    iHks = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
	for spin = 1:3
		iHks[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			iHks[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				iHks[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					iHks[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    MPI_iHks1 = MPI_iHks[1]
    MPI_iHks2 = MPI_iHks[2]
    MPI_iHks3 = MPI_iHks[3]
    hst = 0
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        iHks1 = iHks[1][atom][Rn]
        iHks2 = iHks[2][atom][Rn]
        iHks3 = iHks[3][atom][Rn]
        for ist = 1:NO0, jst = 1:NO1
            hst += 1
            iHks1[ist][jst] = MPI_iHks1[hst]
            iHks2[ist][jst] = MPI_iHks2[hst]
            iHks3[ist][jst] = MPI_iHks3[hst]
        end
    end

    iHks1 = iHks[1]
    iHks2 = iHks[2]
    iHks3 = iHks[3]
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(iHks1[atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iHks2[atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iHks3[atom][Rn][ist], MPI.SUM, comm)
    end

    return iHks
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