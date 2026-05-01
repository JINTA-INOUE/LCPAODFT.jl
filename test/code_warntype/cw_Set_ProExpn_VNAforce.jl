include("../../src/LCPAODFT.jl")
using .LCPAODFT

# Si case

function check_Set_DS_VNAforce_code_warntype()
    verbosity = 0
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "off"
    Init_Atoms_Nspin = [[2.0,2.0],[2.0,2.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    Total_NumOrbs = [13,13]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
    
    system_grid = ucell.system_grid
    Total_Hsize = system_grid.Total_Hsize
    Natom = system_grid.Natom
    Nspecies = 1
    FNAN = system_grid.FNAN
    Total_NumOrbs = system_grid.Total_NumOrbs

    maxL = maximum([pao.Spe_MaxL_Basis for spe = 1:Nspecies]) + BufferL_ProVNA
    VNATotal_Num = (maxL+1)^2 * maxM

    maxFNAN = maximum(FNAN)
	maxTotal_NumOrbs = maximum(Total_NumOrbs)
	DS_VNAforce = Vector{Vector{Vector{Vector{Vector{Float32}}}}}(undef, 4)
	for xyz = 1:4
		DS_VNAforce[xyz] = Vector{Vector{Vector{Vector{Float32}}}}(undef, Natom+1)
		for atom = 1:Natom+1
			if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
			DS_VNAforce[xyz][atom] = Vector{Vector{Vector{Float32}}}(undef, fan)
			for Rn = 1:fan
				DS_VNAforce[xyz][atom][Rn] = Vector{Vector{Float32}}(undef, NO0)
				for ist = 1:NO0
					DS_VNAforce[xyz][atom][Rn][ist] = zeros(Float32, VNATotal_Num)
				end
			end
		end
	end

    @code_warntype Set_DS_VNAforce!(DS_VNAforce, [pao], [pspot], system_grid)
end



function check_Set_HVNA2_3_code_warntype()
    verbosity = 0
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "off"
    Init_Atoms_Nspin = [[2.0,2.0],[2.0,2.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    Total_NumOrbs = [13,13]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
    
    system_grid = ucell.system_grid
    Natom = system_grid.Natom
	FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HVNA2force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
    HVNA3force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
	for xyz = 1:3
		HVNA2force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		HVNA3force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			HVNA2force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			HVNA3force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				HVNA2force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					HVNA2force[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[atom])
				end

                HVNA3force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[natn[atom][Rn]])
                for ist = 1:Total_NumOrbs[natn[atom][Rn]]
					HVNA3force[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    @code_warntype Set_HVNA2_3force!(HVNA2force, HVNA3force, [pao], [pspot], system_grid)
end


# check_Set_DS_VNAforce_code_warntype()
check_Set_HVNA2_3_code_warntype()