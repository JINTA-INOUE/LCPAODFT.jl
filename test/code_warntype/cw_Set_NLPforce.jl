include("../../src/LCPAODFT.jl")
using .LCPAODFT

# Si case

function check_Set_NLPforce_code_warntype()
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
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    VPS_j_Num = zeros(Int64, Natom)
    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot.Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end

        VPS_j_Num[atom] = pspot.VPS_j_dependency
        NLTotal_Num[atom] = tot
    end

    maxFNAN = maximum(FNAN)
	maxTotal_NumOrbs = maximum(Total_NumOrbs)
    maxVPS_j_Num = maximum(VPS_j_Num)
	maxNLTotal_Num = maximum(NLTotal_Num)


    NLPforce = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, 4)
    for xyz = 1:4
        NLPforce[xyz] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Natom+1)
        for atom = 1:Natom+1
            if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
            
            NLPforce[xyz][atom] = Vector{Vector{Vector{Vector{Float64}}}}(undef, fan)
            for Rn = 1:fan
                if atom == Natom+1
                    VPS_j_dependency = maxVPS_j_Num
                    NO1 = maxNLTotal_Num
                else
                    jatom = natn[atom][Rn]
                    VPS_j_dependency = VPS_j_Num[jatom]
                    NO1 = NLTotal_Num[jatom]
                end
                
                NLPforce[xyz][atom][Rn] = Vector{Vector{Vector{Float64}}}(undef, NO0)
                for ist = 1:NO0
                    NLPforce[xyz][atom][Rn][ist] = Vector{Vector{Float64}}(undef, VPS_j_dependency+1)
                    for so = 1:VPS_j_dependency+1
                        NLPforce[xyz][atom][Rn][ist][so] = zeros(Float64, NO1)
                    end
                end
            end
        end
    end

    @code_warntype Set_NLPforce!(NLPforce, [pao], [pspot], system_grid)
end

check_Set_NLPforce_code_warntype()