include("../../src/LCPAODFT.jl")
using .LCPAODFT

# Si case

function check_Set_Nonlocal_code_warntype()
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

    HNL = Vector{Vector{Float64}}(undef, 1)
	HNL[1] = zeros(Float64, Total_Hsize)
	iHNL = nothing
    @code_warntype Set_Nonlocal!("off", HNL, iHNL, [pao], [pspot], system_grid)
end

check_Set_Nonlocal_code_warntype()