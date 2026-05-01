include("../../src/LCPAODFT.jl")
using .LCPAODFT

# Si case

function check_Set_dOrbitals_Grid_code_warntype()
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
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs, verbosity)
    
    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    GridN_Atom = ucell.GridN_Atom

    dOrbs_Grid = Vector{Vector{Vector{Vector{Float64}}}}(undef, 3)
    for xyz = 1:3
        dOrbs_Grid[xyz] = Vector{Vector{Vector{Float64}}}(undef, Natom)
        for atom = 1:Natom
            dOrbs_Grid[xyz][atom] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                dOrbs_Grid[xyz][atom][ist] = zeros(Float64, GridN_Atom[atom])
            end
        end
    end

    @code_warntype Set_dOrbitals_Grid!(dOrbs_Grid, [pao], ucell)
end

check_Set_dOrbitals_Grid_code_warntype()