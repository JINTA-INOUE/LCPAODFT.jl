include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_PCC_Force_code_warntype()
    verbosity = 0
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "off"
    Init_Atoms_Nspin = [[2.0,2.0],[2.0,2.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; verbosity)
    ADensity_Grid, PCCDensity_Grid, _ = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, [pao], [pspot], ucell)
    Vxc_Grid = [zeros(Float64, length(ADensity_Grid))]
    dVHart_Grid = zeros(Float64, length(ADensity_Grid))

    @code_warntype PCC_Force(SpinPol, ADensity_Grid, PCCDensity_Grid, Vxc_Grid, dVHart_Grid, [pao], [pspot], ucell)
end

check_PCC_Force_code_warntype()