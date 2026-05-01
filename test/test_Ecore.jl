include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Ecore Si" begin
    
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    Atom_Core_Charge = [4.0, 4.0]
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    system_grid = ucell.system_grid
    
    Ecore = Calc_Ecore(system_grid, Atom_Core_Charge)
    @test Ecore == 140.5574753035424
end