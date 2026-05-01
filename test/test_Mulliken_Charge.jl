include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "Initial Mulliken_Charge Si" begin

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
    TotalZ = 8.0
    Total_NumOrbs = [3,3]
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs )
    system_grid = ucell.system_grid

    mulliken_charge = Mulliken_Charge(SpinPol, system_grid, TotalZ)

    @test mulliken_charge.Natom == Natom
    @test mulliken_charge.Nspin == 1
    @test mulliken_charge.SpinPol == SpinPol
    @test mulliken_charge.FNAN == system_grid.FNAN
    @test mulliken_charge.natn == system_grid.natn
    @test mulliken_charge.Total_NumOrbs == Total_NumOrbs
    @test mulliken_charge.MulP == [0.0, 0.0, 0.0, 0.0]
    @test mulliken_charge.DecMulP == [[[0.0,0.0,0.0], [0.0,0.0,0.0]]]
    @test mulliken_charge.InitN_USpin == zeros(Float64, Natom)
    @test mulliken_charge.InitN_DSpin == zeros(Float64, Natom)
    @test mulliken_charge.Angle_Spin == [[0.0,0.0], [0.0,0.0]]
    @test mulliken_charge.Total_SpinS == 0.0
    @test mulliken_charge.Total_SpinAngle0 == 0.0
    @test mulliken_charge.Total_SpinAngle1 == 0.0
    @test mulliken_charge.TZ == 8.0
end