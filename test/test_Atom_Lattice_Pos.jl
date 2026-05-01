include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "Atompos" begin
    atompos = [[0.0,0.0,0.0]]*"Frac"

    @test atompos.Natom == 1
    @test atompos.Gxyz == [[0.0,0.0,0.0]]
    @test atompos.Gxyz_frac == [[0.0,0.0,0.0]]
    @test atompos.unit == "frac"

    atompos = [[0.0,0.0,0.0]]*"Ang"

    @test atompos.Natom == 1
    @test atompos.Gxyz == [[0.0,0.0,0.0]]
    @test atompos.Gxyz_frac == [[-99999.0,-99999.0,-99999.0]]
    @test atompos.unit == "ang"

    atompos = [[0.0,0.0,0.0]]*"AU"

    @test atompos.Natom == 1
    @test atompos.Gxyz == [[0.0,0.0,0.0]]
    @test atompos.Gxyz_frac == [[-99999.0,-99999.0,-99999.0]]
    @test atompos.unit == "au"


    atompos = [[0.0,0.0,0.0],[0.5,0.5,0.5]]*"Frac"

    @test atompos.Natom == 2
    @test atompos.Gxyz == [[0.0,0.0,0.0],[0.5,0.5,0.5]]
    @test atompos.Gxyz_frac == [[0.0,0.0,0.0],[0.5,0.5,0.5]]
    @test atompos.unit == "frac"

    atompos = [[0.0,0.0,0.0],[1.0,1.0,1.0]]*"Ang"

    @test atompos.Natom == 2
    @test atompos.Gxyz == [[0.0,0.0,0.0],[1.8897259886,1.8897259886,1.8897259886]]
    @test atompos.Gxyz_frac == [[-99999.0,-99999.0,-99999.0],[-99999.0,-99999.0,-99999.0]]
    @test atompos.unit == "ang"

    atompos = [[0.0,0.0,0.0],[1.0,1.0,1.0]]*"AU"

    @test atompos.Natom == 2
    @test atompos.Gxyz == [[0.0,0.0,0.0],[1.0,1.0,1.0]]
    @test atompos.Gxyz_frac == [[-99999.0,-99999.0,-99999.0],[-99999.0,-99999.0,-99999.0]]
    @test atompos.unit == "au"
end



@testset "Lattice" begin
    Lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]*"Ang"

    @test Lattice.Latvecs == [1.8897259886 0.0 0.0; 0.0 1.8897259886 0.0; 0.0 0.0 1.8897259886]
    @test Lattice.unit == "ang"
    @test Lattice.dim == 3

    Lattice = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]*"AU"

    @test Lattice.Latvecs == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    @test Lattice.unit == "au"
    @test Lattice.dim == 3
end