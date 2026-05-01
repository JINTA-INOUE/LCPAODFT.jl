include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Ngrid Si" begin
    
    Ecut = 150.0
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (24,24,24)

    Ecut = 200.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (27,27,27)

    Ecut = 250.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (30,30,30)

    Ecut = 300.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (32,32,32)

    Ecut = 350.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (35,35,35)

    Ecut = 400.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (36,36,36)

    Ecut = 450.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (40,40,40)

    Ecut = 500.0
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    @test Ngrid == (42,42,42)
end



@testset "Ngrid Graphene" begin
    
    Ecut = 150.0
    Lattice = [ 0.0000000000000 2.489015870 0.0; 
                2.1555509738426 1.244507935 0.0; 
                0.0000000000000 0.000000000 10.0]*"Ang"
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (16,16,75)

    Ecut = 200.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (18,18,84)

    Ecut = 250.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (21,21,96)

    Ecut = 300.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (24,24,105)

    Ecut = 350.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (24,24,112)

    Ecut = 400.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (25,25,120)

    Ecut = 450.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (28,28,128)

    Ecut = 500.0
    Ngrid = Calc_Ngrid(Ecut, Lattice.Latvecs)
    @test Ngrid == (30,30,135)
end