include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Calc_Spe_Vna Si" begin
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)

    Spe_Vna = Calc_Spe_Vna(pao, pspot)
    
    @test length(Spe_Vna) == pspot.Spe_Num_Mesh_VPS
    @test Spe_Vna[1:10] ≈ [-1.8456313944608693, -1.8456314220284078, -1.8456314071900546, -1.8456314120862185, -1.8456313954871812, -1.845631411940559, -1.8456314179089275, -1.845631384629665, -1.8456314220560053, -1.8456313981677068]
end