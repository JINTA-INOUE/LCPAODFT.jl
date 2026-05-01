include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Calc_Atomic_Den2 Si" begin
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)

    Spe_Atomic_Den2 = Calc_Atomic_Den2(pao, pspot)
    
    @test length(Spe_Atomic_Den2) == pao.Spe_Num_Mesh_PAO+4
    @test Spe_Atomic_Den2[1:10] ≈ [0.05876999307780021, 0.05876999334333089, 0.05876999360886157, 0.058769993885316715, 0.0587699941731452, 0.05876999447281635, 0.05876999478481549, 0.058769995109651085, 0.058769995447851116, 0.05876999579996492]
end