include("../src/utils/Associated_Legendre.jl")
using Test

@testset "Associated_Legendre" begin
    l = 1
    m = 1
    x = 0.3
    result = Associated_Legendre(l, m, x)
    @test result ≈ -0.9539392014169457
end


@testset "Associated_Legendre2" begin
    l = 1
    m = 1
    x = 0.3
    result1, result2 = Associated_Legendre2(l, m, x)
    @test result1 ≈ -0.9539392014169457
    @test result2 ≈ -0.3
end
