include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "Spherical_Besselj" begin
    xmin = 1e-2
    xmax = 60.0
    N = 900
    dx = (xmax-xmin)/N
    x = zeros(Float64, N)
    for i = 1:N
        x[i] = xmin + (i-1)*dx
    end

    SphB0 = zeros(Float64, N)
    SphB1 = zeros(Float64, N)
    SphB2 = zeros(Float64, N)
    SphB3 = zeros(Float64, N)
    SphB4 = zeros(Float64, N)

    @. SphB0 = SphericalBesselj(0, x)
    @. SphB1 = SphericalBesselj(1, x)
    @. SphB2 = SphericalBesselj(2, x)
    @. SphB3 = SphericalBesselj(3, x)
    @. SphB4 = SphericalBesselj(4, x)

    @test SphB0 ≈ @. sin(x)/x
    @test SphB1 ≈ @. sin(x)/x^2 - cos(x)/x
    @test SphB2 ≈ @. sin(x)*(3/x^3 - 1/x) - 3*cos(x)/x^2
    @test SphB3 ≈ @. sin(x)*(15/x^4 - 6/x^2) - cos(x)*(15/x^3 - 1/x)
    @test SphB4 ≈ @. sin(x)*(105/x^5 - 45/x^3 + 1/x) - cos(x)*(105/x^4 - 10/x^2)
end