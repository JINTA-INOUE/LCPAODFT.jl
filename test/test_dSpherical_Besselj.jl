include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "dSpherical_Besselj" begin
    
    xmin = 1e-2
    xmax = 60.0
    N = 900
    lmax = 3
    dx = (xmax-xmin)/N
    x = zeros(Float64, N)
    for i = 1:N
        x[i] = xmin + (i-1)*dx
    end

    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*lmax+3)
    dSphB_l = zeros(Float64, 2*lmax+3)


    dSphB = zeros(Float64, lmax+1, N)
    for i = 1:N
        Calc_SphericalBesselj!(lmax, x[i], tsb, SphB_l, dSphB_l)
        for l = 1:lmax+1
            dSphB[l,i] = dSphB_l[l]
        end
    end

    ref_SphB0 = zeros(Float64, N)
    ref_SphB1 = zeros(Float64, N)
    ref_SphB2 = zeros(Float64, N)
    ref_SphB3 = zeros(Float64, N)
    ref_SphB4 = zeros(Float64, N)
    @. ref_SphB0 = sin(x)/x
    @. ref_SphB1 = sin(x)/x^2 - cos(x)/x
    @. ref_SphB2 = sin(x)*(3/x^3 - 1/x) - 3*cos(x)/x^2
    @. ref_SphB3 = sin(x)*(15/x^4 - 6/x^2) - cos(x)*(15/x^3 - 1/x)
    @. ref_SphB4 = sin(x)*(105/x^5 - 45/x^3 + 1/x) - cos(x)*(105/x^4 - 10/x^2)

    ref_dSphB0 = zeros(Float64, N)
    ref_dSphB1 = zeros(Float64, N)
    ref_dSphB2 = zeros(Float64, N)
    ref_dSphB3 = zeros(Float64, N)

    @. ref_dSphB0 = -ref_SphB1
    @. ref_dSphB1 = -ref_SphB2 + ref_SphB1/x
    @. ref_dSphB2 = -ref_SphB3 + 2*ref_SphB2/x
    @. ref_dSphB3 = -ref_SphB4 + 3*ref_SphB3/x

    @test isapprox(dSphB[1,:], ref_dSphB0, atol=1e-10)
    @test isapprox(dSphB[2,:], ref_dSphB1, atol=1e-10)
    @test isapprox(dSphB[3,:], ref_dSphB2, atol=1e-6)
    @test isapprox(dSphB[4,:], ref_dSphB3, atol=1e-6)
end