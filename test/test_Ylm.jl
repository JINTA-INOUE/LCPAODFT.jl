include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Ylm_complex" begin

    atol = 1e-12
    theta = rand()
    phi = rand()
    Ylm = Ylm_table(0,0,theta,phi)
    
    Y00 = Ylm_table(0,0,theta,phi)
    @test isapprox(Ylm, 0.5*sqrt(1.0/pi) + im*0.0, atol=atol)

    # l = 1
    Y10 = Ylm_table(1, 0, theta, phi)
    @test isapprox(Y10, sqrt(3/(4*pi)) * cos(theta) + im*0.0, atol=atol)

    Y11 = Ylm_table(1, 1, theta, phi)
    @test isapprox(Y11, -sqrt(3/(8*pi)) * sin(theta) * exp(im*phi), atol=atol)

    Y1m1 = Ylm_table(1, -1, theta, phi)
    @test isapprox(Y1m1, sqrt(3/(8*pi)) * sin(theta) * exp(-im*phi), atol=atol)

    # l = 2
    Y20 = Ylm_table(2, 0, theta, phi)
    @test isapprox(Y20, sqrt(5/(16*pi)) * (3*cos(theta)^2 - 1) + im*0.0, atol=atol)

    Y21 = Ylm_table(2, 1, theta, phi)
    @test isapprox(Y21, -sqrt(15/(8*pi)) * sin(theta)*cos(theta) * exp(im*phi), atol=atol)

    Y2m1 = Ylm_table(2, -1, theta, phi)
    @test isapprox(Y2m1, sqrt(15/(8*pi)) * sin(theta)*cos(theta) * exp(-im*phi), atol=atol)

    Y22 = Ylm_table(2, 2, theta, phi)
    @test isapprox(Y22, sqrt(15/(32*pi)) * sin(theta)^2 * exp(2im*phi), atol=atol)

    Y2m2 = Ylm_table(2, -2, theta, phi)
    @test isapprox(Y2m2, sqrt(15/(32*pi)) * sin(theta)^2 * exp(-2im*phi), atol=atol)

    # l = 3
    Y30 = Ylm_table(3, 0, theta, phi)
    @test isapprox(Y30, sqrt(7/(16*pi)) * (5*cos(theta)^3 - 3*cos(theta)) + im*0.0, atol=atol)

    Y31 = Ylm_table(3, 1, theta, phi)
    @test isapprox(Y31, -sqrt(21/(64*pi)) * sin(theta) * (5*cos(theta)^2 - 1) * exp(im*phi), atol=atol)

    Y3m1 = Ylm_table(3, -1, theta, phi)
    @test isapprox(Y3m1, sqrt(21/(64*pi)) * sin(theta) * (5*cos(theta)^2 - 1) * exp(-im*phi), atol=atol)

    Y32 = Ylm_table(3, 2, theta, phi)
    @test isapprox(Y32, sqrt(105/(32*pi)) * sin(theta)^2 * cos(theta) * exp(2im*phi), atol=atol)

    Y3m2 = Ylm_table(3, -2, theta, phi)
    @test isapprox(Y3m2, sqrt(105/(32*pi)) * sin(theta)^2 * cos(theta) * exp(-2im*phi), atol=atol)

    Y33 = Ylm_table(3, 3, theta, phi)
    @test isapprox(Y33, -sqrt(35/(64*pi)) * sin(theta)^3 * exp(3im*phi), atol=atol)

    Y3m3 = Ylm_table(3, -3, theta, phi)
    @test isapprox(Y3m3, sqrt(35/(64*pi)) * sin(theta)^3 * exp(-3im*phi), atol=atol)
end