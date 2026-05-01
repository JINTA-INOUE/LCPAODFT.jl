include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "dYlm_complex" begin

    atol = 1e-11
    theta = rand()
    phi = rand()
    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)

    Lmax_Four_Int = 3
    fact = zeros(Float64, 2*Lmax_Four_Int+1, 2*Lmax_Four_Int+1)
    for i = 0:2*Lmax_Four_Int, j = 0:2*Lmax_Four_Int
        tmp0 = sqrt(factorial(big(i)))
        tmp1 = sqrt(factorial(big(j)))
        fact[i+1,j+1] = tmp0/tmp1
    end

    getY(SH) = SH[1] + im*SH[2]
    getd(dSH) = dSH[1] + im*dSH[2]


    
    l = 0; m = 0
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    @test isapprox(Y, 0.5*sqrt(1/pi) + 0.0*im, atol=atol)
    @test isapprox(dYt, 0.0 + 0.0*im, atol=atol)
    @test isapprox(dYp, 0.0 + 0.0*im, atol=atol)

    l = 1; m = 0
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    @test isapprox(Y, sqrt(3/(4*pi))*cos(theta) + 0.0*im, atol=atol)
    @test isapprox(dYt, -sqrt(3/(4*pi))*sin(theta) + 0.0*im, atol=atol)
    @test isapprox(dYp, 0.0 + 0.0*im, atol=atol)

    l = 1; m = 1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = -sqrt(3/(8*pi))*sin(theta)*exp(im*phi)
    dYt_expected = -sqrt(3/(8*pi))*cos(theta)*exp(im*phi)
    dYp_expected = im*1 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 1; m = -1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(3/(8*pi))*sin(theta)*exp(-im*phi)
    dYt_expected = sqrt(3/(8*pi))*cos(theta)*exp(-im*phi)
    dYp_expected = im*(-1) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)


    l = 2; m = 0
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(5/(16*pi))*(3*cos(theta)^2 - 1)
    dYt_expected = sqrt(5/(16*pi))*( -6*cos(theta)*sin(theta) )
    dYp_expected = 0.0*im

    @test isapprox(Y, Y_expected + 0.0*im, atol=atol)
    @test isapprox(dYt, dYt_expected + 0.0*im, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 2; m = 2
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(15/(32*pi))*sin(theta)^2 * exp(2.0*im*phi)
    dYt_expected = sqrt(15/(32*pi))*(2*sin(theta)*cos(theta)) * exp(2.0*im*phi)
    dYp_expected = im*2 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 2; m = -2
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(15/(32*pi)) * sin(theta)^2 * exp(-2.0*im*phi)
    dYt_expected = sqrt(15/(32*pi)) * (2*sin(theta)*cos(theta)) * exp(-2.0*im*phi)
    dYp_expected = im*(-2) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 2; m = 1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = -sqrt(15/(8*pi))*sin(theta)*cos(theta)*exp(im*phi)
    dYt_expected = -sqrt(15/(8*pi))*(cos(theta)^2 - sin(theta)^2)*exp(im*phi)
    dYp_expected = im*1 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 2; m = -1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(15/(8*pi))*sin(theta)*cos(theta)*exp(-im*phi)
    dYt_expected = sqrt(15/(8*pi))*(cos(theta)^2 - sin(theta)^2)*exp(-im*phi)
    dYp_expected = im*(-1) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)


    l = 3; m = 0
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(7/(16*pi))*(5*cos(theta)^3 - 3*cos(theta))
    dYt_expected = sqrt(7/(16*pi))*( -15*cos(theta)^2*sin(theta) + 3*sin(theta) )
    dYp_expected = 0.0*im

    @test isapprox(Y, Y_expected + 0.0*im, atol=atol)
    @test isapprox(dYt, dYt_expected + 0.0*im, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = 3
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = -sqrt(35/(64*pi))*sin(theta)^3 * exp(3.0*im*phi)
    dYt_expected = -sqrt(35/(64*pi))*(3*sin(theta)^2*cos(theta)) * exp(3.0*im*phi)
    dYp_expected = im*3 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = -3
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(35/(64*pi)) * sin(theta)^3 * exp(-3.0*im*phi)
    dYt_expected = sqrt(35/(64*pi)) * (3*sin(theta)^2 * cos(theta)) * exp(-3.0*im*phi)
    dYp_expected = im*(-3) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = 1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = -sqrt(21/(64*pi)) * sin(theta) * (5*cos(theta)^2 - 1) * exp(im*phi)
    dYt_expected = -sqrt(21/(64*pi)) * (cos(theta)*(5*cos(theta)^2 - 1) - 10*sin(theta)^2*cos(theta)) * exp(im*phi)
    dYp_expected = im*1 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = -1
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(21/(64*pi)) * sin(theta) * (5*cos(theta)^2 - 1) * exp(-im*phi)
    dYt_expected = sqrt(21/(64*pi)) * (cos(theta)*(5*cos(theta)^2 - 1) - 10*sin(theta)^2*cos(theta)) * exp(-im*phi)
    dYp_expected = im*(-1) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = 2
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(105/(32*pi)) * sin(theta)^2 * cos(theta) * exp(2.0*im*phi)
    dYt_expected = sqrt(105/(32*pi)) * (2*sin(theta)*cos(theta)^2 - sin(theta)^3) * exp(2.0*im*phi)
    dYp_expected = im*2 * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)

    l = 3; m = -2
    calc_Ylm!(l,m,fact[l-abs(m)+1,l+abs(m)+1],theta,phi,SH,dSHt,dSHp)
    Y = getY(SH)
    dYt = getd(dSHt)
    dYp = getd(dSHp)

    Y_expected = sqrt(105/(32*pi)) * sin(theta)^2 * cos(theta) * exp(-2.0*im*phi)
    dYt_expected = sqrt(105/(32*pi)) * (2*sin(theta)*cos(theta)^2 - sin(theta)^3) * exp(-2.0*im*phi)
    dYp_expected = im*(-2) * Y_expected

    @test isapprox(Y, Y_expected, atol=atol)
    @test isapprox(dYt, dYt_expected, atol=atol)
    @test isapprox(dYp, dYp_expected, atol=atol)
end