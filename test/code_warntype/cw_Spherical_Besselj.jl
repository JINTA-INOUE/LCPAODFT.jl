include("../../src/utils/Spherical_Besselj.jl")
using Bessels: besselj

function check_Spherical_Besselj_code_warntype()
    x = rand(Float64, 100)
    @code_warntype SphericalBesselj(0, x)
end

function check_dSpherical_Besselj_code_warntype()
    lmax = 3
    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    x = rand(Float64, 100)
    SphB = zeros(Float64, lmax+3)
    dSphB = zeros(Float64, lmax+3)

    @code_warntype Calc_SphericalBesselj!(lmax, x, tsb, SphB, dSphB)
    @code_warntype Calc_SphericalBesselj2!(lmax, x, tsb, SphB)
    @code_warntype Calc_SphericalBesselj2!(lmax, x, tsb, SphB, dSphB)
end

# check_Spherical_Besselj_code_warntype()
check_dSpherical_Besselj_code_warntype()