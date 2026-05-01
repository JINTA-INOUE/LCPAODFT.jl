include("../../src/LCPAODFT.jl")
using .LCPAODFT

function check_Ylm_complex_code_warntype()
    theta = rand()
    phi = rand()
    @code_warntype Ylm_table(0,0,theta,phi)
end


function check_Ylm_real_code_warntype()
    theta = rand()
    phi = rand()
    @code_warntype Ylm_real(0,0,theta,phi)
end


function check_calc_Ylm!_code_warntype()
    theta = rand()
    phi = rand()
    fact = rand()
    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp  = zeros(Float64, 2)
    @code_warntype calc_Ylm!(0, 0, fact, theta, phi, SH)
    @code_warntype calc_Ylm!(0, 0, fact, theta, phi, SH, dSHt, dSHp)
end

# check_Ylm_complex_code_warntype()
# check_Ylm_real_code_warntype()
# check_calc_Ylm!_code_warntype()