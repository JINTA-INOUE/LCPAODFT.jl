include("../../src/utils/Calc_Ngrid.jl")

function check_Calc_Ngrid_code_warntype()
    Ecut = 150.0
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    @code_warntype Calc_Ngrid(Ecut, Latvecs)
end

check_Calc_Ngrid_code_warntype()