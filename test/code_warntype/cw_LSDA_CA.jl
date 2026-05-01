include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_LSDA_CA_code_warntype()
    den0 = 0.1
    den0 = 0.2
    @code_warntype LSDA_CA(den0, den1, 0)
end

check_LDA_CA_code_warntype()