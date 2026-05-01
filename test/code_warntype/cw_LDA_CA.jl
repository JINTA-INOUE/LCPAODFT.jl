include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_LDA_CA_code_warntype()
    den = 1.0
    @code_warntype LDA_CA(den, 0)
end

check_LDA_CA_code_warntype()