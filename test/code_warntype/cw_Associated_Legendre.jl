include("../../src/utils/Associated_Legendre.jl")


function check_Associated_Legendre_code_warntype()
    l = 1
    m = 1
    x = 0.3

    @code_warntype Associated_Legendre(l, m, x)
    @code_warntype Associated_Legendre2(l, m, x)
end

check_Associated_Legendre_code_warntype()