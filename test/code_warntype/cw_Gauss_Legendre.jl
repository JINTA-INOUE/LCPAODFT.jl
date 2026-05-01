include("../../src/utils/Gauss_Legendre.jl")

function check_Gauss_Legendre_code_warntype()
    x = zeros(Float64, 3)
    w = zeros(Float64, 3)
    @code_warntype Gauss_Legendre(3)
    # @code_warntype Gauss_Legendre_x(3)
end

check_Gauss_Legendre_code_warntype()