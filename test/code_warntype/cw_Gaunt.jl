include("../../src/utils/Gaunt.jl")


function check_Gaunt_code_warntype()
    f = zeros(Float64, 40)

    @code_warntype s3j(f,0,0,0,0,0,0)
    @code_warntype Clebsch_Gordan(f,0,0,0,0,0,0)
    @code_warntype Gaunt(f,0,0,0,0,0,0)
end

check_Gaunt_code_warntype()