include("../../src/utils/dampingF.jl")

function check_dampingF_code_warntype()
    rcut = 10.0
    r = 1.0
    @code_warntype dampingF(rcut, r)
    @code_warntype deri_dampingF(rcut, r)
end

check_dampingF_code_warntype()