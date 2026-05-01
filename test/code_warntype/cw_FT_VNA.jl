include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_FT_VNA_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    
    Spe_CrudeVNA_Bessel= zeros(Float64, GL_Mesh)
    @code_warntype FT_VNA!(pao, pspot, Spe_CrudeVNA_Bessel)
end

check_FT_VNA_code_warntype()