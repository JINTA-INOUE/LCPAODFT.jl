include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_VNAF_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Vna = Calc_Spe_Vna(pao, pspot)
    
    @show typeof(Spe_Num_Mesh_VPS)
    @code_warntype VNAF(Spe_Num_Mesh_VPS, 1.0, Spe_VPS_RV, Spe_Vna)
end

check_VNAF_code_warntype()