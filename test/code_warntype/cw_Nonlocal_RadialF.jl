include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_Nonlocal_RadialF_code_warntype()
    verbosity = 0
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_VNL = pspot.Spe_VNL

    @show typeof(Spe_Num_Mesh_VPS)
    @code_warntype Nonlocal_RadialF(1, Spe_Num_Mesh_VPS, 1.0, Spe_VPS_RV, Spe_VNL[1][1])
end

check_Nonlocal_RadialF_code_warntype()