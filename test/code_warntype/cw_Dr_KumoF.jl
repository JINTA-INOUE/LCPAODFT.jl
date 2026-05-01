include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_Dr_KumoF_code_warntype()

    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Atomic_PCC = pspot.Spe_Atomic_PCC
    Spe_PAO_RV = pao.Spe_PAO_RV
    
    @show typeof(Spe_Num_Mesh_VPS)
    @code_warntype Dr_KumoF(Spe_Num_Mesh_VPS, log(Spe_PAO_RV[begin]), Spe_VPS_RV, Spe_Atomic_PCC)
end

check_Dr_KumoF_code_warntype()