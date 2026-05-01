include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_RadialF_code_warntype()
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", "")
    Spe_Num_Mesh_PAO = pao.Spe_Num_Mesh_PAO
    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_PAO_RWF = pao.Spe_PAO_RWF
    
    @show typeof(Spe_Num_Mesh_PAO)
    @code_warntype RadialF(0, Spe_Num_Mesh_PAO, 1.0, Spe_PAO_RV, Spe_PAO_RWF[1][1])
end

check_RadialF_code_warntype()