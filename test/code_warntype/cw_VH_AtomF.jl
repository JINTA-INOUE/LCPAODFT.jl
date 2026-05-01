include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_VH_AtomF_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)

    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Core_Charge = pspot.Spe_Core_Charge

    Spe_VH_Atom = Calc_Spe_VH_Atom(pao, pspot)

    @show typeof(Spe_Num_Mesh_VPS)
    r = 1.0
    @code_warntype VH_AtomF(Spe_Core_Charge, Spe_Num_Mesh_VPS, r, Spe_VPS_RV, Spe_VH_Atom)
end

check_VH_AtomF_code_warntype()