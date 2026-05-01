include("../../src/LCPAODFT.jl")
using .LCPAODFT

function check_Spe_VHart_Atom_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Atomic_Den = pao.Spe_Atomic_Den
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Vcore = pspot.Spe_Vcore
    Spe_Core_Charge = pspot.Spe_Core_Charge

    Spe_VH_Atom = zeros(Float64, Spe_Num_Mesh_VPS+2)

    @code_warntype Spe_VHart_Atom!(Spe_VPS_RV, Spe_PAO_RV, Spe_VH_Atom, Spe_Atomic_Den)
end

check_Spe_VHart_Atom_code_warntype()