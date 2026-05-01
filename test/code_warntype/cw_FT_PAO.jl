include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_FT_PAO_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_MaxL_Basis = pao.Spe_MaxL_Basis
    Spe_Num_Basis = pao.Spe_Num_Basis
    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_PAO_RWF = pao.Spe_PAO_RWF
    Spe_RF_Bessel = pao.Spe_RF_Bessel
    
    @code_warntype FT_PAO!(Spe_Atom_Cut1, Spe_MaxL_Basis, Spe_Num_Basis, Spe_PAO_RV, Spe_PAO_RWF, Spe_RF_Bessel)
end

check_FT_PAO_code_warntype()