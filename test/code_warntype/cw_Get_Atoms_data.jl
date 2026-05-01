include("../../src/LCPAODFT.jl")
using .LCPAODFT

function check_Get_Atoms_data_code_warntype()

    Atom_orb = ["Si7.0-s2p2d1"]
    @code_warntype Get_Atoms_data(Atom_orb)
end

check_Get_Atoms_data_code_warntype()