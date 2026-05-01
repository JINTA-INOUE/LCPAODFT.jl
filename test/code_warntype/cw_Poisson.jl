include("../../src/LCPAODFT.jl")
using .LCPAODFT

function check_Poisson_code_warntype()
    verbosity = 0
    Nspin = 1
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; verbosity)
    system_grid = ucell.system_grid

    dft_options = default_DFT_Options()
    dft_mixing = DFT_Mixing(Nspin, dft_options, system_grid)

    ADensity_Grid = [0.007158940964450057, 0.016335944309335652, 0.016335944309335586, 0.016335944309335586, 0.03587097291278265, 0.01633594430933559, 0.01633594430933559, 0.016335944309335586, 0.002226387584455796, 0.016335944309335586, 0.03587097291278265, 0.016335944309335586, 0.03587097291278265, 0.03587097291278265, 0.004172169081670134, 0.016335944309335586, 0.004172169081670134, 0.004172169081670133, 0.01633594430933559, 0.016335944309335586, 0.002226387584455796, 0.016335944309335583, 0.004172169081670134, 0.004172169081670133, 0.002226387584455796, 0.004172169081670134, 0.002226387584455796]
    dVHart_Grid = zeros(Float64, prod(Ngrid))
    @code_warntype Solve_Poisson!(dft_mixing, ADensity_Grid, dVHart_Grid)
end


check_Poisson_code_warntype()