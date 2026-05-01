include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Set_VNA2 Si" begin
    verbosity = 0
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]
    Natom = 2
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    Total_NumOrbs = [13, 13]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs )
    system_grid = ucell.system_grid

    HVNA2 = Set_VNA2([pao], [pspot], system_grid)

    HVNA2_ref1 = [-0.40919332350907706, -0.19627843948785825, -2.981069352412969e-11, 0.0, -1.8253785202344073e-27, -2.7504361624176647e-11, 0.0, -1.6841564212819613e-27, 7.257399316179395e-22, -1.2570184346438333e-21]
    HVNA2_ref2 = [7.362183205178302e-12, -7.787158003236868e-12, 1.4936853809385888e-11, -2.5136027538423742e-11, 2.6599086002885264e-11, -5.100358430907824e-11, 1.948148968679879e-12, -6.269223827405519e-14, -1.1665591570048568e-12, 2.2390799565986307e-12, -2.3661058060269744e-12]

    @test isapprox(HVNA2[1:10], HVNA2_ref1, atol=1e-6)
    @test isapprox(HVNA2[end-10:end], HVNA2_ref2, atol=1e-6)
end

