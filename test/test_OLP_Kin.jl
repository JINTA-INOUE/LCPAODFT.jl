include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Set_OLP_Kin Si" begin
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
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs )
    system_grid = ucell.system_grid
    Total_Hsize = system_grid.Total_Hsize

    OLP = zeros(Float64, Total_Hsize)
    Hkin = zeros(Float64, Total_Hsize)
    Set_OLP_Kin!( OLP, Hkin, [pao], system_grid )

    OLP_ref1 = [0.9999998316102195, 1.8650412318473225e-6, -3.4049549310010006e-11, 0.0, -2.0849335787456865e-27, 7.167927685776964e-12, 0.0, 4.389089848453227e-28, -9.271348434019752e-23, 1.6058446542396352e-22]
    OLP_ref2 = [7.142436118086677e-6, -6.400040090852871e-6, 1.3913674222556454e-5, -2.413173654896234e-5, 2.158296157657883e-5, -4.698908561173294e-5, -1.9272212067442396e-6, -1.1171803252993817e-7, 1.0382029941484499e-6, -2.243983037091807e-6, 2.0253141844785407e-6]
    Hkin_ref1 = [0.23977463666305046, 0.2501343338218208, -9.793986624578283e-12, 0.0, -5.997087185340892e-28, -4.906377658577917e-13, 0.0, -3.0042898474927655e-29, 4.0044564150559517e-22, -6.935921967572031e-22]
    Hkin_ref2 = [-3.584018754160301e-5, 3.0066261895509927e-5, -6.87934122601595e-5, 0.00011339151422733202, -9.439513367512058e-5, 0.00021728483817855834, 1.1427860614026186e-5, 9.863893649851518e-7, -5.940285491779441e-6, 1.3360155031036608e-5, -1.1446175800618182e-5]

    @test isapprox(OLP[1:10], OLP_ref1, atol=1e-6)
    @test isapprox(OLP[end-10:end], OLP_ref2, atol=1e-6)
    @test isapprox(Hkin[1:10], Hkin_ref1, atol=1e-6)
    @test isapprox(Hkin[end-10:end], Hkin_ref2, atol=1e-6)
end

