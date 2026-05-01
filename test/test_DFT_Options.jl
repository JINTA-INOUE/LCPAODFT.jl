include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "DFT_Options default_DFT_Options" begin
    
    dft_options = default_DFT_Options()

    @test dft_options.Mixing_method == "RMM-DIISH"
    @test dft_options.SCF_criterion == 1e-6
    @test dft_options.SCF_max == 10
    @test dft_options.Mixing_weight == 0.3
    @test dft_options.Min_Mixing_weight == 0.001 
    @test dft_options.Max_Mixing_weight == 0.4
    @test dft_options.Max_Mixing_weight2 == 0.4
    @test dft_options.Num_Mixing_Pulay == 5
    @test dft_options.SCF_RENZOKU == -1
    @test dft_options.Start_Pulay_SCF == 6
    @test dft_options.crystal_sym == false
    @test dft_options.time_rev == true
end



@testset "DFT_Options" begin
    
    Mixing_method = "RMM-DIISH"
    SCF_criterion = 1e-6
    SCF_max = 10
    Init_Mixing_weight = 0.3
    Min_Mixing_weight = 0.001
    Max_Mixing_weight = 0.4
    Num_Mixing_Pulay = 5
    Start_Pulay_SCF = 6
    time_rev = true

    dft_options = dft_options = DFT_Options( Mixing_method, SCF_criterion, SCF_max, 
    Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
    Num_Mixing_Pulay, -1, Start_Pulay_SCF, false, time_rev)


    @test dft_options.Mixing_method == "RMM-DIISH"
    @test dft_options.SCF_criterion == 1e-6
    @test dft_options.SCF_max == 10
    @test dft_options.Mixing_weight == 0.3
    @test dft_options.Min_Mixing_weight == 0.001 
    @test dft_options.Max_Mixing_weight == 0.4
    @test dft_options.Max_Mixing_weight2 == 0.4
    @test dft_options.Num_Mixing_Pulay == 5
    @test dft_options.SCF_RENZOKU == -1
    @test dft_options.Start_Pulay_SCF == 6
    @test dft_options.crystal_sym == false
    @test dft_options.time_rev == true
end



@testset "DFT_Options from DFT_Setup" begin
    
    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]*"AU"
    atomorb = ["Si7.0-s2p2d1"]
    atomsymbol = ["Si", "Si"]
    atompos = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]*"Frac"
    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 50
    xc_type = "LDA"
    kmesh = (7,7,7)
    verbose = false

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, verbose)
    
    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    Init_Mixing_weight = dft_setup.Init_Mixing_weight
    Min_Mixing_weight = dft_setup.Min_Mixing_weight
    Max_Mixing_weight = dft_setup.Max_Mixing_weight
    Num_Mixing_Pulay = dft_setup.Num_Mixing_Pulay
    Start_Pulay_SCF = dft_setup.Start_Pulay_SCF
    time_rev = dft_setup.time_rev

    dft_options = dft_options = DFT_Options( Mixing_method, SCF_criterion, SCF_max, 
    Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
    Num_Mixing_Pulay, -1, Start_Pulay_SCF, false, time_rev)


    @test dft_options.Mixing_method == "RMM-DIISH"
    @test dft_options.SCF_criterion == 1e-8
    @test dft_options.SCF_max == 50
    @test dft_options.Mixing_weight == 0.3
    @test dft_options.Min_Mixing_weight == 0.001 
    @test dft_options.Max_Mixing_weight == 0.4
    @test dft_options.Max_Mixing_weight2 == 0.4
    @test dft_options.Num_Mixing_Pulay == 5
    @test dft_options.SCF_RENZOKU == -1
    @test dft_options.Start_Pulay_SCF == 6
    @test dft_options.crystal_sym == false
    @test dft_options.time_rev == true
end