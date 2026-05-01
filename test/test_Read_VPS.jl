include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Read_VPS Si LDA" begin
    verbosity = 0
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)

    pspot_filebase = basename(pspot.psfile)
    @test pspot.Atom_symbol == "Si"
    @test pspot.Atom_extra == ""
    @test pspot.xc_type == "LDA"
    @test pspot.Spe_Core_Charge == 4.0
    @test pspot.VPS_j_dependency == 0
    @test pspot.SO_switch == false
    @test pspot.Spe_Num_Mesh_VPS == 500
    @test pspot.Spe_Num_RVPS == 6
    @test pspot.Spe_VPS_List == [0, 0, 1, 1, 2, 2]
    @test pspot.Spe_VNLE == [1.0 1.0 1.0 1.0 -1.0 1.0; 1.0 1.0 1.0 1.0 -1.0 1.0]
    # @test pspot.Spe_VNL
    # @test pspot.Spe_VPS_RV
    # @test pspot.Spe_Vcore
    @test pspot.is_pcc == true
    # @test pspot.Spe_Atomic_PCC
    @test pspot_filebase == "Si_CA19.jld2"
end