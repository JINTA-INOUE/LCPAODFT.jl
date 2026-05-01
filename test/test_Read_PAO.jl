include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Read_PAO Si LDA" begin
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)

    pao_filebase = basename(pao.paofile)
    @test pao.Atom_symbol == "Si"
    @test pao.Atom_orbital == "s2p2d1"
    @test pao.Spe_MaxL_Basis == 2
    @test pao.Spe_Num_Basis == [2,2,1]
    @test pao.Spe_Total_NumOrbs == 13
    @test pao.Spe_Num_Mesh_PAO == 500
    @test pao.Spe_Atom_Cut1 == 7.0
    # @test pao.Spe_PAO_RV == 
    # @test pao.Spe_Atomic_Den == 
    @test pao.Spe_PAO_Lmax == 3
    @test pao.Spe_PAO_Mul == 15
    # @test pao.Spe_PAO_RWF == 
    # @test pao.Spe_RF_Bessel == 
    @test pao_filebase == "Si7.0.jld2"
end


@testset "Read_PAO Fe GGA-PBE" begin
    verbosity = 0
    pao = Read_PAO(7.0, "Fe", 6.0, "s2p2d1", "H"; verbosity)

    pao_filebase = basename(pao.paofile)
    @test pao.Atom_symbol == "Fe"
    @test pao.Atom_orbital == "s2p2d1"
    @test pao.Spe_MaxL_Basis == 2
    @test pao.Spe_Num_Basis == [2,2,1]
    @test pao.Spe_Total_NumOrbs == 13
    @test pao.Spe_Num_Mesh_PAO == 500
    @test pao.Spe_Atom_Cut1 == 6.0
    # @test pao.Spe_PAO_RV == 
    # @test pao.Spe_Atomic_Den == 
    @test pao.Spe_PAO_Lmax == 3
    @test pao.Spe_PAO_Mul == 15
    # @test pao.Spe_PAO_RWF == 
    # @test pao.Spe_RF_Bessel == 
    @test pao_filebase == "Fe6.0H.jld2"
end