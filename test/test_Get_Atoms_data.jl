include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Get_Atoms_data" begin

    Atom_orb = ["Si7.0-s2p2d1"]
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atom_orb)

    @test Spe_Symbol == ["Si"]
    @test Spe_cutoff == [7.0]
    @test Spe_orb == ["s2p2d1"]
    @test Spe_extra == [""]


    Atom_orb = ["Si7.0-s2p2d1", "Si7.0-s2p2d1"]
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atom_orb)

    @test Spe_Symbol == ["Si", "Si"]
    @test Spe_cutoff == [7.0, 7.0]
    @test Spe_orb == ["s2p2d1", "s2p2d1"]
    @test Spe_extra == ["", ""]


    Atom_orb = ["Bi8.0-s3p2d2f1", "Se7.0-s3p2d2"]
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atom_orb)

    @test Spe_Symbol == ["Bi", "Se"]
    @test Spe_cutoff == [8.0, 7.0]
    @test Spe_orb == ["s3p2d2f1", "s3p2d2"]
    @test Spe_extra == ["", ""]


    Atom_orb = ["Se7.0-s3p2d2", "Bi8.0-s3p2d2f1"]
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atom_orb)

    @test Spe_Symbol == ["Se", "Bi"]
    @test Spe_cutoff == [7.0, 8.0]
    @test Spe_orb == ["s3p2d2", "s3p2d2f1"]
    @test Spe_extra == ["", ""]
end