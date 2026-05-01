include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Set_ProExpn! Si" begin
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
    Total_Hsize = system_grid.Total_Hsize

    HVNA = zeros(Float64, Total_Hsize)
    Set_ProExpn!( HVNA, [pao], [pspot], system_grid )

    HVNA_ref1 = [-0.4789323168738972, -0.13516551891366135, -2.9810172090602664e-11, -1.4055759052388775e-17, -1.770549344712422e-17, -2.7505794525366304e-11, -3.0940025917661047e-19, -1.3681415466436492e-19, 8.43314357268214e-18, -1.1403131422060203e-17]
    HVNA_ref2 = [-9.329159928162941e-7, 8.232268665494124e-7, -2.3592309784017275e-6, 2.635557507602936e-6, -2.4978691401904064e-6, 7.783171511970438e-6, 3.502298131642156e-7, -4.5907942088318304e-8, -3.6188856485647244e-7, 8.357743265506609e-7, -4.594737118997714e-7]

    @test isapprox(HVNA[1:10], HVNA_ref1, atol=1e-6)
    @test isapprox(HVNA[end-10:end], HVNA_ref2, atol=1e-6)
end

