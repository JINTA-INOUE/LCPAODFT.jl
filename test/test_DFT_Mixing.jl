include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "DFT_Mixing Si" begin
    
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
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    system_grid = ucell.system_grid

    dft_options = default_DFT_Options()
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay
    dft_mixing = DFT_Mixing( Nspin, dft_options, system_grid )
    SCF_max = dft_mixing.SCF_max

    @test dft_mixing.Mixing_method == "RMM-DIISH"
    @test dft_mixing.SCF_max == 10
    @test dft_mixing.is_convergence == false
    @test dft_mixing.Conv_iter == 0
    @test dft_mixing.Latvecs == Latvecs
    @test dft_mixing.Recvecs ≈ 2*pi*inv(Latvecs')
    @test dft_mixing.Natom == 2
    @test dft_mixing.Nspin == 1
    @test dft_mixing.FNAN == system_grid.FNAN
    @test dft_mixing.natn == system_grid.natn
    @test dft_mixing.Total_NumOrbs == system_grid.Total_NumOrbs
    @test dft_mixing.Total_Hsize == system_grid.Total_Hsize
    @test dft_mixing.ChemP == 0.0
    @test dft_mixing.NormRD[begin] == 100.0
    @test dft_mixing.NormRD[2:end] == zeros(Float64, SCF_max+1)
    @test dft_mixing.HisEele == zeros(Float64, SCF_max+2)
    @test size(dft_mixing.Density_xyz) == Ngrid
    @test length(dft_mixing.HisH) == Num_Mixing_Pulay+1
    @test size(dft_mixing.HisH[1]) == (system_grid.Total_Hsize, Nspin)
    @test length(dft_mixing.ResH) == Num_Mixing_Pulay+1
    @test size(dft_mixing.ResH[1]) == (system_grid.Total_Hsize, Nspin)
    @test dft_mixing.Ngrid == system_grid.Ngrid
end