include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "KPoints using time_rev = true" begin
    system = "Crystal"
    SpinPol = "off"
    kmesh = (3,3,3)
    time_rev = true
    kpoints = KPoints( system, SpinPol, kmesh, time_rev )

    @test kpoints.AllNkpt == 27
    @test kpoints.Nkpt == 14
    # @test kpoints.MPI_Nkpt == 
    @test kpoints.kmesh == (3,3,3)
    # @test kpoints.MPI_kpts == 
    # @test kpoints.MPI_kweight == 
    @test kpoints.All_kweight == [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]
    # @test kpoints.MPI_krange == 
    # @test kpoints.MPkpts == 
    @test kpoints.crystal_sym == false
    @test kpoints.time_rev == true
end