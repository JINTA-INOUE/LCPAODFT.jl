include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "XC_PBE Si LDA" begin
    verbosity = 0
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Natom = 2
    Nspin = 1
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "off"
    xc_type = "LDA"
    Init_Atoms_Nspin = [[2.0,2.0],[2.0,2.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", xc_type, false; verbosity)
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    _, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid( SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, [pao], [pspot], ucell )
   
    xc_func = XC_Func( xc_type, SpinPol, Nspin, Ngrid )

    Set_XC_Grid!( xc_func, PCCDensity_Grid, Density_Grid )

    @test xc_func.xc_type == xc_type
    @test xc_func.SpinPol == SpinPol
    @test xc_func.Nspin == Nspin
    @test xc_func.Ngrid == Ngrid
    @test xc_func.Vxc_Grid ≈ [[-0.4446008699668833, -0.36738544854127303, -0.36738544854127253, -0.36738544854127253, -0.48556124412337187, -0.36738544854127264, -0.36738544854127264, -0.36738544854127253, -0.20093057937905404, -0.36738544854127253, -0.48556124412337187, -0.36738544854127253, -0.4855612441233718, -0.4855612441233718, -0.24276599919949265, -0.36738544854127253, -0.24276599919949265, -0.24276599919949265, -0.36738544854127264, -0.36738544854127253, -0.20093057937905404, -0.36738544854127253, -0.24276599919949265, -0.24276599919949265, -0.20093057937905404, -0.24276599919949265, -0.20093057937905404]]
end


@testset "XC_PBE Si GGA-PBE" begin

    verbosity = 0
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Natom = 2
    Nspin = 1
    atom2spe = [1, 1]
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Atom_Cut1 = [7.0, 7.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "off"
    xc_type = "GGA-PBE"
    Init_Atoms_Nspin = [[2.0,2.0],[2.0,2.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]]
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", xc_type, false; verbosity)
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    _, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid( SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, [pao], [pspot], ucell )
   
    xc_func = XC_Func( xc_type, SpinPol, Nspin, Ngrid )

    Set_XC_Grid!( xc_func, PCCDensity_Grid, Density_Grid )
    
    # @show xc_func.ex2primitive
    # @show xc_func.dDensity_Grid
    # @show xc_func.Vxc_Grid
    # @show xc_func.dEXC_dGD
    # @show xc_func.igtv


    @test xc_func.xc_type == xc_type
    @test xc_func.SpinPol == SpinPol
    @test xc_func.Nspin == Nspin
    @test xc_func.Ngrid == Ngrid
    # @test xc_func.ex2primitive == 
    # @test xc_func.dDensity_Grid ≈
    # @test xc_func.dEXC_dGD ≈
    # @test xc_func.Vxc_Grid ≈
end


@testset "XC_PBE BCC_Fe LSDA" begin
    verbosity = 0
    Latvecs = [-2.71176  2.71176  2.71176; 2.71176 -2.71176  2.71176; 2.71176  2.71176 -2.71176]
    Natom = 1
    Nspin = 2
    atom2spe = [1]
    Gxyz = [[0.0,0.0,0.0]]
    Atom_Cut1 = [6.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "on"
    xc_type = "LSDA"
    Init_Atoms_Nspin = [[8.0,6.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0]]
    pao = Read_PAO(16.0, "Fe", 6.0, "s2p2d1", "H"; verbosity)
    pspot = Read_VPS("Fe", "H", xc_type, false; verbosity)
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    _, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid( SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, [pao], [pspot], ucell )
   
    xc_func = XC_Func( xc_type, SpinPol, Nspin, Ngrid )

    Set_XC_Grid!( xc_func, PCCDensity_Grid, Density_Grid )

    @test xc_func.xc_type == xc_type
    @test xc_func.SpinPol == SpinPol
    @test xc_func.Nspin == Nspin
    @test xc_func.Ngrid == Ngrid
    @test xc_func.Vxc_Grid ≈ [[-1.2149312126522975, -0.5723999028367164, -0.5723999028367156, -0.5723999028367164, -0.4760876844448474, -0.3795360488809768, -0.5723999028367156, -0.3795360488809768, -0.4760876844448473, -0.5723999028367164, -0.4760876844448474, -0.3795360488809768, -0.4760876844448474, -0.5723999028367164, -0.3795360488809768, -0.3795360488809768, -0.37953604888097675, -0.3795360488809768, -0.5723999028367156, -0.3795360488809768, -0.4760876844448473, -0.3795360488809768, -0.37953604888097675, -0.3795360488809768, -0.4760876844448473, -0.3795360488809768, -0.5723999028367156], [-1.2026869942420644, -0.5378501447465774, -0.5378501447465768, -0.5378501447465774, -0.4490198444695973, -0.3596772073151438, -0.5378501447465768, -0.3596772073151438, -0.4490198444695972, -0.5378501447465774, -0.4490198444695973, -0.3596772073151438, -0.4490198444695973, -0.5378501447465774, -0.3596772073151438, -0.3596772073151438, -0.35967720731514374, -0.3596772073151438, -0.5378501447465768, -0.3596772073151438, -0.4490198444695972, -0.3596772073151438, -0.35967720731514374, -0.3596772073151438, -0.4490198444695972, -0.3596772073151438, -0.5378501447465768]]
end


@testset "XC_PBE BCC_Fe GGA-PBE" begin
    verbosity = 0
    Latvecs = [-2.71176  2.71176  2.71176; 2.71176 -2.71176  2.71176; 2.71176  2.71176 -2.71176]
    Natom = 1
    Nspin = 2
    atom2spe = [1]
    Gxyz = [[0.0,0.0,0.0]]
    Atom_Cut1 = [6.0]
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    SpinPol = "on"
    xc_type = "GGA-PBE"
    Init_Atoms_Nspin = [[8.0,6.0]]
    Init_Atoms_Angle = [[0.0,0.0,0.0,0.0]]
    pao = Read_PAO(16.0, "Fe", 6.0, "s2p2d1", "H"; verbosity)
    pspot = Read_VPS("Fe", "H", xc_type, false; verbosity)
    ucell = UCell( Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin )
    _, PCCDensity_Grid, Density_Grid = Set_AdenPCC_Grid( SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, [pao], [pspot], ucell )
   
    xc_func = XC_Func( xc_type, SpinPol, Nspin, Ngrid )

    Set_XC_Grid!( xc_func, PCCDensity_Grid, Density_Grid )
    
    # @show xc_func.ex2primitive
    # @show xc_func.dDensity_Grid
    # @show xc_func.Vxc_Grid
    # @show xc_func.dEXC_dGD
    # @show xc_func.igtv


    @test xc_func.xc_type == xc_type
    @test xc_func.SpinPol == SpinPol
    @test xc_func.Nspin == Nspin
    @test xc_func.Ngrid == Ngrid
    # @test xc_func.ex2primitive == 
    # @test xc_func.dDensity_Grid ≈
    # @test xc_func.dEXC_dGD ≈
    # @test xc_func.Vxc_Grid ≈
end