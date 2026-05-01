include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Set_Nonlocal! Si without SOC" begin
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

    HNL = Vector{Vector{Float64}}(undef, 1)
	HNL[1] = zeros(Float64, Total_Hsize)
	iHNL = nothing
    Set_Nonlocal!("off", HNL, iHNL, [pao], [pspot], system_grid)

    HNL_ref1 = [0.23591205823917796, 0.014091656737249459, 3.04392146038334e-11, 6.694052382236842e-19, 4.2050608993518914e-19, 3.5774325198221916e-11, 1.1628699516762757e-18, 1.1102484638362625e-18, -5.286241029695551e-21, 3.475202363278469e-18]
    HNL_ref2 = [-5.682448603107019e-7, -2.6045178336656657e-7, -5.389483083407372e-7, 1.488962119459827e-6, 6.850414787797347e-7, 1.630791023974417e-6, 4.854036066332703e-8, 1.6641694032547902e-7, 6.828158225344807e-8, 2.331452135919867e-7, 1.094796111303171e-7]
    iHNL_ref = nothing

    @test isapprox(HNL[1][1:10], HNL_ref1, atol=1e-6)
    @test isapprox(HNL[1][end-10:end], HNL_ref2, atol=1e-6)
    @test iHNL == iHNL_ref
end

