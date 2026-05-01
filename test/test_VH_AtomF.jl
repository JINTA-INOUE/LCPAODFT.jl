include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "VH_AtomF Si" begin
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)

    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Core_Charge = pspot.Spe_Core_Charge
    Spe_VH_Atom = Calc_Spe_VH_Atom(pao, pspot)

    N = 30
    rmin = 1e-3
    rmax = 10.0
    dr = (rmax-rmin)/N
    r = zeros(Float64, N)
    for i = 1:N
        r[i] = rmin + (i-1)*dr
    end
    check_VH_AtomF = zeros(Float64, N)

    for i = 1:N
        r2 = r[i]^2
        check_VH_AtomF[i] = VH_AtomF( Spe_Core_Charge, Spe_Num_Mesh_VPS, r2, Spe_VPS_RV, Spe_VH_Atom )
    end

    @test check_VH_AtomF ≈ [2.1560673457792108, 2.1527952335965774, 2.1378548281376077, 2.097070411463603, 2.0137960947435545, 1.8861356263494096, 1.7329706543234056, 1.5759005963454973, 1.4284466126628328, 1.2964521500103952, 1.1811006864512592, 1.0813235840821955, 0.9952183436295936, 0.9207479864529826, 0.8560434972929749, 0.7994985658339807, 0.7497761809090516, 0.7057788116267107, 0.6666110315451717, 0.6315409839260723, 0.5999699926688167, 0.5714041359069726, 0.5454347815978535, 0.5217232829247238, 0.49998746100905134, 0.4799903160445043, 0.46153144779901173, 0.4444395111471103, 0.42856831173307813, 0.4137917189100345]
end