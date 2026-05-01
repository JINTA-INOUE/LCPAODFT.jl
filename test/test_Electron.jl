include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "Electron Si" begin
    system = "Crystal"
    SpinPol = "off"
    kmesh = (3,3,3)
    time_rev = true
    E_Temp = 300.0
    Atoms_Core_Charge = [4.0, 4.0]
    Total_NumOrbs = [13, 13]
    kpoints = KPoints(system, SpinPol, kmesh, time_rev)
    electron = Electron( kpoints, SpinPol, E_Temp, Atoms_Core_Charge, sum(Total_NumOrbs), system )

    @test electron.system == system
    @test electron.SpinPol == SpinPol
    @test electron.Spindeg == 2
    @test electron.kpoints == kpoints
    @test electron.Nspin == 1
    @test electron.E_Temp == E_Temp
    @test electron.Core_Charge == Atoms_Core_Charge
    @test electron.TotalZ == sum(Atoms_Core_Charge)
    @test electron.Nocc == 4.0
    @test electron.Eele == 0.0
    @test electron.ChemP == 0.0
    @test size(electron.FF) == (1, sum(Total_NumOrbs), 14)
    @test size(electron.Enk) == (1, sum(Total_NumOrbs), 14)
    @test length(electron.Cnk) == 1
    @test size(electron.Cnk[1][1]) == (sum(Total_NumOrbs), sum(Total_NumOrbs))
end


@testset "Electron Fe" begin
    system = "Crystal"
    SpinPol = "on"
    kmesh = (3,3,3)
    time_rev = true
    E_Temp = 300.0
    Atoms_Core_Charge = [14.0]
    Total_NumOrbs = [13]
    kpoints = KPoints(system, SpinPol, kmesh, time_rev)
    electron = Electron( kpoints, SpinPol, E_Temp, Atoms_Core_Charge, sum(Total_NumOrbs), system )

    @test electron.system == system
    @test electron.SpinPol == SpinPol
    @test electron.Spindeg == 1
    @test electron.kpoints == kpoints
    @test electron.Nspin == 2
    @test electron.E_Temp == E_Temp
    @test electron.Core_Charge == Atoms_Core_Charge
    @test electron.TotalZ == sum(Atoms_Core_Charge)
    @test electron.Nocc == 7.0
    @test electron.Eele == 0.0
    @test electron.ChemP == 0.0
    @test size(electron.FF) == (2, sum(Total_NumOrbs), 14)
    @test size(electron.Enk) == (2, sum(Total_NumOrbs), 14)
    @test length(electron.Cnk) == 2
    @test size(electron.Cnk[1][1]) == (sum(Total_NumOrbs), sum(Total_NumOrbs))
end


@testset "Electron GaAs with SOC" begin
    system = "Crystal"
    SpinPol = "nc"
    kmesh = (3,3,3)
    time_rev = false
    E_Temp = 300.0
    Atoms_Core_Charge = [13.0, 15.0]
    Total_NumOrbs = [13, 13]
    kpoints = KPoints(system, SpinPol, kmesh, time_rev)
    electron = Electron( kpoints, SpinPol, E_Temp, Atoms_Core_Charge, sum(Total_NumOrbs), system )

    @test electron.system == system
    @test electron.SpinPol == SpinPol
    @test electron.Spindeg == 1
    @test electron.kpoints == kpoints
    @test electron.Nspin == 1
    @test electron.E_Temp == E_Temp
    @test electron.Core_Charge == Atoms_Core_Charge
    @test electron.TotalZ == sum(Atoms_Core_Charge)
    @test electron.Nocc == 28.0
    @test electron.Eele == 0.0
    @test electron.ChemP == 0.0
    @test size(electron.FF) == (1, 2*sum(Total_NumOrbs), 27)
    @test size(electron.Enk) == (1, 2*sum(Total_NumOrbs), 27)
    @test length(electron.Cnk) == 1
    @test size(electron.Cnk[1][1]) == (2*sum(Total_NumOrbs), 2*sum(Total_NumOrbs))
end