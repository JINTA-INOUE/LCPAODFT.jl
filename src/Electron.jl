abstract type AbstractBloch end
abstract type ClusterBloch <: AbstractBloch end
abstract type CrystalBloch <: AbstractBloch end


mutable struct Cluster_Electron <: ClusterBloch
    system::String
    SpinPol::String
    Spindeg::Int32
    spinsize::Int32
    fsize::Int32
    Nfsize::Int32
    E_Temp::Float64
    Core_Charge::Vector{Float64}
    TotalZ::Float64
    Nocc::Float64
    Eele::Float64
    ChemP::Float64
    FF::Array{Float64,2}
    Enk::Array{Float64,2}
    Cnk::Vector{Matrix{ComplexF64}}
end


mutable struct Crystal_Electron <: CrystalBloch
    system::String
    SpinPol::String
    Spindeg::Int32
    spinsize::Int32
    fsize::Int32
    Nfsize::Int32
    Nkpt::Int32
    MPI_Nkpt::Int32
    E_Temp::Float64
    Core_Charge::Vector{Float64}
    TotalZ::Float64
    Nocc::Float64
    Eele::Float64
    ChemP::Float64
    FF::Array{Float64,3}
    Enk::Array{Float64,3}
    Cnk::Vector{Vector{Matrix{ComplexF64}}}
end



"""
    Electron(...)

Create an instance of `Electron`.

Mandatory arguments:

- `kpoints`: an instance of `KPoints`
- `SpinPol`: which use (`no`, `on`, `nc`)
- `E_Temp`: Electron Temperatue
- `Core_Charge`: Atomic Core charge
- `fsize`: the number of total Orbitals in unit cell
- `system`: system name (`Cluster`, `Crystal`)
"""
function Electron(kpoints::KPoints, SpinPol::String, E_Temp, Core_Charge, fsize::Integer, system::String)

    TotalZ = sum(Core_Charge)

    if SpinPol == "off"
        spinsize = 1
        Spindeg = 2
        Nocc = TotalZ/2
        Nfsize = fsize
    elseif SpinPol == "on"
        spinsize = 2
        Spindeg = 1
        Nocc = TotalZ/2
        Nfsize = fsize
    elseif SpinPol == "nc"
        spinsize = 1
        Spindeg = 1
        Nocc = TotalZ
        Nfsize = 2*fsize
    else
        println("Now SpinPol is $SpinPol")
        error("please check SpinPol")
    end


    if system == "Cluster"

        FF = zeros(Float64, Nfsize, spinsize)
        Enk = zeros(Float64, Nfsize, spinsize)
        Cnk = Vector{Matrix{ComplexF64}}(undef, spinsize)
        for spin = 1:spinsize
            Cnk[spin] = zeros(ComplexF64, Nfsize, Nfsize)
        end
        
        return Cluster_Electron(system, SpinPol, Spindeg, 
                                spinsize, fsize, Nfsize,
                                E_Temp, Core_Charge, TotalZ, Nocc, 0.0, 0.0, 
                                FF, Enk, Cnk)
    elseif system == "Crystal"

        Nkpt = kpoints.Nkpt
        MPI_Nkpt = kpoints.MPI_Nkpt

        FF = zeros(Float64, Nfsize, Nkpt, spinsize)
        Enk = zeros(Float64, Nfsize, Nkpt, spinsize)
        Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
        for spin = 1:spinsize
            Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
            for ik = 1:MPI_Nkpt
                Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
            end
        end
        
        return Crystal_Electron(system, SpinPol, Spindeg, 
                                spinsize, fsize, Nfsize, Nkpt, MPI_Nkpt,
                                E_Temp, Core_Charge, TotalZ, Nocc, 0.0, 0.0, 
                                FF, Enk, Cnk)
    end
end
