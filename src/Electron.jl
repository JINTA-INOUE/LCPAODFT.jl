abstract type AbstractBloch end
abstract type ClusterBloch <: AbstractBloch end
abstract type CrystalBloch <: AbstractBloch end


mutable struct Crystal_Electron <: CrystalBloch
    system::String
    SpinPol::String
    Spindeg::Int32
    kpoints::KPoints
    Nspin::Int32
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
- `system`: system name (`Atom`, `Cluster`, `Crystal`)
"""
function Electron(kpoints::Union{KPoints, Nothing}, SpinPol::String, E_Temp, Core_Charge, fsize::Integer, system::String)

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        println("Now SpinPol is $SpinPol")
        error("please check SpinPol")
    end

    if SpinPol == "off"
        Spindeg = 2
    elseif SpinPol ∈ ("on", "nc")
        Spindeg = 1
    else
        error("please check SpinPol")
    end

    TotalZ = sum(Core_Charge)
    if SpinPol == "nc"
        Nocc = TotalZ
    else
        Nocc = TotalZ/2
    end



    if system ∈ ("Atom", "Cluster")
        error("not support $system.")
    elseif system == "Crystal"

        Nkpt = kpoints.Nkpt
        MPI_Nkpt = kpoints.MPI_Nkpt
        Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
        
        if SpinPol ∈ ("off", "on")
            FF = zeros(Float64, spinsize, fsize, Nkpt)
            Enk = zeros(Float64, spinsize, fsize, Nkpt)
            for spin = 1:spinsize
                Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
                for ik = 1:MPI_Nkpt
                    Cnk[spin][ik] = zeros(ComplexF64, fsize, fsize)
                end
            end
        elseif SpinPol == "nc"
            FF = zeros(Float64, 1, 2*fsize, Nkpt)
            Enk = zeros(Float64, 1, 2*fsize, Nkpt)
            Cnk[1] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
            for ik = 1:MPI_Nkpt
                Cnk[1][ik] = zeros(ComplexF64, 2*fsize, 2*fsize)
            end
        end


        return Crystal_Electron(system, SpinPol, Spindeg, kpoints, 
                                spinsize, 
                                E_Temp, Core_Charge, TotalZ, Nocc, 0.0, 0.0, 
                                FF, Enk, Cnk)
    end
end
