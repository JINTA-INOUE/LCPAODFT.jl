struct MLWF_model
    spinsize::Int32
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    Nwann::Int32
    kmesh::Tuple{Int32,Int32,Int32}
    SpinPol::String
    SO_switch::Bool
    Dis_Energy::Vector{Float64}
    NCell::Int32
    cell_list_ijk::Vector{Vector{Int32}}
    Rdegens::Vector{Int32}
    HmnR::Array{ComplexF64,4}
    ChemP::Float64
    scf_inputfile::Vector{String}
    mlwf_inputfile::Vector{String}
end


function Print_MLWF_model(filepath::String, material::MLWF_model)

    Nwann = material.Nwann
    spinsize = material.spinsize
    SO_switch = material.SO_switch
    ChemP = material.ChemP
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    
    println("<Print_MLWF_model>")
    println("MLWF_model.jld2 read path\t$filepath")
    # println("\tNatom : $(Natom)")
    # println("\tNspecies : $(Nspecies)")
    println("\tspinsize : $(spinsize)")
    println("\tNwann : $(Nwann)")
    println("\tSO_switch : $(SO_switch)")
    println("\tChemP (Hartree) : $(ChemP)")
    println("Latvecs (AU)")
    @printf("\tA : %15.12f  %15.12f  %15.12f\n", Latvecs[1,1], Latvecs[1,2], Latvecs[1,3])
    @printf("\tB : %15.12f  %15.12f  %15.12f\n", Latvecs[2,1], Latvecs[2,2], Latvecs[2,3])
    @printf("\tC : %15.12f  %15.12f  %15.12f\n", Latvecs[3,1], Latvecs[3,2], Latvecs[3,3])

    println("Recvecs (1/AU)")
    @printf("\tA : %15.12f  %15.12f  %15.12f\n", Recvecs[1,1], Recvecs[1,2], Recvecs[1,3])
    @printf("\tB : %15.12f  %15.12f  %15.12f\n", Recvecs[2,1], Recvecs[2,2], Recvecs[2,3])
    @printf("\tC : %15.12f  %15.12f  %15.12f\n", Recvecs[3,1], Recvecs[3,2], Recvecs[3,3])

    println("")
end