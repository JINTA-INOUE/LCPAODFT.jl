struct CWF_model
    # Natom::Int32
    # Nspecies::Int32
    spinsize::Int32
    # Guide_Symbol::Vector{String}
    # Guide_Orbs_scale::Vector{String}
    # Guide_atom2spe::Vector{Int32}
    # Guide_Cutoff::Vector{Float64}
    # Guide_Total_NumOrbs::Vector{Int32}
    # Guide_MP::Vector{Int32}
    # Guide_Gxyz::Vector{Vector{Float64}}
    # Guide_Gxyz_frac::Vector{Vector{Float64}}
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    gsize::Int32
    Ngsize::Int32
    kmesh::Tuple{Int32,Int32,Int32}
    SpinPol::String
    SO_switch::Bool
    Dis_Energy::Vector{Float64}
    DMfunc::Float64
    NCell::Int32
    cell_list::Vector{UnitRange{Int32}}
    cell_list_ijk::Vector{Vector{Int32}}
    HmnR::Array{ComplexF64,4}
    ChemP::Float64
    weight_type::String
    scf_inputfile::Vector{String}
    cwf_inputfile::Vector{String}
end


function Print_CWF_model(filepath::String, material::CWF_model)

    Nwann = material.Ngsize
    spinsize = material.spinsize
    SO_switch = material.SO_switch
    ChemP = material.ChemP
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    weight_type = material.weight_type
    
    println("<Print_CWF_model>")
    println("CWF_model.jld2 read path")
    println("\t$filepath")
    # println("\tNatom : $(Natom)")
    # println("\tNspecies : $(Nspecies)")
    println("\tspinsize : $(spinsize)")
    println("\tNwann : $(Nwann)")
    println("\tSO_switch : $(SO_switch)")
    println("\tweight_type : $(weight_type)")
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