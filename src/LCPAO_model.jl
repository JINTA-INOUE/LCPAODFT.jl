struct LCPAO_model
    Natom::Int32
    Nspecies::Int32
    Nspin::Int32
    atom2spe::Vector{Int32}
    Atoms_symbol::Vector{String}
    Atoms_Cut1::Vector{Float64}
    Atoms_pao::Vector{String}
    Atoms_Core_Charge::Vector{Float64}
    Init_Atoms_Nspin::Vector{Vector{Float64}}
    Init_Atoms_Angle::Vector{Vector{Float64}}
    Atoms_Angle::Vector{Vector{Float64}}
    Total_SpinS::Float64
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    Gxyz::Vector{Vector{Float64}}
    TCpyCell::Int32
    atv::Vector{Vector{Float64}}
    atv_ijk::Vector{Vector{Int32}}
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    ncn::Vector{Vector{Int32}}
    Total_NumOrbs::Vector{Int32}
    MP::Vector{Int32}
    Grid_Origin::Vector{Float64}
    Ngrid::Tuple{Int32,Int32,Int32}
    SO_switch::Bool
    SpinPol::String
    xc_type::String
    time_rev::Bool
    E_Temp::Float64
    kmesh::Tuple{Int32,Int32,Int32}
    SCF_criterion::Float64
    pao_file::Vector{String}
    pspot_file::Vector{String}
    OLP::Vector{Vector{Vector{Vector{Float64}}}}
    Hks::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
    iHks::Union{Vector{Vector{Vector{Vector{Vector{Float64}}}}}, Nothing}
    DM::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
    iDM::Union{Vector{Vector{Vector{Vector{Vector{Float64}}}}}, Nothing}
    ChemP::Float64
    Eele::Float64
    Etot::Float64
    ForceAll::Matrix{Float64}
    scf_inputfile::Vector{String}
end


function Print_LCPAO_model(filepath::String, material::LCPAO_model)

    system = "Crystal"
    Natom = material.Natom
    Nspecies = material.Nspecies
    Nspin = material.Nspin
    SpinPol = material.SpinPol
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    E_Temp = material.E_Temp
    kmesh = material.kmesh
    time_rev = material.time_rev
    Total_SpinS = material.Total_SpinS
    ChemP = material.ChemP
    SCF_criterion = material.SCF_criterion
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Atoms_symbol = material.Atoms_symbol
    Gxyz = material.Gxyz
    Ngrid = material.Ngrid
    Init_Atoms_Nspin = material.Init_Atoms_Nspin
    Init_Atoms_Angle = material.Init_Atoms_Angle
    Atoms_pao = material.Atoms_pao

    println("<Print_LCPAO_model>")
    println("LCPAO_model.jld2 read path")
    println("\t$filepath")
    println("\tSystem : $(system)")
    println("\tNatom : $(Natom)")
    println("\tNspecies : $(Nspecies)")
    println("\tNspin : $(Nspin)")
    println("\tSpinPolarization : $(SpinPol)")
    println("\tSO_switch : $(SO_switch)")
    println("\tExchange Correlation type : $(xc_type)")
    println("\tElectron Temperatue : $(E_Temp)")
    println("\tBrillouin zone sampling : $(kmesh)")
    println("\ttime reversal symmetry : $(time_rev)")
    println("\tTotal_SpinS : $(Total_SpinS)")
    println("\tChemP (Hartree) : $(ChemP)")
    println("\tSCF_criterion : $(SCF_criterion)")
    println("")
    println("Latvecs (AU)")
    @printf("\tA : %15.12f  %15.12f  %15.12f\n", Latvecs[1,1], Latvecs[1,2], Latvecs[1,3])
    @printf("\tB : %15.12f  %15.12f  %15.12f\n", Latvecs[2,1], Latvecs[2,2], Latvecs[2,3])
    @printf("\tC : %15.12f  %15.12f  %15.12f\n", Latvecs[3,1], Latvecs[3,2], Latvecs[3,3])

    println("Recvecs (1/AU)")
    @printf("\tA : %15.12f  %15.12f  %15.12f\n", Recvecs[1,1], Recvecs[1,2], Recvecs[1,3])
    @printf("\tB : %15.12f  %15.12f  %15.12f\n", Recvecs[2,1], Recvecs[2,2], Recvecs[2,3])
    @printf("\tC : %15.12f  %15.12f  %15.12f\n", Recvecs[3,1], Recvecs[3,2], Recvecs[3,3])

    println("Atom Catesian positions (AU)")
    println("\tatom\tAtom Name\t   x\t     y\t       z")
    for atom = 1:Natom
        @printf("\t%d\t%s\t\t%5.6f  %5.6f  %5.6f\n", atom, Atoms_symbol[atom], Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3])
    end
    println("")
    println("Real space Grid number a, b, c : $(Ngrid[1]) $(Ngrid[2]) $(Ngrid[3])")
    println("")

    println("Initial Number of up-spin/dn-spin per atoms")
    for atom = 1:Natom
        println("\t$atom  $(Atoms_symbol[atom])\tNup: $(Init_Atoms_Nspin[atom][1])\tNdown: $(Init_Atoms_Nspin[atom][2])")
    end

    println("Initial angler of spin/orbitals per atoms")
    for atom = 1:Natom
        println("\t$atom  $(Atoms_symbol[atom])\tspin: $(Init_Atoms_Angle[atom][1])\tspin: $(Init_Atoms_Angle[atom][2])\torbital: $(Init_Atoms_Angle[atom][3])\torbital: $(Init_Atoms_Angle[atom][4])")
    end


    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end
    println("")
end