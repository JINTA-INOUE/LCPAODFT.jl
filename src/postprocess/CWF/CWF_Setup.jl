struct CWF_Setup
    filepath::String
    material::LCPAO_model
    spinsize::Int32
    GNatom::Int32
    Wannier_Guide::Vector{Vector{Float64}}
    Guide_index::Vector{Vector{Int32}}
    Guide_Total_NumOrbs::Vector{Int32}
    gsize::Int32
    Ngsize::Int32
    Dis_Energy::Vector{Float64}
    kmesh::Tuple{Int32,Int32,Int32}
    Ecut::Union{Float64,Nothing}
    CWF_Plot_Cube::Union{Vector{Int32},Nothing}
    CWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    weight_type::String
    CWF_HmnR::Bool
    CWF_Wannier::Bool
    CWF_SOC::Bool
    CWF2MLWF::Bool
    write_coef::Bool
    filename::String
    verbose::Bool
end


struct CWF_Setup_MO
    filepath::String
    material::LCPAO_model
    spinsize::Int32
    Num_CWF_Grouped_Atoms::Int32
    CWF_Grouped_Atoms_EachNum::Vector{Int32}
    CWF_Grouped_Atoms::Vector{Vector{Int32}}
    Num_CWF_MOs_Group::Vector{Int32}
    CWF_Total_NumOrbs::Vector{Int32}
    MP3::Vector{Vector{Int32}}
    gsize::Int32
    Ngsize::Int32
    Dis_Energy::Vector{Float64}
    kmesh::Tuple{Int32,Int32,Int32}
    Ecut::Union{Float64,Nothing}
    CWF_Plot_Cube::Union{Vector{Int32},Nothing}
    CWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    weight_type::String
    CWF_HmnR::Bool
    CWF_Wannier::Bool
    CWF_SOC::Bool
    CWF2MLWF::Bool
    write_coef::Bool
    guide_out::Bool
    filename::String
    verbose::Bool
end


function Print_CWF_Setup(cwf_setup::CWF_Setup)
    
    material = cwf_setup.material

    # Print material information
    filepath = cwf_setup.filepath
    Natom = material.Natom
    Nspecies = material.Nspecies
    Nspin = material.Nspin
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    ChemP = material.ChemP
    atom2spe = material.atom2spe
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    Atoms_symbol = material.Atoms_symbol
    Atoms_pao = material.Atoms_pao



    println("<Print_material>")
    println("LCPAODFT jld2 read path")
    println("\t$filepath")
    println("")
    println("\tNatom : $(Natom)")
    println("\tNspecies : $(Nspecies)")
    println("\tNspin : $(Nspin)")
    println("\tSO_switch : $(SO_switch)")
    println("\txc_type : $(xc_type)")
    println("\tChemP : $(ChemP) Hartree")
    println("Latvecs (AU)")
    @printf("\tA : %5.10f  %5.10f  %5.10f\n", Latvecs[1,1], Latvecs[1,2], Latvecs[1,3])
    @printf("\tB : %5.10f  %5.10f  %5.10f\n", Latvecs[2,1], Latvecs[2,2], Latvecs[2,3])
    @printf("\tC : %5.10f  %5.10f  %5.10f\n", Latvecs[3,1], Latvecs[3,2], Latvecs[3,3])

    println("Recvecs (1/AU)")
    @printf("\tA : %5.10f  %5.10f  %5.10f\n", Recvecs[1,1], Recvecs[1,2], Recvecs[1,3])
    @printf("\tB : %5.10f  %5.10f  %5.10f\n", Recvecs[2,1], Recvecs[2,2], Recvecs[2,3])
    @printf("\tC : %5.10f  %5.10f  %5.10f\n", Recvecs[3,1], Recvecs[3,2], Recvecs[3,3])

    println("Atom Catesian positions (AU)")
    println("\tatom\tAtom Name\t   x\t     y\t       z")
    for atom = 1:Natom
        @printf("\t%d\t%s\t\t%5.6f  %5.6f  %5.6f\n", atom, Atoms_symbol[atom], Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3])
    end

    println("Pseudo Atomic Orbitals")
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end
    println("")




    Guide_index = cwf_setup.Guide_index
    Guide_Total_NumOrbs = cwf_setup.Guide_Total_NumOrbs
    Ngsize = cwf_setup.Ngsize
    Dis_Energy = cwf_setup.Dis_Energy
    kmesh = cwf_setup.kmesh
    Ecut = cwf_setup.Ecut
    CWF_Plot_Cube = cwf_setup.CWF_Plot_Cube
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    weight_type = cwf_setup.weight_type
    CWF_HmnR = cwf_setup.CWF_HmnR
    CWF_Wannier = cwf_setup.CWF_Wannier
    CWF_SOC = cwf_setup.CWF_SOC
    CWF2MLWF = cwf_setup.CWF2MLWF
    write_coef = cwf_setup.write_coef
    filename = cwf_setup.filename
    verbosity = cwf_setup.verbose

    # Print CWF_Setup
    println("<Print_CWF_Setup>")
    println("Guide PAO")
    println("\tatom\tAtom Name")
    for atom = 1:Natom
        @printf("\t%d\t%s\t\t", atom, Atoms_symbol[atom])
        for i = 1:Guide_Total_NumOrbs[atom]
            @printf("  %d ", Guide_index[atom][i])
        end
        @printf("\n")
    end
    println("\tCWF Number: $(Ngsize)")
    println("\tweight_type: $(weight_type)")
    println("\tCWF energy region[eV]: $(Dis_Energy*27.2113845)")
    println("\tCWF kpoint sampling mesh: $(kmesh)")
    println("\tCalc HmnR: $(CWF_HmnR)")
    println("\tWrite CWF Cube: $(CWF_Wannier)")
    if CWF_Wannier
        println("\tCWF real space energy cutoff[Ry]: $(Ecut)")
        println("\tCWF_Plot_Cube: $(CWF_Plot_Cube)")
        println("\tCWF_Plot_SuperCells: $(CWF_Plot_SuperCells)")
    end
    println("\tCWF_SOC: $(CWF_SOC)")
    println("\tCWF2MLWF: $(CWF2MLWF)")
    println("\twrite_coef: $(write_coef)")
    println("\tfilename: $(filename)")
    println("\tverbosity: $(verbosity)")
    println("")
end



function Print_CWF_Setup(cwf_setup::CWF_Setup_MO)
    
    material = cwf_setup.material

    # Print material information
    filepath = cwf_setup.filepath
    Natom = material.Natom
    Nspecies = material.Nspecies
    Nspin = material.Nspin
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    atom2spe = material.atom2spe
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    Atoms_symbol = material.Atoms_symbol
    Atoms_pao = material.Atoms_pao



    println("<Print_material>")
    println("LCPAODFT jld2 read path")
    println("\t$filepath")
    println("")
    println("\tNatom : $(Natom)")
    println("\tNspecies : $(Nspecies)")
    println("\tNspin : $(Nspin)")
    println("\tSO_switch : $(SO_switch)")
    println("\txc_type : $(xc_type)")
    println("Latvecs (AU)")
    @printf("\tA : %5.10f  %5.10f  %5.10f\n", Latvecs[1,1], Latvecs[1,2], Latvecs[1,3])
    @printf("\tB : %5.10f  %5.10f  %5.10f\n", Latvecs[2,1], Latvecs[2,2], Latvecs[2,3])
    @printf("\tC : %5.10f  %5.10f  %5.10f\n", Latvecs[3,1], Latvecs[3,2], Latvecs[3,3])

    println("Recvecs (1/AU)")
    @printf("\tA : %5.10f  %5.10f  %5.10f\n", Recvecs[1,1], Recvecs[1,2], Recvecs[1,3])
    @printf("\tB : %5.10f  %5.10f  %5.10f\n", Recvecs[2,1], Recvecs[2,2], Recvecs[2,3])
    @printf("\tC : %5.10f  %5.10f  %5.10f\n", Recvecs[3,1], Recvecs[3,2], Recvecs[3,3])

    println("Atom Catesian positions (AU)")
    println("\tatom\tAtom Name\t   x\t     y\t       z")
    for atom = 1:Natom
        @printf("\t%d\t%s\t\t%5.6f  %5.6f  %5.6f\n", atom, Atoms_symbol[atom], Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3])
    end

    println("Pseudo Atomic Orbitals")
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end
    println("")




    Ngsize = cwf_setup.Ngsize
    Dis_Energy = cwf_setup.Dis_Energy
    kmesh = cwf_setup.kmesh
    Ecut = cwf_setup.Ecut
    CWF_Plot_Cube = cwf_setup.CWF_Plot_Cube
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    weight_type = cwf_setup.weight_type
    CWF_HmnR = cwf_setup.CWF_HmnR
    CWF_Wannier = cwf_setup.CWF_Wannier
    CWF_SOC = cwf_setup.CWF_SOC
    CWF2MLWF = cwf_setup.CWF2MLWF
    write_coef = cwf_setup.write_coef
    filename = cwf_setup.filename
    verbosity = cwf_setup.verbose

    # Print CWF_Setup
    println("<Print_CWF_Setup>")
    println("\tCWF Number: $(Ngsize)")
    println("\tweight_type: $(weight_type)")
    println("\tCWF energy region[eV]: $(Dis_Energy*27.2113845)")
    println("\tCWF kpoint sampling mesh: $(kmesh)")
    println("\tCalc HmnR: $(CWF_HmnR)")
    println("\tWrite CWF Cube: $(CWF_Wannier)")
    if CWF_Wannier
        println("\tCWF real space energy cutoff[Ry]: $(Ecut)")
        println("\tCWF_Plot_Cube: $(CWF_Plot_Cube)")
        println("\tCWF_Plot_SuperCells: $(CWF_Plot_SuperCells)")
    end
    println("\tCWF_SOC: $(CWF_SOC)")
    println("\tCWF2MLWF: $(CWF2MLWF)")
    println("\twrite_coef: $(write_coef)")
    println("\tfilename: $(filename)")
    println("\tverbosity: $(verbosity)")
    println("")
end


"""
```
    CWF_Setup(...)

Setup for Generating Closest Wannier Functions

Mandatory arguments:

- `filepath`
- `Guide_Orbs`
- `Guide_Symbol`
- `Guide_Pos`
- `ε`
- `kBT`

Thw following is the most commonly used optional arguments:


function CWF_Setup(
    filepath::String,
    Guide_Orbs::Vector{String},
    Guide_Symbol::Vector{String},
    Guide_Pos::Atompos,
    ε::Vector{Float64},
    kBT::Vector{Float64};
    kmesh::Union{Nothing,Tuple{Int64,Int64,Int64}}=(1,1,1),
    Ecut = 150.0,
    CWF_Plot_atom=nothing,
    CWF_Plot_SuperCells=nothing,
    RotMat = nothing,
    PAO_scale::Vector{Float64} = ones(Float64,length(Guide_Symbol)),
    Plot_line::String = "",
    gtype::String="AO",
    filename=splitext(basename(filepath))[1],
    CWF_HmnR = false,
    CWF_Wannier = false,
    verbose = true
)


# Examples
```
Example case please check https://github.com/JINTA-INOUE/juOpenMX
```
"""
function CWF_Setup(
    filepath::String,
    Guide_index,
    Dis_Energy::Vector{Float64};
    kmesh::Tuple{Signed,Signed,Signed}=(1,1,1),
    Ecut = 150.0,
    CWF_Plot_Cube=nothing,
    CWF_Plot_SuperCells=nothing,
    weight_type::String = "Fermi",
    filename=splitext(basename(filepath))[1],
    CWF_HmnR::Bool = false,
    CWF_Wannier::Bool = false,
    CWF_SOC::Bool = false,
    CWF2MLWF::Bool = false,
    write_coef::Bool = false,
    verbose::Bool = true)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()


    # Input Error check
    weight_type = lowercase(weight_type)
    if weight_type ∉ ("fermi", "poly")
        error("please check weight_type")
    end
    
    if length(Dis_Energy) ≠ 4
        error("please check Dis_Energy")
    end

    if weight_type=="fermi"
        if Dis_Energy[3] < 0.0 || Dis_Energy[4] < 0.0
            error("please check kBT")
        end
    end
    @. Dis_Energy = Dis_Energy/27.2113845

    if !isnothing(kmesh)
        if kmesh[1] <= 0 || kmesh[2] <= 0 || kmesh[3] <= 0
            error("please check kmesh")
        end

        if any(iseven.(kmesh))
            error("must set kmesh is odd number")
        end
    end

    if !isnothing(Ecut)
        if Ecut < 0.0
            error("please check Ecut")
        end
    end

    if CWF_Wannier
        if isnothing(Ecut)
            error("please set Ecut")
        end

        if nprocs >= 2
            MPI.Finalized()
            error("please run serial.")
        end

        if isnothing(CWF_Plot_Cube)
            error("please set CWF_Plot_Cube")
        end

        if isnothing(CWF_Plot_SuperCells)
            error("please set CWF_Plot_SuperCells")
        end

        if length(CWF_Plot_SuperCells) ≠ 3 || all(x->x>0, CWF_Plot_SuperCells)
            error("please check CWF_Plot_SuperCells")
        end
    end


    if CWF2MLWF
        if nprocs >= 2
            MPI.Finalized()
            error("please run serial.")
        end

        if weight_type == "fermi"
            MPI.Finalized()
            error("please use weight_type = \"poly\"")
        end

        if lowercase(weight_type) == "fermi"
            error("please set weight_type = Poly.")
        end
    end





    material = Load_LCPAODFT_model(filepath)

    filepath = pwd()*"/"*filepath



    Natom = material.Natom
    Latvecs = material.Latvecs
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol



    GNatom = length(Guide_index)
    Guide_Total_NumOrbs = zeros(Int32, GNatom)
    for atom = 1:GNatom
        Guide_Total_NumOrbs[atom] = length(Guide_index[atom])
    end

    
    fsize = sum(Total_NumOrbs)
    gsize = sum(Guide_Total_NumOrbs)
    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
        Ngsize = gsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
        Ngsize = 2*gsize
    end

    if CWF_Wannier
        for p in CWF_Plot_Cube
            if 1 > p || p > Ngsize
                error("please check CWF_Ploy_Cube")
            end
        end
    end


    Wannier_Guide = Vector{Vector{Float64}}(undef, Ngsize)
    for p = 1:Ngsize
        Wannier_Guide[p] = zeros(Float64, 3)
    end 

    p = 0
    for atom = 1:Natom, _ = 1:length(Guide_index[atom])
        p += 1
        Wannier_Guide[p][1] = Gxyz[atom][1]
        Wannier_Guide[p][2] = Gxyz[atom][2]
        Wannier_Guide[p][3] = Gxyz[atom][3]
        if SpinPol == "nc"
            Wannier_Guide[p+gsize][1] = Gxyz[atom][1]
            Wannier_Guide[p+gsize][2] = Gxyz[atom][2]
            Wannier_Guide[p+gsize][3] = Gxyz[atom][3]
        end
    end

    
    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        error("please check SpinPol")
    end



    cwf_setup = CWF_Setup(
        filepath, material, spinsize, 
        GNatom, Wannier_Guide, Guide_index, Guide_Total_NumOrbs, gsize, Ngsize,
        Dis_Energy, kmesh,
        Ecut, CWF_Plot_Cube, CWF_Plot_SuperCells,
        weight_type,
        CWF_HmnR, CWF_Wannier, CWF_SOC, CWF2MLWF, write_coef, 
        filename, verbose)

    if myrank == 0
        Print_CWF_Setup(cwf_setup)
    end
    MPI.Barrier(comm)


    return cwf_setup
end


# MO type
function CWF_Setup(
    filepath::String,
    Num_CWF_MOs_Group::Vector{Int64},
    CWF_MO_Grouped_Atoms::Vector{Int64},
    Dis_Energy::Vector{Float64};
    kmesh::Union{Nothing,Tuple{Int64,Int64,Int64}}=(1,1,1),
    Ecut = 100.0,
    CWF_Plot_Cube=nothing,
    CWF_Plot_SuperCells=nothing,
    weight_type::String = "Fermi",
    filename = nothing,
    CWF_HmnR::Bool = false,
    CWF_Wannier::Bool = false,
    CWF_SOC::Bool = false,
    CWF2MLWF::Bool = false,
    write_coef::Bool = false,
    guide_out::Bool = false,
    verbose::Bool = true)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()


    # Input Error check
    weight_type = lowercase(weight_type)
    if weight_type ∉ ("fermi", "poly")
        myrank == 0 && @show weight_type
        myrank == 0 && error("please check weight_type")
    end
    
    if length(Dis_Energy) ≠ 4
        myrank == 0 && @show Dis_Energy
        error("please check Dis_Energy")
    end

    if weight_type=="fermi"
        if Dis_Energy[3] < 0.0 || Dis_Energy[4] < 0.0
            error("please check kBT")
        end
    end
    @. Dis_Energy = Dis_Energy/27.2113845

    if !isnothing(kmesh)
        if kmesh[1] <= 0 || kmesh[2] <= 0 || kmesh[3] <= 0
            error("please check kmesh")
        end

        if any(iseven.(kmesh))
            error("must set kmesh is odd number")
        end
    end

    if !isnothing(Ecut)
        if Ecut < 0.0
            error("please check Ecut")
        end
    end


    if CWF_Wannier

        if nprocs >= 2
            MPI.Finalized()
            error("please run serial.")
        end
        
        if isnothing(Ecut)
            error("please set Ecut")
        end

        if isnothing(CWF_Plot_Cube)
            error("please set CWF_Plot_Cube")
        end

        if isnothing(CWF_Plot_SuperCells)
            error("please set CWF_Plot_SuperCells")
        end

        if length(CWF_Plot_SuperCells) ≠ 3 || all(x->x>0, CWF_Plot_SuperCells)
            error("please check CWF_Plot_SuperCells")
        end
    end


    if CWF2MLWF
        if nprocs >= 2
            MPI.Finalized()
            error("please run serial.")
        end

        if weight_type == "fermi"
            MPI.Finalized()
            error("please use weight_type = \"poly\"")
        end
        
        if lowercase(weight_type) == "fermi"
            error("please set weight_type = Poly.")
        end
    end

    


    material = Load_LCPAODFT_model(filepath)
    filepath = pwd()*"/"*filepath

    Natom = material.Natom
    Latvecs = material.Latvecs
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol

    if SpinPol == "nc"
        error("not support yet.")
    end

    tmp_flag = zeros(Int32, Natom)
    for atom = 1:Natom
        k = CWF_MO_Grouped_Atoms[atom]
        tmp_flag[k] += 1
    end

    Num_CWF_Grouped_Atoms = 0
    for atom = 1:Natom
        if tmp_flag[atom] > 0
            Num_CWF_Grouped_Atoms += 1
        end
    end


    CWF_Grouped_Atoms_EachNum = zeros(Int32, Num_CWF_Grouped_Atoms)
    for atom = 1:Num_CWF_Grouped_Atoms
        CWF_Grouped_Atoms_EachNum[atom] = tmp_flag[atom]
    end


    CWF_Grouped_Atoms = Vector{Vector{Int32}}(undef, Num_CWF_Grouped_Atoms)
    for atom = 1:Num_CWF_Grouped_Atoms
        CWF_Grouped_Atoms[atom] = zeros(Int32, CWF_Grouped_Atoms_EachNum[atom])
    end    


    fill!(tmp_flag, 1)
    for atom = 1:Natom
        k = CWF_MO_Grouped_Atoms[atom]
        CWF_Grouped_Atoms[k][tmp_flag[k]] = atom
        tmp_flag[k] += 1
    end


    CWF_Total_NumOrbs = zeros(Int32, Num_CWF_Grouped_Atoms)
    MP3 = Vector{Vector{Int32}}(undef, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        MP3[gidx] = zeros(Int32, CWF_Grouped_Atoms_EachNum[gidx])
    end

    for gidx = 1:Num_CWF_Grouped_Atoms
        dim = 0
        for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            MP3[gidx][Lidx] = dim
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            dim += Total_NumOrbs[atom]
        end
        CWF_Total_NumOrbs[gidx] = dim
    end



    gsize = 0
    for gidx = 1:Num_CWF_Grouped_Atoms, Lidx = 1:Num_CWF_MOs_Group[gidx]
        gsize += 1
    end

    fsize = sum(Total_NumOrbs)
    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
        Ngsize = gsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
        Ngsize = 2*gsize
    end

    if CWF_Wannier
        for p in CWF_Plot_Cube
            if 1 > p || p > Ngsize
                error("plase check CWF_Ploy_Cube")
            end
        end
    end



    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        error("please check SpinPol")
    end

    guide_out = false

    cwf_setup = CWF_Setup_MO(
        filepath, material, spinsize, 
        Num_CWF_Grouped_Atoms, CWF_Grouped_Atoms_EachNum, CWF_Grouped_Atoms,
        Num_CWF_MOs_Group, CWF_Total_NumOrbs, MP3,
        gsize, Ngsize,
        Dis_Energy, kmesh, Ecut, 
        CWF_Plot_Cube, CWF_Plot_SuperCells,
        weight_type, CWF_HmnR, CWF_Wannier, CWF_SOC, CWF2MLWF, write_coef,
        guide_out, filename, verbose)

    if myrank == 0
        Print_CWF_Setup(cwf_setup)
    end


    return cwf_setup
end

