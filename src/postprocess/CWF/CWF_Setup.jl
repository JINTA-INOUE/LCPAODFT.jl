struct CWF_Setup_temp2
    filepath::String
    material::LCPAO_model
    kpoints::KPoints
    spinsize::Int32
    Guide_Symbol::Vector{String}
    Guide_OrbsString::Vector{String}
    GNatom::Int32
    Guide_Atom::Vector{Int32}
    pao2guide::Vector{Vector{Int32}}
    Guide_OrbsNumber::Vector{Int32}
    # Guide_Atoms_Cut1::Vector{Float64}
    gsize::Int32
    Ngsize::Int32
    Guide_atom2spe::Vector{Int32}
    Guide_Gxyz::Vector{Vector{Float64}}
    Guide_Gxyz_frac::Vector{Vector{Float64}}
    # CWF_FNAN::Vector{Int32}
    # CWF_natn::Vector{Vector{Int32}}
    # CWF_ncn::Vector{Vector{Int32}}
    ε::Vector{Float64}
    kBT::Vector{Float64}
    # NCell::Int32
    # cell_list::Vector{UnitRange{Int32}}
    # cell_list_ijk::Vector{Vector{Int32}}
    Ecut::Union{Float64,Nothing}
    CWF_Plot_atom::Union{Vector{Int32},Nothing}
    CWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    gtype::String
    CWF_HmnR::Bool
    CWF_Wannier::Bool
    filename::String
    verbose::Bool
end



struct CWF_Setup
    filepath::String
    material::LCPAO_model
    kpoints::KPoints
    spinsize::Int32
    Guide_Symbol::Vector{String}
    GNatom::Int32
    GNspecies::Int32
    Guide_atom2spe::Vector{Int32}
    Guide_Cutoff::Vector{Float64}
    Guide_Total_NumOrbs::Vector{Int32}
    Guide_MP::Vector{Int32}
    gsize::Int32
    Ngsize::Int32
    Guide_Gxyz::Vector{Vector{Float64}}
    Guide_Gxyz_frac::Vector{Vector{Float64}}
    CWF_FNAN::Vector{Int32}
    CWF_natn::Vector{Vector{Int32}}
    CWF_ncn::Vector{Vector{Int32}}
    Atom_Cut1::Vector{Float64}
    Grid_Origin::Vector{Float64}
    ε::Vector{Float64}
    kBT::Vector{Float64}
    Ecut::Union{Float64,Nothing}
    CWF_Plot_atom::Union{Vector{Int32},Nothing}
    CWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    Plot_Natom::Int32
    Plot_NCell::Int32
    Plot_cell_ijk::Vector{Vector{Int32}}
    gtype::String
    CWF_HmnR::Bool
    CWF_Wannier::Bool
    filename::String
    verbose::Bool
end



struct CWF_Setup_MO
    filepath::String
    material::LCPAO_model
    kpoints::KPoints
    spinsize::Int32
    Num_CWF_Grouped_Atoms::Int32
    CWF_Grouped_Atoms_EachNum::Vector{Int32}
    CWF_Grouped_Atoms::Vector{Vector{Int32}}
    Num_CWF_MOs_Group::Vector{Int32}
    CWF_MO_Selection::Vector{Vector{Int32}}
    CWF_Total_NumOrbs::Vector{Int32}
    MP3::Vector{Vector{Int32}}
    CWF_Guiding_MOs::Vector{Vector{Vector{Float64}}}
    gsize::Int32
    Ngsize::Int32
    Grid_Origin::Vector{Float64}
    ε::Vector{Float64}
    kBT::Vector{Float64}
    NCell::Int32
    cell_list::Vector{UnitRange{Int32}}
    cell_list_ijk::Vector{Vector{Int32}}
    Ecut::Union{Float64,Nothing}
    CWF_Plot_atom::Union{Vector{Int32},Nothing}
    CWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    Plot_Natom::Int32
    Plot_NCell::Int32
    Plot_cell_ijk::Vector{Vector{Int32}}
    CWF_HmnR::Bool
    CWF_Wannier::Bool
    filename::String
    verbose::Bool
end


function Print_CWF_Setup(cwf_setup::CWF_Setup)
    
    material = cwf_setup.material

    # Print material information
    filepath = cwf_setup.filepath
    Natom = material.Natom
    Nspecies = material.Nspecies
    Nspin = material.SpinP_switch+1
    atom2spe = material.atom2spe
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    Atoms_symbol = material.Atoms_symbol
    Grid_Origin = material.Grid_Origin
    GridVol = material.GridVol



    println("<Print_material>")
    println("OpenMX scfout read path")
    println("\t$filepath")
    println("")
    println("\tNatom : $(Natom)")
    println("\tNspecies : $(Nspecies)")
    println("\tNspin : $(Nspin)")
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
    println("")
    println("Real space Grid number a, b, c : $(Ngrid[1]) $(Ngrid[2]) $(Ngrid[3])")
	println("Grid_Origin:   $(Grid_Origin[1]) $(Grid_Origin[2]) $(Grid_Origin[3])")
	println("GridVol    :   $(GridVol)")
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
    println("Pseudo Atomic Orbitals : $(pao_type) type\tscale : $(scale)")
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end




    # Print CWF_Setup
    println("Guide Pseudo Atomic Orbitals : $(pao_type) type\tscale : $(scale)")
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end

    println("Guide functions Total number")

    println("CWF energy region[eV]: ")
    println("CWF kpoint sampling mesh: ")
    println("CWF real space energy cutoff[Ry]: ")
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
    Guide_Orbs::Vector{String},
    Guide_Symbol::Vector{String},
    Guide_Pos::Atompos,
    ε::Vector{Float64},
    kBT::Vector{Float64};
    kmesh::Tuple{Signed,Signed,Signed}=(1,1,1),
    Ecut = 150.0,
    CWF_Plot_atom=nothing,
    CWF_Plot_SuperCells=nothing,
    gtype::String="AO",
    filename=splitext(basename(filepath))[1],
    CWF_HmnR = false,
    CWF_Wannier = false,
    verbose = true)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()


    if kBT[1] < 0 || kBT[2] < 0
        error("please check kBT")
    end

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


    gtype = uppercase(gtype)
    if gtype ∉ ("AO", "HO")
        error("not suppot $(gtype) guide type")
    end


    if CWF_Wannier
        if isnothing(Ecut)
            error("please set Ecut")
        end

        if isnothing(CWF_Plot_atom)
            error("please set CWF_Plot_atom")
        end

        if length(Guide_Orbs) < maximum(CWF_Plot_atom)
            error("please check CWF_Plot_atom")
        end

        if isnothing(CWF_Plot_SuperCells)
            error("please set CWF_Plot_SuperCells")
        end
    end


    if length(Guide_Symbol) ≠ length(Guide_Orbs)
        error("please check Guide input")
    end


    GNatom = Guide_Pos.Natom
    if !isnothing(CWF_Plot_atom)
        if maximum(CWF_Plot_atom) > GNatom
            error("please check CWF_Plot_atom")
        end
    end


    if GNatom ≠ length(Guide_Orbs)
        error("please check input")
    end
    



    filename, ext = splitext(basename(filepath))
    if ext == ".jld2"
        material = Load_JLD2(filepath)
    else
        error("please check input file")
    end
    filepath = pwd()*"/"*filepath


    KP_flag = "Gcenter"
    kpoints = KPoints(kmesh, false; KP_flag)


    Natom = material.Natom
    Nspecies = material.Nspecies
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    Atoms_Cut1 = material.Atoms_Cut1


    Atoms_pao = material.Atoms_pao
    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    @show Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra
    for spe = 1:Nspecies
        if Spe_Symbol[spe] ∈ ("Fe", "Co", "Ni", "Cu", "Zn") 
            if Spe_extra[spe] == ""
                println("Atoms orbital = $(Spe_orb[spe])")
                error("please check input Atoms orbital")
            end
        else
            if Spe_extra[spe] ≠ ""
                println("Atoms orbital = $(Spe_orb[spe])")
                error("please check input Atoms orbital")
            end
        end
    end


    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    @show atompos_unit = Guide_Pos.unit

    Grid_Origin = zeros(Float64, 3)
    if lowercase(atompos_unit) ≠ "frac"
        xc = yc = zc = 0.0
        for atom = 1:Natom
            xc += Gxyz[atom][1]
            yc += Gxyz[atom][2]
            zc += Gxyz[atom][3]
        end
        xc = xc/Natom
        yc = yc/Natom
        zc = zc/Natom
        sn1 = 0.5*mod(Ngrid1+1, 2)
        sn2 = 0.5*mod(Ngrid2+1, 2)
        sn3 = 0.5*mod(Ngrid3+1, 2)

        xm = (div(Ngrid1,2)-sn1)*gLatvecs[1,1] + (div(Ngrid2,2)-sn2)*gLatvecs[2,1] + (div(Ngrid3,2)-sn3)*gLatvecs[3,1]
        ym = (div(Ngrid1,2)-sn1)*gLatvecs[1,2] + (div(Ngrid2,2)-sn2)*gLatvecs[2,2] + (div(Ngrid3,2)-sn3)*gLatvecs[3,2]
        zm = (div(Ngrid1,2)-sn1)*gLatvecs[1,3] + (div(Ngrid2,2)-sn2)*gLatvecs[2,3] + (div(Ngrid3,2)-sn3)*gLatvecs[3,3]
        @. Grid_Origin = [xc-xm, yc-ym, zc-zm]
    end


    
    if atompos_unit == "frac"
        Guide_Gxyz_frac = Guide_Pos.Gxyz_frac
        Guide_Gxyz_AU = Vector{Vector{Float64}}(undef, GNatom)
        for atom = 1:GNatom
            Guide_Gxyz_AU[atom] = zeros(Float64, 3)
            x = Guide_Gxyz_frac[atom][1]*Latvecs[1,1] + Guide_Gxyz_frac[atom][2]*Latvecs[2,1] + Guide_Gxyz_frac[atom][3]*Latvecs[3,1]
            y = Guide_Gxyz_frac[atom][1]*Latvecs[1,2] + Guide_Gxyz_frac[atom][2]*Latvecs[2,2] + Guide_Gxyz_frac[atom][3]*Latvecs[3,2]
            z = Guide_Gxyz_frac[atom][1]*Latvecs[1,3] + Guide_Gxyz_frac[atom][2]*Latvecs[2,3] + Guide_Gxyz_frac[atom][3]*Latvecs[3,3]
            Guide_Gxyz_AU[atom][1] = x
            Guide_Gxyz_AU[atom][2] = y
            Guide_Gxyz_AU[atom][3] = z
        end
    else
        
        Guide_Gxyz_AU = Guide_Pos.Gxyz
        Guide_Gxyz_frac = Vector{Vector{Float64}}(undef, GNatom)
        for atom = 1:GNatom
            Guide_Gxyz_frac[atom] = zeros(Float64, 3)
            Guide_Gxyz_frac[atom][1] = dot(Guide_Gxyz_AU[atom], Recvecs[1,:])*0.5/pi
            Guide_Gxyz_frac[atom][2] = dot(Guide_Gxyz_AU[atom], Recvecs[2,:])*0.5/pi
            Guide_Gxyz_frac[atom][3] = dot(Guide_Gxyz_AU[atom], Recvecs[3,:])*0.5/pi
            
            for i = 1:3
                tmp = floor(Int64, Guide_Gxyz_frac[atom][i])
                if Guide_Gxyz_frac[atom][i] > 1.0
                    Guide_Gxyz_frac[atom][i] = abs(Guide_Gxyz_frac[atom][i]-tmp)
                elseif Guide_Gxyz_frac[atom][i] < -1e-13
                    Guide_Gxyz_frac[atom][i] = abs(Guide_Gxyz_frac[atom][i]+abs(tmp)+1)
                end
            end
        end
    end

    Guide_Gxyz_AU = Gxyz




    
    GNspecies = 1
    Guide_atom2spe = zeros(Int32, 1)
    #=
    Spe_Total_NumOrbs = zeros(Int32, GNspecies)
    for spe = 1:GNspecies
        @show join(Spe_orb[spe])
        pao2string = Calc_pao2string(join(Spe_orb[spe]))
        Spe_Total_NumOrbs[spe] = length(pao2string)
    end=#


    Guide_Cutoff = zeros(Float64, GNatom)
    Guide_Total_NumOrbs = zeros(Int32, GNatom)
    #=
    for atom = 1:GNatom
        spe = Guide_atom2spe[atom]
        Guide_Cutoff[atom] = Spe_cutoff[spe]
        Guide_Total_NumOrbs[atom] = Spe_Total_NumOrbs[spe]
    end=#


    Sum = 0
    Guide_MP = zeros(Int32, GNatom+1)
    for atom = 1:Natom
        Guide_MP[atom+1] = Sum + Guide_Total_NumOrbs[atom]
        Sum += Guide_Total_NumOrbs[atom]
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


    
    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        error("please check SpinPol")
    end




    # CWF_FNAN, CWF_natn, CWF_ncn = Find_NN_Projectors(Atoms_Cut1, GNatom, Guide_Gxyz_AU, Guide_Cutoff, PAO_scale, material)
    CWF_FNAN = [0]
    CWF_natn = [[0]]
    CWF_ncn = [[0]]




    if CWF_Wannier
        Plot_Natom = length(CWF_Plot_atom)
        Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
        Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
        for cell = 1:Plot_NCell
            Plot_cell_ijk[cell] = zeros(Int32, 3)
        end

        cell = 0
        for l1 = -CWF_Plot_SuperCells[1]:CWF_Plot_SuperCells[1], l2 = -CWF_Plot_SuperCells[2]:CWF_Plot_SuperCells[2], l3 = -CWF_Plot_SuperCells[3]:CWF_Plot_SuperCells[3]
            cell += 1
            Plot_cell_ijk[cell] = [l1, l2, l3]
        end
    else
        Plot_Natom = 0
        Plot_NCell = 0
        Plot_cell_ijk = [[0, 0, 0]]
    end


    
    if verbose
        println("\n")
        @show filepath
        @show SpinPol
        @show spinsize
        # @show atom2spe
        # @show Grid_Origin
        # @show FNAN
        @show GNatom
        # @show GNspecies
        @show Guide_Symbol
        # @show Guide_atom2spe
        @show Guide_Cutoff
        @show Guide_Total_NumOrbs
        @show Guide_MP
        @show gsize
        @show Ngsize
        @show Guide_Gxyz_AU
        @show Guide_Gxyz_frac
        @show Guide_Symbol
        @show CWF_FNAN
        @show natn == CWF_natn
        @show ncn == CWF_ncn
        @show ε
        @show kBT
        @show kmesh
        # @show NCell
        @show Ecut
        @show CWF_Plot_atom
        @show CWF_Plot_SuperCells
        @show Plot_Natom
        @show Plot_NCell
        @show gtype
        @show CWF_HmnR
        @show CWF_Wannier
        @show filename
    end



    return CWF_Setup(
        filepath, material, kpoints, spinsize, 
        Guide_Symbol, GNatom, GNspecies, Guide_atom2spe, 
        Guide_Cutoff, Guide_Total_NumOrbs, Guide_MP, gsize, Ngsize,
        Guide_Gxyz_AU, Guide_Gxyz_frac,
        CWF_FNAN, CWF_natn, CWF_ncn, Atoms_Cut1, Grid_Origin,
        ε, kBT,
        Ecut, CWF_Plot_atom, CWF_Plot_SuperCells, Plot_Natom, Plot_NCell, Plot_cell_ijk,
        gtype,
        CWF_HmnR, CWF_Wannier, filename, verbose)
end


# MO type
function CWF_Setup(
    filepath::String,
    CWF_MO_Selection::Vector{Vector{Int64}},
    CWF_Grouped_Atoms::Vector{Vector{Int64}},
    ε::Vector{Float64},
    kBT::Vector{Float64};
    kmesh::Union{Nothing,Tuple{Int64,Int64,Int64}}=(1,1,1),
    Ecut = 150.0,
    CWF_Plot_atom=nothing,
    CWF_Plot_SuperCells=nothing,
    filename=splitext(basename(filepath))[1],
    CWF_HmnR = false,
    CWF_Wannier = false,
    verbose = true
)


    if kBT[1] < 0 || kBT[2] < 0
        error("please check kBT")
    end

    if !isnothing(kmesh)
        if kmesh[1] < 0 || kmesh[2] < 0 || kmesh[3] < 0
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


    if length(CWF_MO_Selection) ≠ length(CWF_Grouped_Atoms)
        error("please check input")
    end


    if CWF_Wannier
        if isnothing(Ecut)
            error("please set Ecut")
        end

        if isnothing(CWF_Plot_atom)
            error("please set CWF_Plot_atom")
        end

        if isnothing(CWF_Plot_SuperCells)
            error("please set CWF_Plot_SuperCells")
        end
    end


    



    filename, ext = splitext(basename(filepath))
    if ext == ".jld2"
        material = Load_JLD2(filepath)
    else
        error("please check input file")
    end
    filepath = pwd()*"/"*filepath

    Natom = material.Natom
    Nspecies = material.Nspecies
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn



    Ngrid = Calc_Ngrid( Ecut, Latvecs )
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3



    Grid_Origin = zeros(Float64, 3)
    xc = yc = zc = 0.0
    for atom = 1:Natom
        xc += Gxyz[atom][1]
        yc += Gxyz[atom][2]
        zc += Gxyz[atom][3]
    end
    xc = xc/Natom
    yc = yc/Natom
    zc = zc/Natom
    sn1 = 0.5*mod(Ngrid1+1, 2)
    sn2 = 0.5*mod(Ngrid2+1, 2)
    sn3 = 0.5*mod(Ngrid3+1, 2)

    xm = (div(Ngrid1,2)-sn1)*gLatvecs[1,1] + (div(Ngrid2,2)-sn2)*gLatvecs[2,1] + (div(Ngrid3,2)-sn3)*gLatvecs[3,1]
    ym = (div(Ngrid1,2)-sn1)*gLatvecs[1,2] + (div(Ngrid2,2)-sn2)*gLatvecs[2,2] + (div(Ngrid3,2)-sn3)*gLatvecs[3,2]
    zm = (div(Ngrid1,2)-sn1)*gLatvecs[1,3] + (div(Ngrid2,2)-sn2)*gLatvecs[2,3] + (div(Ngrid3,2)-sn3)*gLatvecs[3,3]
    @. Grid_Origin = [xc-xm, yc-ym, zc-zm]



    Num_CWF_Grouped_Atoms = length(CWF_MO_Selection)
    CWF_Grouped_Atoms_EachNum = zeros(Int32, Num_CWF_Grouped_Atoms)
    Num_CWF_MOs_Group = zeros(Int32, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        CWF_Grouped_Atoms_EachNum[gidx] = length(CWF_Grouped_Atoms[gidx])
        Num_CWF_MOs_Group[gidx] = length(CWF_MO_Selection[gidx])
    end

    
    
    MP3 = Vector{Vector{Int32}}(undef, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        MP3[gidx] = zeros(Int32, CWF_Grouped_Atoms_EachNum[gidx])
    end


    CWF_Total_NumOrbs = zeros(Int32, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        dim = 0
        for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            MP3[gidx][Lidx] = dim
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            dim += Total_NumOrbs[atom]
        end
        CWF_Total_NumOrbs[gidx] = dim
    end



    CWF_Guiding_MOs = Set_CWF_Guiding_MOs(material, Num_CWF_Grouped_Atoms, CWF_Grouped_Atoms_EachNum, Num_CWF_MOs_Group, CWF_Grouped_Atoms, CWF_MO_Selection)



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



    
    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        error("please check SpinPol")
    end


    
    kpoints = KPoints(kmesh, false)
    NCell, cell_list, cell_list_ijk = Get_cell_list(kmesh)



    if CWF_Wannier
        Plot_Natom = length(CWF_Plot_atom)
        Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
        Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
        for cell = 1:Plot_NCell
            Plot_cell_ijk[cell] = zeros(Int32, 3)
        end

        cell = 0
        for l1 = -CWF_Plot_SuperCells[1]:CWF_Plot_SuperCells[1], l2 = -CWF_Plot_SuperCells[2]:CWF_Plot_SuperCells[2], l3 = -CWF_Plot_SuperCells[3]:CWF_Plot_SuperCells[3]
            cell += 1
            Plot_cell_ijk[cell] = [l1, l2, l3]
        end
    else
        Plot_Natom = 0
        Plot_NCell = 0
        Plot_cell_ijk = [[0, 0, 0]]
    end



    
    if verbose
        println("\n")
        @show filepath
        @show SpinPol
        @show spinsize
        @show Grid_Origin
        @show FNAN
        @show Num_CWF_Grouped_Atoms
        @show CWF_Grouped_Atoms_EachNum
        @show CWF_Grouped_Atoms
        @show Num_CWF_MOs_Group
        @show CWF_MO_Selection
        @show MP3
        @show CWF_Total_NumOrbs
        @show gsize
        @show Ngsize
        @show ε
        @show kBT
        @show kmesh
        @show NCell
        @show Ecut
        @show CWF_Plot_atom
        @show CWF_Plot_SuperCells
        @show Plot_Natom
        @show Plot_NCell
        @show CWF_HmnR
        @show CWF_Wannier
        @show filename
    end


    
    return CWF_Setup_MO(
        filepath, material, kpoints, spinsize, 
        Num_CWF_Grouped_Atoms, CWF_Grouped_Atoms_EachNum, CWF_Grouped_Atoms,
        Num_CWF_MOs_Group, CWF_MO_Selection, CWF_Total_NumOrbs, MP3, CWF_Guiding_MOs,
        gsize, Ngsize, Grid_Origin,
        ε, kBT, NCell, cell_list, cell_list_ijk,
        Ecut, CWF_Plot_atom, CWF_Plot_SuperCells, Plot_Natom, Plot_NCell, Plot_cell_ijk,
        CWF_HmnR, CWF_Wannier, filename, verbose)
end

