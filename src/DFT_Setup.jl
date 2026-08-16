struct DFT_Setup
    Natom::Int32
    Nspecies::Int32
    Nspin::Int32
    spinsize::Int32
    atom2spe::Vector{Int32}
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    gLatvecs::Matrix{Float64}
    gRecvecs::Matrix{Float64}
    Gxyz::Vector{Vector{Float64}}
    Gxyz_frac::Vector{Vector{Float64}}
    GridVol::Float64
    Grid_Origin::Vector{Float64}
    Atoms_symbol::Vector{String}
    Atoms_cutoff::Vector{Float64}
    Atoms_pao::Vector{String}
    system::String
    Init_Atoms_Nspin::Vector{Vector{Float64}}
    Init_Atoms_Angle::Vector{Vector{Float64}}
    SpinPol::String
    SO_switch::Bool
    xc_type::String
    pao_file::Vector{String}
    pspot_file::Vector{String}
    Ngrid::Tuple{Int32,Int32,Int32}
    Mixing_method::String
    SCF_criterion::Float64
    SCF_max::Int32
    Init_Mixing_weight::Float64
    Min_Mixing_weight::Float64
    Max_Mixing_weight::Float64
    Num_Mixing_Pulay::Int32
    Start_Pulay_SCF::Int32
    E_Temp::Float64
    kmesh::Tuple{Int32,Int32,Int32}
    symmetry::Symmetry
    Hub_U::Bool
    Hub_U_atom::Vector{Vector{Float64}}
    Hub_U_orbpol::Vector{Bool}
    Hub_U_occ::String
    Hub_Type::String
    dc_Type::String
    time_rev::Bool
    cal_force::Bool
    cal_mode::Int32
    fileout::Bool
    filename::String
    restart::Bool
    filepath::String
    send_email::Bool
    verbosity::Int64
end


function Print_DFT_Setup(dft_setup::DFT_Setup)

    Natom = dft_setup.Natom
    Nspecies = dft_setup.Nspecies
    Nspin = dft_setup.Nspin
    Latvecs = dft_setup.Latvecs
    Recvecs = dft_setup.Recvecs
    Gxyz = dft_setup.Gxyz
    Atoms_symbol = dft_setup.Atoms_symbol
    Init_Atoms_Nspin = dft_setup.Init_Atoms_Nspin
    Init_Atoms_Angle = dft_setup.Init_Atoms_Angle
    Atoms_pao = dft_setup.Atoms_pao
    Grid_Origin = dft_setup.Grid_Origin
    GridVol = dft_setup.GridVol
    Ngrid = dft_setup.Ngrid


    println("<Print_DFT_Setup>")
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
    for spe = 1:Nspecies
        println("\t$spe  $(Spe_Symbol[spe])\tcutoff: $(Spe_cutoff[spe])\torbitals $(Spe_orb[spe]*Spe_extra[spe])")
    end


    system = dft_setup.system
    SpinPol = dft_setup.SpinPol
    SO_switch = dft_setup.SO_switch
    xc_type = dft_setup.xc_type
    E_Temp = dft_setup.E_Temp
    kmesh = dft_setup.kmesh
    time_rev = dft_setup.time_rev
    cal_force = dft_setup.cal_force
    cal_mode = dft_setup.cal_mode
    println("")
    println("<SCF Setup>")
    println("\tSystem : $(system)")
    println("\tSpinPolarization : $(SpinPol)")
    println("\tSpinOrbitalCoupling : $(SO_switch)")
    println("\tExchange Correlation type : $(xc_type)")
    println("\tElectron Temperatue : $(E_Temp)")
    println("\tBrillouin zone sampling : $(kmesh)")
    println("\ttime reversal symmetry : $(time_rev)")
    println("\tcal_force : $(cal_force)")
    println("\tcal_mode : $(cal_mode)")



    Mixing_method = dft_setup.Mixing_method
    SCF_criterion = dft_setup.SCF_criterion
    SCF_max = dft_setup.SCF_max
    Init_Mixing_weight = dft_setup.Init_Mixing_weight
    Min_Mixing_weight = dft_setup.Min_Mixing_weight
    Max_Mixing_weight = dft_setup.Max_Mixing_weight
    Num_Mixing_Pulay = dft_setup.Num_Mixing_Pulay
    Start_Pulay_SCF = dft_setup.Start_Pulay_SCF
    println("")
    println("<Mixing Parameters>")
    println("\tMixing_method : $(Mixing_method)")
    println("\tSCF_criterion : $(SCF_criterion)")
    println("\tSCF_max : $(SCF_max)")
    println("\tInit_Mixing_weight : $(Init_Mixing_weight)")
    println("\tMin_Mixing_weight : $(Min_Mixing_weight)")
    println("\tMax_Mixing_weight : $(Max_Mixing_weight)")
    println("\tNum_Mixing_Pulay : $(Num_Mixing_Pulay)")
    println("\tStart_Pulay_SCF : $(Start_Pulay_SCF)")
    

    Hub_U = dft_setup.Hub_U
    Hub_U_atom = dft_setup.Hub_U_atom
    Hub_U_orbpol = dft_setup.Hub_U_orbpol
    Hub_U_occ = dft_setup.Hub_U_occ
    Hub_Type = dft_setup.Hub_Type
    dc_Type = dft_setup.dc_Type
    if Hub_U
        println("")
        println("<Hubbard U>")
        println("\tHub_U_occ : $(Hub_U_occ)")
        println("\tHubbard_Type : $(Hub_Type)")
        println("\tdc_Type : $(dc_Type)")
        println("\tHubbard Orbital Polarization")
        for atom = 1:Natom
            if Hub_U_orbpol[atom]
                println("\t\t$atom  $(Atoms_symbol[atom])\ton")
            else
                println("\t\t$atom  $(Atoms_symbol[atom])\toff")
            end
        end
        println("\tHubbard U Energy(eV)")
        for spe = 1:Nspecies
            Spe_MaxL_Basis, Spe_Num_Basis = get_ialpha_index(Spe_orb[spe])
            @printf("\t\t%d\t%s\t", spe, Spe_Symbol[spe])
            counts = 0
            for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
                counts += 1
                if l == 0
                    @printf("%ds:%5.2f  ", p, Hub_U_atom[spe][counts])
                elseif l == 1
                    @printf("%dp:%5.2f  ", p, Hub_U_atom[spe][counts])
                elseif l == 2
                    @printf("%dd:%5.2f  ", p, Hub_U_atom[spe][counts])
                elseif l == 3
                    @printf("%df:%5.2f  ", p, Hub_U_atom[spe][counts])
                end
            end
            @printf("\n")
        end
    end


    println("")
    Print_Symmetry(dft_setup.symmetry)
    
    fileout = dft_setup.fileout
    filename = dft_setup.filename
    restart = dft_setup.restart
    filepath = dft_setup.filepath
    verbosity = dft_setup.verbosity
    send_email = dft_setup.send_email
    println("")
    println("<Outputs>")
    println("\tfileout : $(fileout)")
    println("\tfilename : $(filename)")
    println("\trestart : $(restart)")
    if restart
        println("\tfilepath : $(filepath)")
    end
    println("\tverbosity : $(verbosity)")
    println("\tsend_email : $(send_email)")
end


function RestartFile_check(Natom, Nspin, Latvecs, Gxyz, Grid_Origin, SpinPol, SO_switch, xc_type, filepath)

    material = Load_LCPAODFT_model(filepath)
    Read_Natom = material.Natom
    Read_Nspin = material.Nspin
    Read_Latvecs = material.Latvecs
    Read_Gxyz = material.Gxyz
    Read_Grid_Origin = material.Grid_Origin
    Read_SpinPol = material.SpinPol
    Read_SO_switch = material.SO_switch
    Read_xc_type = material.xc_type

    function check_same(data1, data2, val_name::String)
        if data1 ≠ data2
            error("Failed Restart file ($val_name)")
        end
    end

    check_same(Natom, Read_Natom, "Natom")
    check_same(Nspin, Read_Nspin, "Nspin")
    check_same(Latvecs, Read_Latvecs, "Latvecs")
    check_same(Gxyz, Read_Gxyz, "Gxyz")
    check_same(Grid_Origin, Read_Grid_Origin, "Grid_Origin")
    check_same(SpinPol, Read_SpinPol, "SpinPol")
    check_same(SO_switch, Read_SO_switch, "SO_switch")
    check_same(xc_type, Read_xc_type, "xc_type")
end


"""
```
    DFT_Setup(...)

Setup for solving Kohn-Sham equations using Self Consistent calculation.

Mandatory arguments:

- `lattice`: an instance of `Lattice`
- `Atoms_orb`: 
- `Atoms_symbol`: 
- `Atoms_pos`:  an instance of `Atompos`
- `system`: system name (`Atom`, `Cluster`, `Crystal`)

Thw following is the most commonly used optional arguments:
- `Atoms_Nspin` : output file name
- `Atoms_Angle` : output file name
- `SpinPol` : output file name
- `SO_switch` : output file name
- `xc_type` : output file name
- `Ngrid` : output file name
- `Mixing_method` : SCF Mixing method [`RMM-DIISH`]
- `SCF_criterion`: SCF criterion
- `SCF_max` : SCF iteration max
- `Init_Mixing_weight` : Initial Mixing weight
- `Min_Mixing_weight` : minimum Mixing weight
- `Max_Mixing_weight` : minimum Mixing weight
- `Num_Mixing_Pulay` : start SCF number of DMM-DIIS mixing
- `Start_Pulay_SCF` : start SCF number of DMM-DIIS mixing
- `E_Temp`: Temperatue with eV units
- `kmesh` : kpoint smapling mesh in first BZ
- `time_rev` : which use time reversal symmetry
- `verbosity` : terminal print
- `fileout` : output bool
- `filename` : output file name
- `send_email` : send results email


function DFT_Setup(
    lattice::Lattice,
    Atoms_orb::Vector{String},
    Atoms_symbol::Vector{String},
    Atoms_pos::Atompos,
    system::AbstractString;
    Atoms_Nspin = nothing,
    Atoms_Angle = nothing,
    SpinPol::String = "off",
    SO_switch::Bool = false,
    xc_type::String = "LDA",
    Ecut::AbstractFloat = 150.0,
    Ngrid = nothing,
    Mixing_method::AbstractString = "RMM-DIISH",
    SCF_criterion::AbstractFloat = 1e-6,
    SCF_max::Signed = 10,
    Init_Mixing_weight::AbstractFloat = 0.3,
    Min_Mixing_weight::AbstractFloat = 0.001,
    Max_Mixing_weight::AbstractFloat = 0.4,
    Num_Mixing_Pulay::Signed = 5,
    Start_Pulay_SCF::Signed = 6,
    E_Temp::Union{AbstractFloat,Signed} = 300.0,
    kmesh::Tuple{Signed,Signed,Signed} = (1,1,1),
    cal_mode::Singed = 1,
    Hub_U::Bool = false,
    Hub_U_atom::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    Hub_U_occ::AbstractString = "dual",
    Hub_Type::AbstractString = "Dudarev",
    dc_Type::AbstractString = "sFLL",
    verbosity::Int = 1,
    fileout::Bool = false,
    filename::AbstractString = PROGRAM_FILE,
    filepath::AbstractString = nothing,
    send_email::Bool = false
)
"""
function DFT_Setup(
    lattice::Lattice,
    Atoms_orb::Vector{String},
    Atoms_symbol::Vector{String},
    Atoms_pos::Atompos,
    system::AbstractString;
    Atoms_Nspin = nothing,
    Atoms_Angle = nothing,
    SpinPol::String = "off",
    SO_switch::Bool = false,
    xc_type::String = "LDA",
    Ecut::AbstractFloat = 150.0,
    Ngrid = nothing,
    SCF_criterion::AbstractFloat = 1e-6,
    SCF_max::Signed = 10,
    Init_Mixing_weight::AbstractFloat = 0.3,
    Min_Mixing_weight::AbstractFloat = 0.001,
    Max_Mixing_weight::AbstractFloat = 0.4,
    Num_Mixing_Pulay::Signed = 5,
    Start_Pulay_SCF::Signed = 6,
    E_Temp::Union{AbstractFloat,Signed} = 300.0,
    kmesh::Tuple{Signed,Signed,Signed} = (1,1,1),
    cal_force::Bool = false,
    cal_mode::Signed = 1,
    Hub_U::Bool = false,
    Hub_U_atom::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    Hub_U_orbpol::Union{Nothing,Vector{Bool}} = nothing,
    Hub_U_occ::AbstractString = "dual",
    Hub_Type::AbstractString = "Dudarev",
    dc_Type::AbstractString = "sFLL",
    verbosity::Int64 = 1,
    fileout::Bool = true,
    filename::AbstractString = PROGRAM_FILE,
    filepath::Union{Nothing,AbstractString} = nothing,
    send_email::Bool = false)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()

    if myrank == 0 && verbosity >= 1
        # cpu_info = Sys.cpu_info()[1]
        # println("\t$cpu_info")
        # println("")
        println("<DFT MPI process/BLAS>")
        println("\t$nprocs MPI processes and $nthreads threads, $nblas BLAS threads")
        println("\t$(now())")
        println("")
    end
    MPI.Barrier(comm)

    if Hub_U
        error("not support yet.")
    end

    Mixing_method = "RMM-DIISH"


    Latvecs = lattice.Latvecs
    Natom = Atoms_pos.Natom
    Nspecies = length(unique(Atoms_symbol))


    if Natom ≠ length(Atoms_symbol)
        println("please check atompos and atomsymbol, Atoms_orb")
        error("please check input files")
    end


    system = lowercase(strip(system))
    if system ∈ ("atom", "atoms", "cluster")
        system = "Cluster"
    elseif system ∈ ("band", "bands", "crystal", "crystals")
        system = "Crystal"
    else
        println("only support system = Crystal")
		error("please check system")
    end


    if !isnothing(Atoms_Nspin)
        if length(Atoms_Nspin) ≠ Natom
            error("please check Atoms_Nspin")
        end
    end

    if !isnothing(Atoms_Angle)
        if length(Atoms_Angle) ≠ Natom
            error("please check Atoms_Angle")
        end
    end

    SpinPol = lowercase(strip(SpinPol))
	if SpinPol ∉ ("off", "on", "nc")
        println("Now SpinPol is $SpinPol")
        error("please check SpinPol")
    end

    if SpinPol ∈ ("off", "on") && SO_switch
        error("please check SpinPol and SO_switch")
    end

    if SpinPol ∈ ("on", "nc") && xc_type == "LDA"
        error("please check SpinPol and xc_type")
    end

    xc_type = lowercase(strip(xc_type))
	if xc_type == "lda"
		xc_type = "LDA"
	elseif xc_type == "lsda"
		xc_type = "LSDA"
	elseif xc_type ∈ ("gga-pbe", "gga_pbe", "gga", "pbe")
		xc_type = "GGA_PBE"
	else
        println("not support $xc_type")
		error("please check xc_type")
	end

    if Ecut < 0.0
        error("please check Ecut")
    end

    if !isnothing(Ngrid)
        if length(Ngrid) ≠ 3
            error("please check Ngrid")
        end

        for i = 1:3
            if Ngrid[i] <= 0
                error("please check Ngrid")
            end
        end
    end

    if SCF_criterion < 0.0
        error("please check SCF_criterion")
    end

    if SCF_max < 1
        error("please check SCF_max")
    end

    if Init_Mixing_weight <= 0.0
        error("please check Init_Mixing_weight")
    end

    if Min_Mixing_weight <= 0.0
        error("please check Min_Mixing_weight")
    end

    if Max_Mixing_weight <= 0.0
        error("please check Max_Mixing_weight")
    end

    if Num_Mixing_Pulay <= 0
        error("please check Num_Mixing_Pulay")
    end

    if Start_Pulay_SCF <= 0
        error("please check Start_Pulay_SCF")
    end

    if E_Temp <= 0.0
        error("please check Max_Mixing_weight")
    end

    if kmesh[1] <= 0 || kmesh[2] <= 0 || kmesh[3] <= 0
        error("please check kmesh")
    end

    if cal_mode ∉ (1, 2)
        throw(ArgumentError("cal_mode must be either 1 or 2 (got $cal_mode)"))
    end

    if SpinPol ∈ ("off", "on") && system == "Crystal" && nprocs > div(prod(kmesh),2)
        error("not support number of process > number of Total kmesh points")
    elseif SpinPol == "nc" && system == "Crystal" && nprocs > prod(kmesh)
        error("not support number of process > number of Total kmesh points")
    end


    if Hub_U
        if isnothing(Hub_U_atom)
            error("please input Hubbard_U_atom")
        end

        if isnothing(Hub_U_orbpol)
            error("please input Hubbard_U_atom")
        end

        if length(Hub_U_atom) ≠ length(Atoms_symbol)
            error("please check Hub_U_atom")
        end

        if length(Hub_U_orbpol) ≠ Natom
            error("please check Hub_U_orbpol")
        end

        Hub_U_occ = lowercase(Hub_U_occ)
        if Hub_U_occ ≠ "dual"
            error("please check Hub_U_occ")
        end

        Hub_Type = lowercase(Hub_Type)
        if Hub_Type ≠ "dudarev"
            error("please check Hub_Type")
        end

        dc_Type = lowercase(dc_Type)
        if dc_Type ≠ "sfll"
            error("please check dc_Type")
        end
    else
        Hub_U_atom = [[0.0]]
        Hub_U_orbpol = [false]
    end


    restart = ifelse(isnothing(filepath), false, true)
    if restart
        restart_filename, ext = splitext(basename(filepath))
        if ext ≠ ".jld2"
            error("please check filepath.")
        end
        filepath2 = filepath
    else
        filepath2 = ""
    end
    
    

    if SpinPol == "off"
        Nspin = 1
        spinsize = 1
    elseif SpinPol == "on"
        Nspin = 2
        spinsize = 2
    elseif SpinPol == "nc"
        Nspin = 4
        spinsize = 1
    end



    Recvecs = 2*pi*inv(Latvecs')

    atompos_unit = Atoms_pos.unit
    if atompos_unit == "frac"
        Gxyz_frac = Atoms_pos.Gxyz_frac
        Gxyz_AU = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            Gxyz_AU[atom] = zeros(Float64, 3)
            x = Gxyz_frac[atom][1]*Latvecs[1,1] + Gxyz_frac[atom][2]*Latvecs[2,1] + Gxyz_frac[atom][3]*Latvecs[3,1]
            y = Gxyz_frac[atom][1]*Latvecs[1,2] + Gxyz_frac[atom][2]*Latvecs[2,2] + Gxyz_frac[atom][3]*Latvecs[3,2]
            z = Gxyz_frac[atom][1]*Latvecs[1,3] + Gxyz_frac[atom][2]*Latvecs[2,3] + Gxyz_frac[atom][3]*Latvecs[3,3]
            Gxyz_AU[atom][1] = x
            Gxyz_AU[atom][2] = y
            Gxyz_AU[atom][3] = z
        end
    else
        Gxyz_AU = Atoms_pos.Gxyz
        Gxyz_frac = Calc_Gxyz_frac(Natom, Gxyz_AU, Recvecs)
    end
    


    Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_orb)
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


    if length(unique(Spe_Symbol)) ≠ Nspecies
        error("not match Atom Symbol and PAO Orbitals, please check Atoms input, or not support empty atom method")
    end



    atom2spe = zeros(Int32, Natom)
    for atom = 1:Natom, spe = 1:Nspecies
        if Atoms_symbol[atom] == Spe_Symbol[spe]
            atom2spe[atom] = spe
        end
    end


    if Hub_U
        for atom = 1:Natom
            spe = atom2spe[atom]
            _, Spe_Num_Basis = get_ialpha_index(Spe_orb[spe])
            Npao = sum(Spe_Num_Basis)
            if length(Hub_U_atom[atom]) ≠ Npao
                @show spe, Spe_orb, Spe_Num_Basis, length(Hub_U_atom[atom]), Npao
                error("please check Hub_U_atom")
            end
        end
    end


    Atoms_cutoff = zeros(Float64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        Atoms_cutoff[atom] = Spe_cutoff[spe]
    end




    if isnothing(Atoms_Nspin)
        Atoms_Nspin = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            spe = atom2spe[atom]
            Atomsymbol = Atoms_symbol[atom]*Spe_extra[spe]
            Atoms_Nspin[atom] = zeros(Float64, 2)
            Atoms_Nspin[atom][1] = Atom_Core_Charge[Atomsymbol]/2
            Atoms_Nspin[atom][2] = Atom_Core_Charge[Atomsymbol]/2
        end
    else
        for atom = 1:Natom
            spe = atom2spe[atom]
            Atomsymbol = Atoms_symbol[atom]*Spe_extra[spe]
            if Atom_Core_Charge[Atomsymbol] ≠ sum(Atoms_Nspin[atom])
                println("valence electron number of $(Atoms_symbol[atom]) is $(Atom_Core_Charge[Atomsymbol])")
                error("please check Atoms_Nspin")
            end
        end

        if xc_type == "LDA"
            error("please check xc_type")
        end
    end


    
    if isnothing(Atoms_Angle)
        Atoms_Angle = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            Atoms_Angle[atom] = zeros(Float64, 4)
        end
    else
        if !SO_switch
            println("ignore initial Angle")
            println("Set all zero putting")
            Atoms_Angle = Vector{Vector{Float64}}(undef, Natom)
            for atom = 1:Natom
                Atoms_Angle[atom] = zeros(Float64, 4)
            end
        end

        if xc_type == "LDA"
            error("please check xc_type")
        end
    end




    if isnothing(Ngrid)
        Ngrid = Calc_Ngrid(Ecut, Latvecs)
    else
        println("Ecut is ignored. Ngrid fixed $(Ngrid)")
    end



    # calculate Latvecs/Ngrid
    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid[1]
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid[2]
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid[3]


    # reciprocal grid 
    gRecvecs = 2*pi*inv(gLatvecs')

    # Find grids overlaping to each atom
    GridVol = abs(det(gLatvecs))

    # Setting the center of unit cell and grids
    Grid_Origin = Calc_Grid_Origin(Natom, Gxyz_AU, gLatvecs, Ngrid, atompos_unit)


    time_rev = ifelse(SpinPol=="nc", false, true)



    symmetry = Get_Symmetry_Spglib(Latvecs, Gxyz_frac, atom2spe)



    pao_file = Vector{String}(undef, Nspecies)
    pspot_file = Vector{String}(undef, Nspecies)
    for spe = 1:Nspecies
        pao_file[spe] = PAO_File_path*Spe_Symbol[spe]*string(Spe_cutoff[spe])*Spe_extra[spe]*".pao" 

        vps_filename = Spe_Symbol[spe]*"_"
        if xc_type ∈ ["LDA", "LSDA"]
            vps_filename = vps_filename*"CA19"*Spe_extra[spe]*".vps"
        elseif xc_type ∈ ["GGA_PBE"]
            vps_filename = vps_filename*"PBE19"*Spe_extra[spe]*".vps"
        else
            println("xc_type is $xc_type")
            error("please check xc_type")
        end

        pspot_file[spe] = VPS_File_path*vps_filename
    end


    filename2, _ = splitext(basename(filename))
    if myrank == 0 && fileout
        Write_CIFfile(filename2, Natom, Latvecs, Gxyz_frac, Atoms_symbol)
        Write_xyzfile(filename2, Natom, Gxyz_AU, Atoms_symbol)
    end


    if restart
        myrank == 0 && RestartFile_check(Natom, Nspin, Latvecs, Gxyz_AU, Grid_Origin, SpinPol, SO_switch, xc_type, filepath)
        myrank == 0 && println("RestartFile_check pass.\n")
        filename2 = restart_filename*"_restart"
    end


    dft_setup = DFT_Setup(
        Natom, Nspecies, Nspin, spinsize, atom2spe,
        Latvecs, Recvecs, gLatvecs, gRecvecs,
        Gxyz_AU, Gxyz_frac, GridVol, Grid_Origin, 
        Atoms_symbol, Atoms_cutoff, Atoms_orb, system, 
        Atoms_Nspin, Atoms_Angle, 
        SpinPol, SO_switch, xc_type, 
        pao_file, pspot_file, Ngrid,
        Mixing_method, SCF_criterion, SCF_max, 
        Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Num_Mixing_Pulay,
        Start_Pulay_SCF, E_Temp, kmesh, symmetry,
        Hub_U, Hub_U_atom, Hub_U_orbpol, Hub_U_occ, Hub_Type, dc_Type,
        time_rev, cal_force, cal_mode, fileout, filename2, restart, filepath2, send_email, verbosity
    )


    if myrank == 0 && verbosity >= 1
        Print_DFT_Setup(dft_setup)
    end
    MPI.Barrier(comm)


    return dft_setup
end
