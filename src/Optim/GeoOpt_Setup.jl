struct GeoOpt_Setup
    dft_setup::DFT_Setup
    Natom::Int32
    Atoms_symbol::Vector{String}
    Geo_Opt_Max::Int32
    Geo_Opt_criterion::Float64
    M_GDIIS_HISTORY::Int32
    OptStartDIIS::Int32
    OptEveryDIIS::Int32
    Extra_CHistory::Int32
    atom_Fixed_XYZ::Vector{Vector{Int32}}
end


function Print_GeoOpt_Setup(geoopt_setup::GeoOpt_Setup)

    Natom = geoopt_setup.Natom
    Atoms_symbol = geoopt_setup.Atoms_symbol
    Geo_Opt_Max = geoopt_setup.Geo_Opt_Max
    Geo_Opt_criterion = geoopt_setup.Geo_Opt_criterion
    M_GDIIS_HISTORY = geoopt_setup.M_GDIIS_HISTORY
    OptStartDIIS = geoopt_setup.OptStartDIIS
    OptEveryDIIS = geoopt_setup.OptEveryDIIS
    Extra_CHistory = geoopt_setup.Extra_CHistory
    atom_Fixed_XYZ = geoopt_setup.atom_Fixed_XYZ


    println("<Print_GeoOpt_Setup>")
    println("\tGeo_Opt_Max : $(Geo_Opt_Max)")
    println("\tGeo_Opt_criterion : $(Geo_Opt_criterion)")
    println("\tM_GDIIS_HISTORY : $(M_GDIIS_HISTORY)")
    println("\tOptStartDIIS : $(OptStartDIIS)")
    println("\tOptEveryDIIS : $(OptEveryDIIS)")
    println("\tExtra_CHistory : $(Extra_CHistory)")

    println("")
    println("Fix atom positions for Geometry optimization")
    println("\tatom\tAtom Name\t x\t y\t z")
    for atom = 1:Natom
        print("\t$(atom)\t$(Atoms_symbol[atom])\t")
        for xyz = 1:3
            if atom_Fixed_XYZ[atom][xyz] == 0
                print("\trelax  ")
            elseif atom_Fixed_XYZ[atom][xyz] == 1
                print("\tfix  ")
            else
                error("please check atom_Fixed_XYZ")
            end
        end
        print("\n")
    end
end


function GeoOpt_Setup(
    dft_setup::DFT_Setup;
    atom_Fixed_XYZ::Union{Vector{Vector{Int64}},Nothing} = nothing,
    Geo_Opt_Max = 1,
    Geo_Opt_criterion = 0.0003,
    M_GDIIS_HISTORY = 2,
    OptStartDIIS = 1,
    OptEveryDIIS = 200000,
    Extra_CHistory = 3)


    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)


    Nspin = dft_setup.Nspin
    Natom = dft_setup.Natom
    Atoms_symbol = dft_setup.Atoms_symbol
    Ngrid = dft_setup.Ngrid
    filename = dft_setup.filename
    restart = dft_setup.restart
    verbosity = dft_setup.verbosity

    if restart
        error("Geometrix Optimization restart not support.")
    end

    if Geo_Opt_Max <= 0
        error("please check Geo_Opt_Max\n")
    end

    if Geo_Opt_criterion < 0.0 || Geo_Opt_criterion > 1.0
        error("please check Geo_Opt_criterion\n")
    end

    if !isnothing(atom_Fixed_XYZ)    
        atom_Fixed_XYZ2 = atom_Fixed_XYZ
        if length(atom_Fixed_XYZ2) ≠ Natom
            error("please check atom_Fixed_XYZ\n")
        end

        for atom = 1:Natom
            if length(atom_Fixed_XYZ2[atom]) ≠ 3
                error("please check atom_Fixed_XYZ\n")
            end

            for xyz = 1:3
                if atom_Fixed_XYZ2[atom][xyz] ≠ 0 && atom_Fixed_XYZ2[atom][xyz] ≠ 1
                    error("please check atom_Fixed_XYZ\n")
                end
            end
        end
    else
        atom_Fixed_XYZ2 = Vector{Vector{Int64}}(undef, Natom)
        for atom = 1:Natom
            atom_Fixed_XYZ2[atom] = zeros(Int64, 3)
        end
    end

    if M_GDIIS_HISTORY > 19
        error("please check M_GDIIS_HISTORY")
    end


    if myrank == 0
        work_dirname = pwd()*"/"*filename*"_work"
        mkpath(work_dirname)
        
        NN = prod(Ngrid)
        tmp_array = Vector{Vector{Float64}}(undef, Nspin)
        for spin = 1:Nspin
            tmp_array[spin] = zeros(Float64, NN)
        end
        for i = 1:Extra_CHistory
            file = work_dirname*"/"*filename*"_work_rho$i.jld2"
            Generate_rhoFile(file, Ngrid, tmp_array)
        end
    end


    geoopt_setup = GeoOpt_Setup(
        dft_setup, 
        Natom, Atoms_symbol, 
        Geo_Opt_Max, Geo_Opt_criterion, M_GDIIS_HISTORY, 
        OptStartDIIS, OptEveryDIIS, Extra_CHistory,
        atom_Fixed_XYZ2)

    if myrank == 0 && verbosity >= 1
        Print_GeoOpt_Setup(geoopt_setup)
    end

    
    return geoopt_setup
end
