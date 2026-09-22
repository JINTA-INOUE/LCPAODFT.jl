struct MLWF_Setup
    filepath::String
    material::LCPAO_model
    # mlwf_kpoints::MLWF_KPoints
    WANNUM::Int32
    Guide_index::Vector{Vector{Int32}}
    MLWF_Outer_Window_Bottom::Float64
    MLWF_Outer_Window_Top::Float64
    MLWF_Inner_Window_Bottom::Float64
    MLWF_Inner_Window_Top::Float64
    MLWF_kmesh::Tuple{Int32,Int32,Int32}
    MLWF_Plot_Cube::Union{Vector{Int32},Nothing}
    MLWF_Plot_SuperCells::Union{Vector{Int32},Nothing}
    MAXSHELL::Int32
    write_coef::Bool
    filename::String
    verbosity::Int32
end


function Print_MLWF_Setup(mlwf_setup::MLWF_Setup)

    filepath = mlwf_setup.filepath
    WANNUM = mlwf_setup.WANNUM
    MLWF_Outer_Window_Bottom = mlwf_setup.MLWF_Outer_Window_Bottom
    MLWF_Outer_Window_Top = mlwf_setup.MLWF_Outer_Window_Top
    MLWF_Inner_Window_Bottom = mlwf_setup.MLWF_Inner_Window_Bottom
    MLWF_Inner_Window_Top = mlwf_setup.MLWF_Inner_Window_Top
    MLWF_kmesh = mlwf_setup.MLWF_kmesh
    MAXSHELL = mlwf_setup.MAXSHELL
    filename = mlwf_setup.filename


    println("<MLWF Setup>")
    println("\tscf_inputfilepath: $filepath")
    println("\tOuter_Window: $MLWF_Outer_Window_Bottom, $MLWF_Outer_Window_Top")
    println("\tInner_Window: $MLWF_Inner_Window_Bottom, $MLWF_Inner_Window_Top")
    println("\tWANNUM: $WANNUM")
    println("\tMLWF_kmesh: $MLWF_kmesh")
    println("\tMAXSHELL: $MAXSHELL")
    println("\tfilename: $(filename)")
    println("\t")
end


function MLWF_Setup(
    filepath::AbstractString,
    Guide_index,
    Dis_Energy::Vector{Float64},
    kmesh::Tuple{Signed,Signed,Signed}=(1,1,1);
    MLWF_Plot_Cube = nothing,
    MLWF_Plot_SuperCells = nothing,
    MAXSHELL::Integer = 30,
    filename::String = splitext(basename(filepath))[1],
    write_coef::Bool = false,
    verbosity::Integer = 1)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()
    MKL.set_num_threads(1)

    LCPAODFT.reset_timer!(LCPAODFT.timer)


    if kmesh[1] <= 0 || kmesh[2] <= 0 || kmesh[3] <= 0
        error("please check kmesh")
    end

    MLWF_Outer_Window_Bottom = Dis_Energy[1]
    MLWF_Outer_Window_Top = Dis_Energy[4]
    MLWF_Inner_Window_Bottom = Dis_Energy[2]
    MLWF_Inner_Window_Top = Dis_Energy[3]

    if (MLWF_Outer_Window_Top-MLWF_Outer_Window_Bottom)<=0.0
        error("Error:WF For OUTER window, its top should be higher than bottom.")
    end

    if (MLWF_Inner_Window_Top-MLWF_Inner_Window_Bottom)<0.0
        error("Error:WF For OUTER window, its top should be higher than bottom.")
    end

    if MLWF_Inner_Window_Bottom<MLWF_Outer_Window_Bottom || MLWF_Inner_Window_Top>MLWF_Outer_Window_Top
        error("Error:WF INNER window must be inside of OUTER window.")
    end

    if abs(MLWF_Inner_Window_Bottom-MLWF_Inner_Window_Top) < 1e-6
        MLWF_Inner_Window_Top = MLWF_Outer_Window_Bottom - 999999.0
        MLWF_Inner_Window_Bottom = MLWF_Inner_Window_Top
    end    



    model = select_model(filepath)
    if model == 1
        material = Load_LCPAODFT_model(filepath)
    else
        error("please check filepath.")
    end
    
    if myrank == 0
        Print_LCPAO_model(filepath, material)
    end
    MPI.Barrier(comm)


    SpinPol = material.SpinPol

    GNatom = length(Guide_index)
    WANNUM = 0
    for atom = 1:GNatom
        WANNUM += length(Guide_index[atom])
    end
    WANNUM = ifelse(SpinPol=="nc", 2*WANNUM, WANNUM)


    mlwf_setup = MLWF_Setup(
        filepath, material, WANNUM, Guide_index,
        MLWF_Outer_Window_Bottom, MLWF_Outer_Window_Top, MLWF_Inner_Window_Bottom, MLWF_Inner_Window_Top,
        kmesh, MLWF_Plot_Cube, MLWF_Plot_SuperCells, MAXSHELL,
        write_coef, filename, verbosity
    )


    if myrank == 0 && verbosity >= 1
        Print_MLWF_Setup(mlwf_setup)
    end
    MPI.Barrier(comm)


    return mlwf_setup
end
