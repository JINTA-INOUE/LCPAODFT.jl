struct Boltz_Setup
    filepath::String
    filename::String
    material::Union{CWF_model,LCPAO_model}
    mat_type::String
    kmesh::Tuple{Int32,Int32,Int32}
    tau::Float64
    plane_type::Bool
    Temp::Vector{Float64}
    decomp::Bool
    muE::Vector{Float64}
    TDF_Erange::Vector{Float64}
    TDF_dE::Float64
    Write_TDF::Bool
end


function Print_Boltz_Setup(boltz_setup::Boltz_Setup)
    println("<Print_Boltz_Setup>")
    println("\tkmesh : $(boltz_setup.kmesh)")
    println("\tmat_type : $(boltz_setup.mat_type)")
    println("\ttau : $(boltz_setup.tau)")
    println("\tplane_type : $(boltz_setup.plane_type)")
    println("\tTemperature : $(boltz_setup.Temp)")
    println("\tdecomp : $(boltz_setup.decomp)")
    println("\tmuE : $(boltz_setup.muE)")
    println("\tTDF_Erange : $(boltz_setup.TDF_Erange)")
    println("\tTDF_dE : $(boltz_setup.TDF_dE)")
    println("\tWrite_TDF : $(boltz_setup.Write_TDF)")
    println("\tfilename : $(boltz_setup.filename)")
end


function Boltz_Setup(
    filepath::String,
    kmesh::Tuple{Signed,Signed,Signed},
    TDF_Erange::Vector{Float64},
    _Temp::Union{Real,AbstractVector{<:Real}};
    TDF_dE=0.01,
    tau=10.0,
    plane_type::Bool=false,
    decomp::Bool=false,
    muE::Union{Real,AbstractVector{<:Real}}=[1000.0],
    Write_TDF::Bool=false,
    filename=nothing)

    Threads.nthreads() == 1 || error(
        "MPI-flat Boltz requires exactly one Julia thread per MPI process; " *
        "start Julia with --threads=1",
    )
    provided_thread_level = MPI.Init(; threadlevel=:single)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    BLAS.set_num_threads(1)
    LCPAODFT.reset_timer!(LCPAODFT.timer)

    all(>(0), kmesh) || error("all kmesh dimensions must be positive")
    length(TDF_Erange) == 2 || error("TDF_Erange must contain its lower and upper limits")
    TDF_Erange[1] < TDF_Erange[2] || error("TDF_Erange must be strictly increasing")
    TDF_dE > 0 || error("TDF_dE must be positive")
    tau > 0 || error("tau must be positive")
    plane_type && kmesh[3] != 1 && error("plane_type=true requires kmesh[3] == 1")

    Temp = _Temp isa AbstractVector ? Float64.(_Temp) : [Float64(_Temp)]
    isempty(Temp) && error("at least one temperature is required")
    all(>(0), Temp) || error("all temperatures must be positive")
    mu_values = muE isa AbstractVector ? Float64.(muE) : [Float64(muE)]

    if decomp
        mu_values[1] == 1000.0 && error("please check muE")
        for mu in mu_values
            TDF_Erange[1] < mu < TDF_Erange[2] || error("each muE value must lie inside TDF_Erange")
        end
    end

    model = select_model(filepath)
    model == 1 && error("transport is not supported for LCPAO models")

    if myrank == 0
        println("<Boltz MPI-flat configuration>")
        println("\t$nprocs MPI processes × 1 Julia thread")
        println("\t$(BLAS.get_num_threads()) BLAS thread per process")
        println("\tMPI thread level: $provided_thread_level")
    end

    if model == 2
        mat_type = "CWF"
        material = Load_CWF_model(filepath)
        myrank == 0 && Print_CWF_model(filepath, material)
    else
        error("unsupported model file: $filepath")
    end
    MPI.Barrier(comm)

    output_filename = isnothing(filename) ? splitext(filepath)[1] : String(filename)
    boltz_setup = Boltz_Setup(
        filepath, output_filename, material, mat_type,
        Tuple(Int32.(kmesh)), Float64(tau), plane_type, Temp, decomp,
        mu_values, TDF_Erange, Float64(TDF_dE), Write_TDF,
    )

    myrank == 0 && Print_Boltz_Setup(boltz_setup)
    MPI.Barrier(comm)

    
    return boltz_setup
end
