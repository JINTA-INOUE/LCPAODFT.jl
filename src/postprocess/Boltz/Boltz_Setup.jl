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
    TDF_dE::Float64         # TDF Enegry size
    Write_TDF::Bool
end


function Print_Boltz_Setup(boltz_setup::Boltz_Setup)

    filename = boltz_setup.filename    
    kmesh = boltz_setup.kmesh
    mat_type = boltz_setup.mat_type
    tau = boltz_setup.tau
    plane_type = boltz_setup.plane_type
    Temp = boltz_setup.Temp
    decomp = boltz_setup.decomp
    muE = boltz_setup.muE
    TDF_Erange = boltz_setup.TDF_Erange
    TDF_dE = boltz_setup.TDF_dE
    Write_TDF = boltz_setup.Write_TDF

    println("<Print_Boltz_Setup>")
    println("\tkmesh : $(kmesh)")
    println("\tmat_type : $(mat_type)")
    println("\ttau : $(tau)")
    println("\tplane_type : $(plane_type)")
    println("\tTemperature : $(Temp)")
    println("\tdecomp : $(decomp)")
    println("\tmuE : $(muE)")
    println("\tTDF_Erange : $(TDF_Erange)")
    println("\tTDF_dE : $(TDF_dE)")
    println("\tWrite_TDF : $(Write_TDF)")
    println("\tfilename : $(filename)")
end


function Boltz_Setup(
    filepath::String,
    kmesh::Tuple{Signed,Signed,Signed},
    TDF_Erange::Vector{Float64},
    _Temp::Union{AbstractFloat,Vector{AbstractFloat}};
    TDF_dE = 0.01,
    tau = 10.0,
    plane_type::Bool = false,
    decomp::Bool = false,
    muE::Union{Float64,Vector{Float64}} = [1000.0],
    Write_TDF::Bool = true,
    filename = nothing)
    

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)


    nthreads = Threads.nthreads()
    # BLAS.set_num_threads(1)
    nblas = BLAS.get_num_threads()
    # println("\t$nthreads threads and $nblas BLAS threads")
    # println("\t$(now())")
    # println("")


    LCPAODFT.reset_timer!(LCPAODFT.timer)


    model = select_model(filepath)


    if nprocs > 1
        error("please run serial.")
    end

    NTemp = length(_Temp)
    Temp = zeros(Float64, NTemp)
    for i = 1:NTemp
        Temp[i] = _Temp[i]
    end


    if decomp
        if muE[1] == 1000.0
            error("please check muE")
        end

        for mu in muE
            if !(TDF_Erange[1] < mu < TDF_Erange[2])
                error("please check muE.")
            end
        end
    end


    if model == 1 && decomp
        error("not support PAO decomp")
    end



    if model == 1
        mat_type = "LCPAO"
        material = Load_LCPAODFT_model(filepath)
        myrank == 0 && Print_LCPAO_model(filepath, material)
    elseif model == 2
        mat_type = "CWF"
        material = Load_CWF_model(filepath)
        myrank == 0 && Print_CWF_model(filepath, material)
    else
        error("please check filepath.")
    end
    MPI.Barrier(comm)


    if isnothing(filename)
        filename = split(filepath, ".")[begin]
    end


    boltz_setup = Boltz_Setup(
        filepath, filename, 
        material, mat_type, 
        kmesh, tau, plane_type, Temp, decomp,
        muE, TDF_Erange, TDF_dE, Write_TDF
    ) 


    if myrank == 0
        Print_Boltz_Setup(boltz_setup)
    end
    MPI.Barrier(comm)
    

    return boltz_setup
end