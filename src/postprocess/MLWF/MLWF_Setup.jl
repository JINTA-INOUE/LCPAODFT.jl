struct MLWF_Setup
    filepath::String
    material::LCPAO_model
    mlwf_kpoints::MLWF_KPoints
    proj2MLWF::Vector{Int32}
    MLWF_Orbs::Vector{String}
    MLWF_Symbol::Vector{String}
    MLWF_Pro_lName::Vector{String}
    MLWF_ProNumBasis::Vector{Int32}
    MLWF_Func_Num::Int32
    MLWF_ProjName::Vector{String}
    MLWF_Outer_Window_Bottom::Float64
    MLWF_Outer_Window_Top::Float64
    MLWF_Inner_Window_Bottom::Float64
    MLWF_Inner_Window_Top::Float64
    MLWF_Initial_Guess::Bool
    MLWF_Num_Kinds_Projectors::Int32
    MLWF_Atom_Cut1::Vector{Float64}
    MLWF_ProjOrbitals::Vector{String}
    MLWF_Total_NumOrbs::Vector{Int32}
    MLWF_Guide::Vector{Vector{Float64}}
    MLWF_Gxyz_AU::Vector{Vector{Float64}}
    MLWF_Gxyz_frac::Vector{Vector{Float64}}
    MLWF_Z_Direction::Array{Float64,2}
    MLWF_X_Direction::Array{Float64,2}
    MLWF_kmesh::Tuple{Int32,Int32,Int32}
    MLWF_MaxShells::Int32
    MLWF_Dis_Mixing_Para::Float64
    MLWF_Dis_Conv_Criterion::Float64
    MLWF_Dis_SCF_Max_Steps::Int32
    MLWF_Min_Scheme::Int32
    MLWF_Minimizing_Max_Steps::Int32
    MLWF_Min_StepLength::Float64
    MLWF_Min_Secant_Steps::Int32
    MLWF_Min_Secant_StepLength::Float64
    MLWF_Min_Conv_Criterion::Float64
    filename::String
    verbosity::Bool
end


function Print_MLWF_Setup(mlwf_setup::MLWF_Setup)

    filepath = mlwf_setup.filepath
    proj2MLWF = mlwf_setup.proj2MLWF
    MLWF_Pro_lName = mlwf_setup.MLWF_Pro_lName
    MLWF_ProjOrbitals = mlwf_setup.MLWF_ProjOrbitals
    MLWF_ProNumBasis = mlwf_setup.MLWF_ProNumBasis
    MLWF_Symbol = mlwf_setup.MLWF_Symbol
    MLWF_Func_Num = mlwf_setup.MLWF_Func_Num
    MLWF_ProjName = mlwf_setup.MLWF_ProjName
    MLWF_Outer_Window_Bottom = mlwf_setup.MLWF_Outer_Window_Bottom
    MLWF_Outer_Window_Top = mlwf_setup.MLWF_Outer_Window_Top
    MLWF_Inner_Window_Bottom = mlwf_setup.MLWF_Inner_Window_Bottom
    MLWF_Inner_Window_Top = mlwf_setup.MLWF_Inner_Window_Top
    MLWF_Initial_Guess = mlwf_setup.MLWF_Initial_Guess
    MLWF_Num_Kinds_Projectors = mlwf_setup.MLWF_Num_Kinds_Projectors
    MLWF_Guide = mlwf_setup.MLWF_Guide
    MLWF_Gxyz_AU = mlwf_setup.MLWF_Gxyz_AU
    MLWF_Gxyz_frac = mlwf_setup.MLWF_Gxyz_frac
    MLWF_Z_Direction = mlwf_setup.MLWF_Z_Direction
    MLWF_X_Direction = mlwf_setup.MLWF_X_Direction
    MLWF_kmesh = mlwf_setup.MLWF_kmesh
    MLWF_MaxShells = mlwf_setup.MLWF_MaxShells
    MLWF_Dis_SCF_Max_Steps = mlwf_setup.MLWF_Dis_SCF_Max_Steps
    MLWF_Dis_Conv_Criterion = mlwf_setup.MLWF_Dis_Conv_Criterion
    MLWF_Dis_Mixing_Para = mlwf_setup.MLWF_Dis_Mixing_Para
    MLWF_Min_Secant_Steps = mlwf_setup.MLWF_Min_Secant_Steps
    MLWF_Min_StepLength = mlwf_setup.MLWF_Min_StepLength
    MLWF_Min_Scheme = mlwf_setup.MLWF_Min_Scheme
    MLWF_Minimizing_Max_Steps = mlwf_setup.MLWF_Minimizing_Max_Steps
    MLWF_Min_Conv_Criterion = mlwf_setup.MLWF_Min_Conv_Criterion
    MLWF_Min_Secant_StepLength = mlwf_setup.MLWF_Min_Secant_StepLength
    filename = mlwf_setup.filename


    println("<MLWF Setup>")
    println("\tscf_inputfilepath: $filepath")
    println("\tOuter_Window: $MLWF_Outer_Window_Bottom, $MLWF_Outer_Window_Top")
    println("\tInner_Window: $MLWF_Inner_Window_Bottom, $MLWF_Inner_Window_Top")
    println("")
    println("\tMLWF_Initial_Guess: $MLWF_Initial_Guess")
    if MLWF_Initial_Guess
        println("\tMLWF_ProjOrbitals: $MLWF_ProjOrbitals")
        println("\tMLWF_Num_Kinds_Projectors: $MLWF_Num_Kinds_Projectors")
        for p = 1:MLWF_Num_Kinds_Projectors
            @printf("\t%d\t%s\t%d%s\t%5.6f  %5.6f  %5.6f\n", p, MLWF_Symbol[p], MLWF_ProNumBasis[p], MLWF_Pro_lName[p], MLWF_Gxyz_AU[p][1], MLWF_Gxyz_AU[p][2], MLWF_Gxyz_AU[p][3])
        end

        println("\tMLWF_Func_Num: $MLWF_Func_Num")

        if MLWF_Func_Num ≠ MLWF_Num_Kinds_Projectors
            for proj = 1:MLWF_Func_Num
                ist = proj2MLWF[proj]
                @printf("\t\t%d\t%d\t%s\t%d  %s\t%5.6f  %5.6f  %5.6f\n", proj, ist, MLWF_Symbol[ist], MLWF_ProNumBasis[ist], MLWF_ProjName[proj], MLWF_Guide[proj][1], MLWF_Guide[proj][2], MLWF_Guide[proj][3])
            end
        end
    end
    println("")

    @show MLWF_Z_Direction
    @show MLWF_X_Direction
    println("\tMLWF_kmesh: $MLWF_kmesh")
    println("\tMLWF_MaxShells: $MLWF_MaxShells")
    println("")
    println("<MLWF Disentangle Parameters>")
    println("\tMLWF_Dis_SCF_Max_Steps:  $MLWF_Dis_SCF_Max_Steps")
    println("\tMLWF_Dis_Conv_Criterion:  $MLWF_Dis_Conv_Criterion")
    println("\tMLWF_Dis_Mixing_Para:  $MLWF_Dis_Mixing_Para")
    println("")
    println("<MLWF Minimizing Parameters>")
    println("\tMLWF_Min_Secant_Steps:  $MLWF_Min_Secant_Steps")
    println("\tMLWF_Min_StepLength:  $MLWF_Min_StepLength")
    println("\tMLWF_Min_Scheme:  $MLWF_Min_Scheme")
    println("\tMLWF_Min_Secant_StepLength:  $MLWF_Min_Secant_StepLength")
    println("\tMLWF_Min_Conv_Criterion:  $MLWF_Min_Conv_Criterion")
    println("\tfilename: $(filename)")
    println("\t")
end


function MLWF_Setup(
    filepath::String,
    MLWF_Orbs::Vector{String},
    MLWF_Pos::Atompos,
    kmesh::Tuple{Int64,Int64,Int64}, 
    Eouter::Vector{Float64},
    Einner::Vector{Float64};
    MLWF_Angle = nothing,
    MLWF_MaxShells = 30,
    MLWF_Minimizing_Max_Steps = 200,
    MLWF_Dis_Mixing_Para = 0.5,
    MLWF_Dis_Conv_Criterion = 1e-10,
    MLWF_Dis_SCF_Max_Steps = 10000,
    MLWF_Min_Scheme = 2,
    MLWF_Min_StepLength = 2.0,
    MLWF_Min_Secant_Steps = 5,
    MLWF_Min_Secant_StepLength = 2.0,
    MLWF_Min_Conv_Criterion = 1e8,
    MLWF_Initial_Guess::Bool = true,
    filename::String = splitext(basename(filepath))[1],
    verbosity::Bool = true)


    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    BLAS.set_num_threads(1)
    nthreads = Threads.nthreads()
    nblas = BLAS.get_num_threads()

    LCPAODFT.reset_timer!(LCPAODFT.timer)


    if kmesh[1] <= 0 || kmesh[2] <= 0 || kmesh[3] <= 0
        error("please check kmesh")
    end

    if length(MLWF_Orbs) ≠ MLWF_Pos.Natom
        error("please check input")
    end



    MLWF_Outer_Window_Bottom, MLWF_Outer_Window_Top = Eouter
    MLWF_Inner_Window_Bottom, MLWF_Inner_Window_Top = Einner

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





    MLWF_Num_Kinds_Projectors = length(MLWF_Orbs)
    MLWF_Pro_lName, MLWF_ProNumBasis, MLWF_Symbol, MLWF_EachNum, MLWF_Func_Num = Set_MLWF_Projector(MLWF_Orbs)
    @show MLWF_Pro_lName 
    @show MLWF_ProNumBasis
    @show MLWF_Symbol
    @show MLWF_EachNum
    @show MLWF_Func_Num
    
    MLWF_ProjName = Set_MLWF_ProjName(MLWF_Num_Kinds_Projectors, MLWF_Pro_lName)
    @show MLWF_ProjName
    



    if isnothing(MLWF_Angle)
        _MLWF_Angle = Vector{Vector{Float64}}(undef, MLWF_Num_Kinds_Projectors)
        for p = 1:MLWF_Num_Kinds_Projectors
            _MLWF_Angle[p] = [0.0, 0.0, 1.0, 1.0, 0.0,  0.0]
        end
    else
        _MLWF_Angle = MLWF_Angle
    end


    MLWF_Spe_Symbol, MLWF_Spe_cutoff, MLWF_Spe_orb, MLWF_Spe_extra = Get_Atoms_data(MLWF_Orbs)
    @show MLWF_Spe_Symbol, MLWF_Spe_cutoff, MLWF_Spe_orb, MLWF_Spe_extra


    MLWF_Atom_Cut1 = MLWF_Spe_cutoff
    MLWF_ProjOrbitals = Vector{String}(undef, MLWF_Num_Kinds_Projectors)
    for proj = 1:MLWF_Num_Kinds_Projectors
        # MLWF_ProjOrbitals[proj] = MLWF_Spe_Symbol[proj]*string(MLWF_Atom_Cut1[proj])*string(MLWF_Spe_extra[proj])*"-s1p1d1f1"
        MLWF_ProjOrbitals[proj] = MLWF_Spe_Symbol[proj]*string(MLWF_Atom_Cut1[proj])*string(MLWF_Spe_extra[proj])*"-s2p2d1f1"
    end

    MLWF_Total_NumOrbs = zeros(Int32, MLWF_Num_Kinds_Projectors)
    for proj = 1:MLWF_Num_Kinds_Projectors
        MLWF_Total_NumOrbs[proj] = 1+3+5+7
    end

    
    MLWF_Z_Direction = zeros(Float64, MLWF_Num_Kinds_Projectors, 3)
    MLWF_X_Direction = zeros(Float64, MLWF_Num_Kinds_Projectors, 3)
    for i = 1:MLWF_Num_Kinds_Projectors
        @. MLWF_Z_Direction[i,:] = _MLWF_Angle[i][1:3]
        @. MLWF_X_Direction[i,:] = _MLWF_Angle[i][4:6]
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


    Latvecs = material.Latvecs
    Recvecs = material.Recvecs


    MLWF_Natom = MLWF_Pos.Natom
    MLWF_atompos_unit = MLWF_Pos.unit
    if MLWF_atompos_unit == "frac"
        MLWF_Gxyz_frac = MLWF_Pos.Gxyz_frac
        MLWF_Gxyz_AU = Vector{Vector{Float64}}(undef, MLWF_Natom)
        for atom = 1:MLWF_Natom
            MLWF_Gxyz_AU[atom] = zeros(Float64, 3)
            x = MLWF_Gxyz_frac[atom][1]*Latvecs[1,1] + MLWF_Gxyz_frac[atom][2]*Latvecs[2,1] + MLWF_Gxyz_frac[atom][3]*Latvecs[3,1]
            y = MLWF_Gxyz_frac[atom][1]*Latvecs[1,2] + MLWF_Gxyz_frac[atom][2]*Latvecs[2,2] + MLWF_Gxyz_frac[atom][3]*Latvecs[3,2]
            z = MLWF_Gxyz_frac[atom][1]*Latvecs[1,3] + MLWF_Gxyz_frac[atom][2]*Latvecs[2,3] + MLWF_Gxyz_frac[atom][3]*Latvecs[3,3]
            MLWF_Gxyz_AU[atom][1] = x
            MLWF_Gxyz_AU[atom][2] = y
            MLWF_Gxyz_AU[atom][3] = z
        end
    else
        MLWF_Gxyz_frac = Calc_Gxyz_frac(MLWF_Natom, MLWF_Gxyz_AU, Recvecs)
    end


    proj2MLWF = zeros(Int32, MLWF_Func_Num)
    proj = 1
    for ist = 1:MLWF_Num_Kinds_Projectors, _ = 1:MLWF_EachNum[ist]
        proj2MLWF[proj] = ist
        proj += 1
    end

    MLWF_Guide = Vector{Vector{Float64}}(undef, MLWF_Func_Num)
    for ist = 1:MLWF_Func_Num
        MLWF_Guide[ist] = zeros(Float64, 3)
    end

    for proj = 1:MLWF_Func_Num
        ist = proj2MLWF[proj]
        @. MLWF_Guide[proj] = MLWF_Gxyz_AU[ist]
    end



    myrank == 0 && println("<Set_MLWF_kgrid>")
    mlwf_kpoints = Set_MLWF_kgrid(Latvecs, kmesh, MLWF_MaxShells)
    myrank == 0 && Print_MLWF_kpoints(mlwf_kpoints)
    myrank == 0 && println("")



    mlwf_setup = MLWF_Setup(
        filepath, material, mlwf_kpoints,
        proj2MLWF, MLWF_Orbs, MLWF_Symbol, MLWF_Pro_lName, MLWF_ProNumBasis,
        MLWF_Func_Num, MLWF_ProjName,
        MLWF_Outer_Window_Bottom, MLWF_Outer_Window_Top, MLWF_Inner_Window_Bottom, MLWF_Inner_Window_Top,
        MLWF_Initial_Guess, 
        MLWF_Num_Kinds_Projectors, MLWF_Atom_Cut1, MLWF_ProjOrbitals, MLWF_Total_NumOrbs,
        MLWF_Guide, MLWF_Gxyz_AU, MLWF_Gxyz_frac, 
        MLWF_Z_Direction, MLWF_X_Direction,
        kmesh, MLWF_MaxShells,
        MLWF_Dis_Mixing_Para, MLWF_Dis_Conv_Criterion, MLWF_Dis_SCF_Max_Steps,
        MLWF_Min_Scheme, MLWF_Minimizing_Max_Steps, MLWF_Min_StepLength, MLWF_Min_Secant_Steps, MLWF_Min_Secant_StepLength, MLWF_Min_Conv_Criterion,
        filename, verbosity
    )


    if myrank == 0 && verbosity
        Print_MLWF_Setup(mlwf_setup)
    end


    return mlwf_setup
end
