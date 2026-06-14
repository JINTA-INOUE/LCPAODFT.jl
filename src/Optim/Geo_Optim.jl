mutable struct Geo_Optim
    Geo_Opt_Max::Int32
    Geo_Opt_criterion::Float64
    Geo_Opt_convergence::Bool
    Geo_Opt_Conv_iter::Int32
    M_GDIIS_HISTORY::Int32
    OptStartDIIS::Int32
    OptEveryDIIS::Int32
    Gxyz::Vector{Vector{Float64}}
    Gxyz_fix::Vector{Vector{Float64}}
    GxyzHisIn::Vector{Vector{Vector{Float64}}}
    GxyzHisR::Vector{Vector{Vector{Float64}}}
    Geo_Opt_Max_Force::Vector{Float64}
    atom_Fixed_XYZ::Vector{Vector{Int32}}
    Utot0::Float64
    scaling_factor::Float64
    SD_scaling_user::Float64
    SD_scaling::Float64
    UtotHis::Vector{Float64}
    Geo_Hessian::Matrix{Float64}
    local_iter::Int32
    SD_iter::Int32
    GDIIS_iter::Int32
    flag::Int32
    Every_iter::Int32
    Natom::Int32
    Atoms_Core_Charge::Vector{Float64}
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    ncn::Vector{Vector{Int32}}
    Dis::Vector{Vector{Float64}}
    atv::Vector{Vector{Float64}}
    verbosity::Int32
end


function Init_Geo_Optim(Atoms_symbol::Vector{String}, geoopt_setup::GeoOpt_Setup)

    Natom = length(Atoms_symbol)
    Gxyz_fix = geoopt_setup.dft_setup.Gxyz
    Gxyz = [[0.0, 0.0, 0.0]]
    FNAN = [0]
    natn = [[0]]
    ncn = [[0]]
    Dis = [[0.0]]
    atv = [[0.0, 0.0, 0.0]]

    Geo_Opt_convergence = false
    Geo_Opt_Max = geoopt_setup.Geo_Opt_Max
    Geo_Opt_Conv_iter = Geo_Opt_Max
    Geo_Opt_criterion = geoopt_setup.Geo_Opt_criterion
    M_GDIIS_HISTORY = geoopt_setup.M_GDIIS_HISTORY
    OptStartDIIS = geoopt_setup.OptStartDIIS
    OptEveryDIIS = geoopt_setup.OptEveryDIIS
    atom_Fixed_XYZ = geoopt_setup.atom_Fixed_XYZ
    verbosity = geoopt_setup.dft_setup.verbosity


    GxyzHisIn = Vector{Vector{Vector{Float64}}}(undef, M_GDIIS_HISTORY+1)
    GxyzHisR = Vector{Vector{Vector{Float64}}}(undef, M_GDIIS_HISTORY+1)
    for i = 1:M_GDIIS_HISTORY+1
        GxyzHisIn[i] = Vector{Vector{Float64}}(undef, Natom)
        GxyzHisR[i] = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            GxyzHisIn[i][atom] = zeros(Float64, 4)
            GxyzHisR[i][atom] = zeros(Float64, 4)
        end
    end
    Geo_Opt_Max_Force = zeros(Float64, Geo_Opt_Max+3)

    Utot0 = 0.0
    scaling_factor = 0.0
    SD_scaling_user = 0.0
    SD_scaling = 0.0
    UtotHis = zeros(Float64, 10)
    Geo_Hessian = zeros(Float64, 3*Natom+2, 3*Natom+2)
    local_iter = 1
    SD_iter = 0
    GDIIS_iter = 0
    flag = 0
    GDIIS_EF = 0
    Every_iter = 0


    Atoms_Zcharge = zeros(Float64, Natom)
    for atom = 1:Natom
        Atoms_Zcharge[atom] = Atom_Znumber[Atoms_symbol[atom]]
    end


    return Geo_Optim(
        Geo_Opt_Max, Geo_Opt_criterion, Geo_Opt_convergence, Geo_Opt_Conv_iter,
        M_GDIIS_HISTORY, OptStartDIIS, OptEveryDIIS,
        Gxyz, Gxyz_fix, GxyzHisIn, GxyzHisR, Geo_Opt_Max_Force, atom_Fixed_XYZ,
        Utot0, scaling_factor, SD_scaling_user, SD_scaling, UtotHis, Geo_Hessian,
        local_iter, SD_iter, GDIIS_iter, flag, Every_iter,
        Natom, Atoms_Zcharge, 
        FNAN, natn, ncn, Dis, atv, verbosity
    )
end



function Set_Geo_Optim!(geo_optim::Geo_Optim, system_grid::System_Grid)

    geo_optim.Natom = system_grid.Natom
    geo_optim.Gxyz = system_grid.Gxyz
    geo_optim.FNAN = system_grid.FNAN
    geo_optim.natn = system_grid.natn
    geo_optim.ncn = system_grid.ncn
    geo_optim.Dis = system_grid.Dis
    geo_optim.atv = system_grid.atv
end


function Geo_Optim!(iter, geo_optim::Geo_Optim, energy::Energy, force::Force)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Natom = geo_optim.Natom
    M_GDIIS_HISTORY = geo_optim.M_GDIIS_HISTORY
    OptStartDIIS = geo_optim.OptStartDIIS
    OptEveryDIIS = geo_optim.OptEveryDIIS
    Gxyz = geo_optim.Gxyz
    GxyzHisIn = geo_optim.GxyzHisIn
    GxyzHisR = geo_optim.GxyzHisR
    atom_Fixed_XYZ = geo_optim.atom_Fixed_XYZ
    flag = geo_optim.flag
    Every_iter = geo_optim.Every_iter
    local_iter = geo_optim.local_iter
    GDIIS_iter = geo_optim.GDIIS_iter
    SD_iter = geo_optim.SD_iter
    Etot = energy.Etot
    ForceAll = force.ForceAll


    myrank == 0 && @show flag, Every_iter, local_iter, GDIIS_iter, SD_iter


    Every_iter = OptEveryDIIS
    diis_iter = ifelse(iter < M_GDIIS_HISTORY, iter, M_GDIIS_HISTORY)

    if iter < OptStartDIIS
        flag = 0
    elseif iter == OptStartDIIS
        flag = 1
        GDIIS_iter += 1
    elseif flag == 0
        SD_iter += 1
    else flag == 1
        GDIIS_iter += 1
    end


    if flag == 0
        if SD_iter == 1
            Steepest_Descent(iter, geo_optim, Etot, ForceAll, 0)
        else
            Steepest_Descent(iter, geo_optim, Etot, ForceAll, 1)
        end


        for i = diis_iter-1:-1:1
            for atom = 1:Natom, k = 1:3
                GxyzHisIn[i+1][atom][k] = GxyzHisIn[i][atom][k]
                GxyzHisR[i+1][atom][k] = GxyzHisR[i][atom][k]
            end
        end

        for atom = 1:Natom, xyz = 1:3
            if atom_Fixed_XYZ[atom][xyz] == 0
                GxyzHisIn[1][atom][xyz] = Gxyz[atom][xyz]
                GxyzHisR[1][atom][xyz] = ForceAll[atom,xyz]
            else
                GxyzHisIn[1][atom][xyz] = Gxyz[atom][xyz]
                GxyzHisR[1][atom][xyz] = 0.0
            end
        end

        local_iter = 1
    else
        GDIIS_EF(local_iter, geo_optim, Etot, ForceAll)
        local_iter += 1
    end


    if Every_iter <= SD_iter
        flag = 1
        SD_iter = 0
        GDIIS_iter = 0
    elseif Every_iter <= GDIIS_iter
        flag = 0
        SD_iter = 0
        GDIIS_iter = 0
    end
    

    geo_optim.GxyzHisIn = GxyzHisIn
    geo_optim.GxyzHisR = GxyzHisR
    geo_optim.flag = flag
    geo_optim.Every_iter = Every_iter
    geo_optim.local_iter = local_iter
    geo_optim.GDIIS_iter = GDIIS_iter
    geo_optim.SD_iter = SD_iter
end


function Steepest_Descent(GeoOpt_iter, geo_optim::Geo_Optim, Etot, ForceAll, SD_scaling_flag)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Natom = geo_optim.Natom
    Gxyz = geo_optim.Gxyz
    Geo_Opt_Max = geo_optim.Geo_Opt_Max
    Geo_Opt_convergence = geo_optim.Geo_Opt_convergence
    Geo_Opt_criterion = geo_optim.Geo_Opt_criterion
    Geo_Opt_Max_Force = geo_optim.Geo_Opt_Max_Force
    atom_Fixed_XYZ = geo_optim.atom_Fixed_XYZ
    SD_scaling_user = geo_optim.SD_scaling_user
    SD_scaling = geo_optim.SD_scaling
    UtotHis = geo_optim.UtotHis
    verbosity = geo_optim.verbosity


    Max_ForceNorm = 0.0
    for atom = 1:Natom
        tmp = 0.0
        for xyz = 1:3
            if atom_Fixed_XYZ[atom][xyz] == 0
                tmp += ForceAll[atom,xyz]*ForceAll[atom,xyz]
            end
        end
        Max_ForceNorm = max(sqrt(tmp), Max_ForceNorm)
    end
    Geo_Opt_Max_Force[GeoOpt_iter] = Max_ForceNorm


    if Max_ForceNorm < Geo_Opt_criterion
        Geo_Opt_convergence = true
        myrank == 0 && println("")
        myrank == 0 && println("\nThe geometry optimization was achieved.\n")
    end
    MPI.Barrier(comm)


    BohrR = 0.529177249
    unified_atomic_mass_unit = 1.660538921
    electron_mass = 0.000910938291
    Wscale = unified_atomic_mass_unit/electron_mass
    dt = 41.3411*2.0
    SD_init = dt*dt/Wscale
    SD_max = SD_init*15.0
    SD_min = SD_init*0.04
    Atom_W = 12.0

   
    
    if GeoOpt_iter == 1 || SD_scaling_flag == 0
        SD_scaling_user = Max_ForceNorm/BohrR/1.5
        SD_scaling = SD_scaling_user/(Max_ForceNorm+1e-10)

        SD_scaling = ifelse(SD_scaling>SD_max, SD_max, SD_scaling)
        SD_scaling = ifelse(SD_scaling<SD_min, SD_min, SD_scaling)
    else
        if UtotHis[1] < Etot
            SD_scaling = SD_scaling/2
        elseif UtotHis[1]<UtotHis[2] && Etot<UtotHis[1] && GeoOpt_iter%4==1
            SD_scaling = SD_scaling*2.5
        end

        SD_scaling = ifelse(SD_scaling>SD_max, SD_max, SD_scaling)
        SD_scaling = ifelse(SD_scaling<SD_min, SD_min, SD_scaling)

        UtotHis[5] = UtotHis[4]
        UtotHis[4] = UtotHis[3]
        UtotHis[3] = UtotHis[2]
        UtotHis[2] = UtotHis[1]
        UtotHis[1] = Etot
    end

    if myrank == 0
        println("")
        println("<Steepest_Descent>  SD_scaling=$(SD_scaling)")
    end
    MPI.Barrier(comm)


    Criterion_Max_Step = 0.2
    if !Geo_Opt_convergence && GeoOpt_iter ≠ Geo_Opt_Max
        if Criterion_Max_Step < Max_ForceNorm*SD_scaling
            scale = Criterion_Max_Step/Max_ForceNorm/SD_scaling
        else
            scale = 1.0
        end


        for atom = 1:Natom, xyz = 1:3
            if atom_Fixed_XYZ[atom][xyz] == 0
                Gxyz[atom][xyz] = Gxyz[atom][xyz] - scale*SD_scaling*ForceAll[atom,xyz]
            end
        end
    end


    if myrank == 0 && verbosity >= 1
        @printf("<Steepest_Descent>  |Maximum force| (Hartree/Bohr) =%15.12f\n", Max_ForceNorm)
        @printf("<Steepest_Descent>  Criterion       (Hartree/Bohr) =%15.12f\n", Geo_Opt_criterion)

        @printf("\n")
        @printf("<ForceAll>\n")
        for atom = 1:Natom
            @printf("  Fxyz(a.u.)  atom = %3d  %15.12f  %15.12f  %15.12f\n", atom, ForceAll[atom,1], ForceAll[atom,2], ForceAll[atom,3])
        end

        @printf("<Atomic positions>  New atomic positions using Steepest_Descent\n")
        for atom = 1:Natom
            @printf("  XYZ(ang)    atom = %3d  %15.12f  %15.12f  %15.12f\n",  atom, Gxyz[atom][1]/Ang_to_bohr, Gxyz[atom][2]/Ang_to_bohr, Gxyz[atom][3]/Ang_to_bohr)
        end
    end
    MPI.Barrier(comm)


    geo_optim.Gxyz = Gxyz
    geo_optim.Geo_Opt_convergence = Geo_Opt_convergence
    geo_optim.SD_scaling_user = SD_scaling_user
    geo_optim.SD_scaling = SD_scaling
    geo_optim.UtotHis = UtotHis
    geo_optim.Geo_Opt_Max_Force = Geo_Opt_Max_Force
end


function GDIIS_EF(GeoOpt_iter, geo_optim::Geo_Optim, Etot, ForceAll)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Natom = geo_optim.Natom
    Geo_Opt_Max = geo_optim.Geo_Opt_Max
    Geo_Opt_convergence = geo_optim.Geo_Opt_convergence
    Geo_Opt_criterion = geo_optim.Geo_Opt_criterion
    M_GDIIS_HISTORY = geo_optim.M_GDIIS_HISTORY
    OptStartDIIS = geo_optim.OptStartDIIS
    Gxyz = geo_optim.Gxyz
    Gxyz_fix = geo_optim.Gxyz_fix
    GxyzHisIn = geo_optim.GxyzHisIn
    GxyzHisR = geo_optim.GxyzHisR
    Geo_Opt_Max_Force = geo_optim.Geo_Opt_Max_Force
    atom_Fixed_XYZ = geo_optim.atom_Fixed_XYZ
    Utot0 = geo_optim.Utot0
    scaling_factor = geo_optim.scaling_factor
    SD_scaling = geo_optim.SD_scaling
    UtotHis = geo_optim.UtotHis
    verbosity = geo_optim.verbosity


    Max_Force = 0.0
    for atom = 1:Natom, xyz = 1:3
        if atom_Fixed_XYZ[atom][xyz] == 0
            Max_Force = max(Max_Force, abs(ForceAll[atom,xyz]))
        end
    end
    Geo_Opt_Max_Force[GeoOpt_iter] = Max_Force



    diis_iter = ifelse(GeoOpt_iter < M_GDIIS_HISTORY, GeoOpt_iter, M_GDIIS_HISTORY)

    for i = diis_iter-1:-1:1
        for atom = 1:Natom, xyz = 1:3
            GxyzHisIn[i+1][atom][xyz] = GxyzHisIn[i][atom][xyz]
            GxyzHisR[i+1][atom][xyz] = GxyzHisR[i][atom][xyz]
        end
    end

    for atom = 1:Natom, xyz = 1:3
        GxyzHisIn[1][atom][xyz] = Gxyz[atom][xyz]
        if atom_Fixed_XYZ[atom][xyz] == 0
            GxyzHisR[1][atom][xyz] = ForceAll[atom,xyz]
        else
            GxyzHisR[1][atom][xyz] = 0.0
        end
    end


    if GeoOpt_iter == 1
        scaling_factor = 2.0
        Estimate_Initial_Hessian(geo_optim)
    end

    Geo_Hessian = geo_optim.Geo_Hessian



    A = zeros(Float64, diis_iter+1, diis_iter+1)
    B = zeros(Float64, diis_iter+1)
    ko = zeros(Float64, 3*Natom+2)
    U = zeros(Float64, 3*Natom+2, 3*Natom+2)
    B[end] = 1.0


    sMD_TimeStep = 0.05/(0.01*41.341105)

    for i = 1:diis_iter, j = 1:diis_iter
        RR = 0.0
        for atom = 1:Natom, xyz = 1:3
            RR += GxyzHisR[i][atom][xyz]*GxyzHisR[j][atom][xyz]
        end

        A[i,j] = RR
    end

    MaxA = 0.0
    for i = 1:diis_iter, j = 1:diis_iter
        MaxA = max(MaxA, abs(A[i,j]))
    end

    for i = 1:diis_iter, j = 1:diis_iter
        A[i,j] = A[i,j]/MaxA
    end
    for j = 1:diis_iter
        A[j,end] = 1.0
        A[end,j] = 1.0
    end
    A[end,end] = 0.0


    if myrank == 0 && verbosity >= 1
        @printf("<EF>  DIIS matrix\n")
        for i = 1:diis_iter+1
            @printf("<EF> ")
            for j = 1:diis_iter+1
                @printf("%6.3f ", A[i,j])
            end
            @printf("\n")
        end

        @printf("M_GDIIS_HISTORY = %2d  diis_iter = %2d\n", M_GDIIS_HISTORY, diis_iter)
    end
    MPI.Barrier(comm)

    sA = Symmetric(A, :L)
    BB = sA\B


    if myrank == 0 && verbosity >= 1
        @printf("<EF> diis alpha = ")
        for i = 1:diis_iter+1
            @printf("%f ", BB[i])
        end
        @printf("\n")
    end


    if Max_Force < Geo_Opt_criterion
        Geo_Opt_convergence = true

        myrank == 0 && println("\nThe geometry optimization was achieved.\n")
    end


    Max_Step = 0.0

    if !Geo_Opt_convergence && GeoOpt_iter ≠ Geo_Opt_Max

        for atom = 1:Natom, xyz = 1:3
            Gxyz[atom][xyz] = 0.0
        end

        for atom = 1:Natom, xyz = 1:3
            if atom_Fixed_XYZ[atom][xyz] == 0
                for i = 1:diis_iter
                    Gxyz[atom][xyz] += GxyzHisR[i][atom][xyz]*BB[i]
                end
            end
        end

        m1 = 1
        for atom = 1:Natom, xyz = 1:3
            m1 += 1
            Geo_Hessian[m1,1] = Gxyz[atom][xyz]
        end


        SumB = 0.0
        for atom = 1:Natom, xyz = 1:3
            SumB += Gxyz[atom][xyz]*Gxyz[atom][xyz]
        end
        SumB = sqrt(SumB)/Natom

        if myrank == 0 && verbosity >= 1
            @printf("<EF> |tilde{R}|=%E\n", SumB)
        end



        if GeoOpt_iter ≠ 1
            for i = 1:3*Natom
                Sum1 = 0.0
                for atom = 1:Natom
                    Sum1 += Geo_Hessian[i+1,3*atom-1]*(GxyzHisIn[1][atom][1]-GxyzHisIn[2][atom][1])
                    Sum1 += Geo_Hessian[i+1,3*atom]*(GxyzHisIn[1][atom][2]-GxyzHisIn[2][atom][2])
                    Sum1 += Geo_Hessian[i+1,3*atom+1]*(GxyzHisIn[1][atom][3]-GxyzHisIn[2][atom][3])
                end
                Geo_Hessian[1,i+1] = Sum1
            end


            tmp1 = 0.0
            tmp2 = 0.0
            for atom = 1:Natom
                tmp1 += (GxyzHisIn[1][atom][1]-GxyzHisIn[2][atom][1])*(GxyzHisR[1][atom][1]-GxyzHisR[2][atom][1])
                tmp1 += (GxyzHisIn[1][atom][2]-GxyzHisIn[2][atom][2])*(GxyzHisR[1][atom][2]-GxyzHisR[2][atom][2])
                tmp1 += (GxyzHisIn[1][atom][3]-GxyzHisIn[2][atom][3])*(GxyzHisR[1][atom][3]-GxyzHisR[2][atom][3])

                tmp2 += Geo_Hessian[1,3*atom-1]*(GxyzHisIn[1][atom][1]-GxyzHisIn[2][atom][1])
                tmp2 += Geo_Hessian[1,3*atom]*(GxyzHisIn[1][atom][2]-GxyzHisIn[2][atom][2])
                tmp2 += Geo_Hessian[1,3*atom+1]*(GxyzHisIn[1][atom][3]-GxyzHisIn[2][atom][3])
            end

            c0 = 1/tmp1
            c1 = 1/tmp2

            if c0 > 0.0
                m1 = 1
                for atom1 = 1:Natom, xyz1 = 1:3
                    m1 += 1
                    m2 = 1
                    for atom2 = 1:Natom, xyz2 = 1:3
                        m2 += 1
                        Geo_Hessian[m1,m2] += c0*(GxyzHisR[1][atom1][xyz1]-GxyzHisR[2][atom1][xyz1])*(GxyzHisR[1][atom2][xyz2]-GxyzHisR[2][atom2][xyz2]) - c1*Geo_Hessian[1,m1]*Geo_Hessian[1,m2]
                    end
                end
            end
        end


        ko[2:3*Natom+1], U[2:3*Natom+1,2:3*Natom+1] = eigen(Symmetric(Geo_Hessian[2:3*Natom+1,2:3*Natom+1]))


        MinKo = ifelse(Natom<=4, 0.10, 0.005)

        for i = 1:3*Natom
            if ko[i+1] < MinKo
                ko[i+1] = MinKo
            end
        end

        for i = 1:3*Natom
            Sum = 0.0
            for j = 1:3*Natom
                Sum += U[j+1,i+1]*Geo_Hessian[j+1,1]
            end
            U[i+1,1] = Sum
        end

        for i = 1:3*Natom
            Sum = 0.0
            for j = 1:3*Natom
                Sum += U[i+1,j+1]*U[j+1,1]/ko[j+1]
            end
            U[1,i+1] = Sum
        end


        if Utot0<Etot && GeoOpt_iter ≠ 1
            scaling_factor = 0.95*scaling_factor
        end


        for atom = 1:Natom, xyz = 1:3
            Gxyz[atom][xyz] = 0.0
        end

        m1 = 1
        for atom = 1:Natom, xyz = 1:3
            for i = 1:diis_iter
                Gxyz[atom][xyz] += GxyzHisIn[i][atom][xyz]*BB[i]
            end

            m1 += 1
            Gxyz[atom][xyz] -= scaling_factor*U[1,m1]
        end

        


        Max_Step = 0.0
        for atom = 1:Natom, xyz = 1:3
            diff = abs(Gxyz[atom][xyz] - GxyzHisIn[1][atom][xyz])
            if Max_Step < diff
                Max_Step = diff
            end
        end

        Criterion_Max_Step = 0.5
        if Max_Step > Criterion_Max_Step
            for atom = 1:Natom, xyz = 1:3
                diff = Gxyz[atom][xyz] - GxyzHisIn[1][atom][xyz]
                Gxyz[atom][xyz] = GxyzHisIn[1][atom][xyz] + diff/Max_Step*Criterion_Max_Step
            end
        end


        Max_Step = 0.0
        for atom = 1:Natom, xyz = 1:3
            diff = abs(Gxyz[atom][xyz] - GxyzHisIn[1][atom][xyz])
            if Max_Step < diff
                Max_Step = diff
            end
        end
    end


    SD_scaling_user = SD_scaling*Max_Force*0.2
    Utot0 = Etot

    for atom = 1:Natom, xyz = 1:3
        if atom_Fixed_XYZ[atom][xyz] == 1
            Gxyz[atom][xyz] = Gxyz_fix[atom][xyz]
        end
    end


    if myrank == 0 && verbosity >= 1
        @printf("<EF>  diff_x = %5.12f , |dE| = %5.12f\n", Max_Step, abs(Etot-UtotHis[1]))
        @printf("<EF>  |Maximum force| (Hartree/Bohr) = %5.12f\n", Max_Force)
        @printf("<EF>  Criterion       (Hartree/Bohr) = %5.12f\n", Geo_Opt_criterion)

        @printf("\n")
        @printf("<ForceAll>\n")
        for atom = 1:Natom
            @printf("  Fxyz(a.u.)  atom = %3d  %15.12f  %15.12f  %15.12f\n", atom, ForceAll[atom,1], ForceAll[atom,2], ForceAll[atom,3])
        end

        @printf("<Atomic positions>  New atomic positions using DIIS_EF\n")
        for atom = 1:Natom
            @printf("  XYZ(ang)    atom = %3d  %15.12f  %15.12f  %15.12f\n",  atom, Gxyz[atom][1]/Ang_to_bohr, Gxyz[atom][2]/Ang_to_bohr, Gxyz[atom][3]/Ang_to_bohr)
        end
    end
    MPI.Barrier(comm)


    geo_optim.Gxyz = Gxyz
    geo_optim.GxyzHisIn = GxyzHisIn
    geo_optim.GxyzHisR = GxyzHisR
    geo_optim.Geo_Opt_Max_Force = Geo_Opt_Max_Force
    geo_optim.Geo_Opt_convergence = Geo_Opt_convergence
    geo_optim.Geo_Hessian = Geo_Hessian
    geo_optim.SD_scaling_user = SD_scaling_user
    geo_optim.SD_scaling = SD_scaling
    geo_optim.Utot0 = Utot0
    geo_optim.scaling_factor = scaling_factor
    geo_optim.UtotHis = UtotHis
end


function Estimate_Initial_Hessian(geo_optim::Geo_Optim)
    
    Natom = geo_optim.Natom
    Gxyz = geo_optim.Gxyz
    FNAN = geo_optim.FNAN
    natn = geo_optim.natn
    ncn = geo_optim.ncn
    Dis = geo_optim.Dis
    atv = geo_optim.atv
    Atoms_Core_Charge = geo_optim.Atoms_Core_Charge
    Geo_Hessian = geo_optim.Geo_Hessian


    B = zeros(Float64,8,8)
    B[1,1] =  0.5000; B[1,2] =  0.5000; B[1,3] =  0.5000; B[1,4] =  0.5000; B[1,5] =  0.5000; B[1,6] =  0.5000; B[1,7] =  0.5000; B[1,8] =  0.5000;
    B[2,1] =  0.5000; B[2,2] = -0.2573; B[2,3] =  0.3401; B[2,4] =  0.6937; B[2,5] =  0.7126; B[2,6] =  0.8335; B[2,7] =  0.9491; B[2,8] =  1.0000;
    B[3,1] =  0.5000; B[3,2] =  0.3401; B[3,3] =  0.9652; B[3,4] =  1.2843; B[3,5] =  1.4625; B[3,6] =  1.6549; B[3,7] =  1.7190; B[3,8] =  2.0000;
    B[4,1] =  0.5000; B[4,2] =  0.6937; B[4,3] =  1.2843; B[4,4] =  1.6925; B[4,5] =  1.8238; B[4,6] =  2.1164; B[4,7] =  2.3185; B[4,8] =  2.5000;
    B[5,1] =  0.5000; B[5,2] =  0.7126; B[5,3] =  1.4625; B[5,4] =  1.8238; B[5,5] =  2.0203; B[5,6] =  2.2137; B[5,7] =  2.5206; B[5,8] =  2.7000;
    B[6,1] =  0.5000; B[6,2] =  0.8335; B[6,3] =  1.6549; B[6,4] =  2.1164; B[6,5] =  2.2137; B[6,6] =  2.3718; B[6,7] =  2.5110; B[6,8] =  2.7000;
    B[7,1] =  0.5000; B[7,2] =  0.9491; B[7,3] =  1.7190; B[7,4] =  2.3185; B[7,5] =  2.5206; B[7,6] =  2.5110; B[7,7] =  2.5200; B[7,8] =  2.7000;
    B[8,1] =  0.5000; B[8,2] =  1.0000; B[8,3] =  2.0000; B[8,4] =  2.5000; B[8,5] =  2.7000; B[8,6] =  2.7000; B[8,7] =  2.7000; B[8,8] =  2.9000;


    for i = 1:3*Natom
        for j = 1:3*Natom
            Geo_Hessian[i+1,j+1] = 0.0
        end
        Geo_Hessian[i+1,i+1] = 0.1
    end


    g = zeros(Float64, 3)

    for atom = 1:Natom

        Zcore = Atoms_Core_Charge[atom]

        if Zcore == 1
            n1 = 1
        elseif Zcore <= 10
            n1 = 2
        elseif Zcore <= 18
            n1 = 3
        elseif Zcore <= 36
            n1 = 4
        elseif Zcore <= 54
            n1 = 5
        elseif Zcore <= 86
            n1 = 6
        elseif Zcore <= 103
            n1 = 7
        end

        if Zcore==2 || Zcore==10 || Zcore==36 || Zcore==54 || Zcore==86
            n1 = 0
        end

        for Rn = 2:FNAN[atom]+1

            jatom = natn[atom][Rn]
            cell = ncn[atom][Rn]
            R = Dis[atom][Rn]

            Zcore = Atoms_Core_Charge[jatom]

            if Zcore == 1
                n2 = 1
            elseif Zcore <= 10
                n2 = 2
            elseif Zcore <= 18
                n2 = 3
            elseif Zcore <= 36
                n2 = 4
            elseif Zcore <= 54
                n2 = 5
            elseif Zcore <= 86
                n2 = 6
            elseif Zcore <= 103
                n2 = 7
            end
    
            if Zcore==2 || Zcore==10 || Zcore==36 || Zcore==54 || Zcore==86
                n2 = 0
            end

            g[1] = (Gxyz[atom][1] - Gxyz[jatom][1] - atv[cell+1][1])/R
            g[2] = (Gxyz[atom][2] - Gxyz[jatom][2] - atv[cell+1][2])/R
            g[3] = (Gxyz[atom][3] - Gxyz[jatom][3] - atv[cell+1][3])/R

            d = R - B[n1+1,n2+1]
            gr = 1.734/d^3

            for i = 1:3, j = 1:3
                Geo_Hessian[3*(atom-1)+i+1,3*(atom-1)+j+1] += g[i]*g[j]*gr
            end

            for i = 1:3, j = 1:3
                Geo_Hessian[3*(atom-1)+i+1,3*(jatom-1)+j+1] += -g[i]*g[j]*gr
            end
        end
    end


    geo_optim.Geo_Hessian = Geo_Hessian
end