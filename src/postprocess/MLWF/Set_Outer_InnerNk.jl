function Set_Outer_InnerNk(Enk, kpoints::KPoints, mlwf_setup::MLWF_Setup)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    ChemP = material.ChemP

    WANNUM = mlwf_setup.WANNUM
    MLWF_Outer_Window_Bottom = mlwf_setup.MLWF_Outer_Window_Bottom
    MLWF_Outer_Window_Top = mlwf_setup.MLWF_Outer_Window_Top
    MLWF_Inner_Window_Bottom = mlwf_setup.MLWF_Inner_Window_Bottom
    MLWF_Inner_Window_Top = mlwf_setup.MLWF_Inner_Window_Top
    MLWF_kmesh = mlwf_setup.MLWF_kmesh
    Nkpt = prod(MLWF_kmesh)
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    MLWF_Outer_Window_Bottom = MLWF_Outer_Window_Bottom/eV2Hartree + ChemP
    MLWF_Outer_Window_Top = MLWF_Outer_Window_Top/eV2Hartree + ChemP
    MLWF_Inner_Window_Bottom = MLWF_Inner_Window_Bottom/eV2Hartree + ChemP
    MLWF_Inner_Window_Top = MLWF_Inner_Window_Top/eV2Hartree + ChemP

    BANDNUM = 0
    MkNUM = 0
    Nk = zeros(Int32, spinsize, Nkpt, 2)
    innerNk = zeros(Int32, spinsize, Nkpt, 2)

    for ik = 1:MPI_Nkpt, spin = 1:spinsize

        Nk[spin,ik+knum,1] = 0
        Nk[spin,ik+knum,2] = 0
        innerNk[spin,ik+knum,1] = 0
        innerNk[spin,ik+knum,2] = 0

        for μ = 1:Nfsize
            if Enk[spin][ik][μ] < MLWF_Outer_Window_Bottom
                Nk[spin,ik+knum,1] = μ
            end

            if Enk[spin][ik][μ] < MLWF_Outer_Window_Top
                Nk[spin,ik+knum,2] = μ
            end

            if Enk[spin][ik][μ] < MLWF_Inner_Window_Bottom
                innerNk[spin,ik+knum,1] = μ
            end

            if Enk[spin][ik][μ] < MLWF_Inner_Window_Top
                innerNk[spin,ik+knum,2] = μ
            end
        end


        if BANDNUM < (Nk[spin,ik+knum,2]-Nk[spin,ik+knum,1])
            BANDNUM = Nk[spin,ik+knum,2]-Nk[spin,ik+knum,1]
        end

        if MkNUM < (innerNk[spin,ik+knum,2]-innerNk[spin,ik+knum,1])
            MkNUM = innerNk[spin,ik+knum,2] - innerNk[spin,ik+knum,1]
        end

        if WANNUM > Nk[spin,ik+knum,2]-Nk[spin,ik+knum,1]
            if myrank == 0
                @show spin, ik
                @show Nk[spin,ik+knum,2], Nk[spin,ik+knum,1]
                println("**************************ERROR**************************")
                println("* Bands number within OUTER window [$MLWF_Outer_Window_Bottom, $MLWF_Outer_Window_Top]")
                println("* is less than MLWF Function number $WANNUM.")
                println("* Please check them and try again.")
                println("**************************ERROR**************************")
            end
	        MPI.Finalized()
	        error("")
        end

        if WANNUM < innerNk[spin,ik+knum,2]-innerNk[spin,ik+knum,1]
            if myrank == 0
                @show spin, ik
                @show innerNk[spin,ik+knum,2], innerNk[spin,ik+knum,1]
                println("**************************ERROR**************************")
                println("* Bands number within INNER window [$MLWF_Inner_Window_Bottom, $MLWF_Inner_Window_Top]")
                println("* is larger than MLWF Function number $WANNUM.")
                println("* Please check them and try again.")
                println("**************************ERROR**************************")
            end
            MPI.Finalized()
	        error("")
        end
    end

    MPI.Allreduce!(Nk, MPI.SUM, comm)
    MPI.Allreduce!(innerNk, MPI.SUM, comm)
    BANDNUM = MPI.Allreduce(BANDNUM, max, comm)

    if myrank == 0
        println("Selected bands within the Outer and Inner Windows at each k point:")
        println(" spin| kpt | Outer Window (Nk) |  Inner Window (Mk)  |")
        for ik = 1:Nkpt, spin = 1:spinsize
            @printf("  %1d  |%5d|", spin-1, ik)
            @printf("  (%3d ,%3d]  %3d  |  (%3d ,%3d]   %3d   |\n", 
                    Nk[spin,ik,1], Nk[spin,ik,2], Nk[spin,ik,2]-Nk[spin,ik,1], 
                    innerNk[spin,ik,1], innerNk[spin,ik,2], innerNk[spin,ik,2] - innerNk[spin,ik,1])
        end
    end
    MPI.Barrier(comm)


    return BANDNUM, Nk
end