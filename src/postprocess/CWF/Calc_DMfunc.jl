@timeit timer "Calc_DMfunc" function Calc_DMfunc(spinsize, Ngsize, Σmk, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Nkpt = kpoints.Nkpt
    MPI_Nkpt = kpoints.MPI_Nkpt


    DM = 0.0
    for spin = 1:spinsize, ik = 1:MPI_Nkpt, p = 1:Ngsize
        if Σmk[spin][ik][p] < 0.0
            error("singler values is negative")
        end
        DM += (Σmk[spin][ik][p]-1)^2
    end
    DM = MPI.Allreduce(DM, MPI.SUM, comm)
    DM = DM/Nkpt/Ngsize/spinsize

    myrank == 0 && println("\tDM per CWF is $(DM)")


    return DM
end