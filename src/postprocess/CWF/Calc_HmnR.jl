function Calc_HmnR!(HmnR, spinsize, Nfsize, NCell, Ngsize, cell_list_ijk, Umnk, Enk, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    Nkpt = kpoints.Nkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    for spin = 1:spinsize, cell = 1:NCell, pst = 1:Ngsize, qst = 1:Ngsize
                
        Sum = ComplexF64(0.0, 0.0)
        for ik = 1:MPI_Nkpt
            tmp = ComplexF64(0.0, 0.0)
            for μ = 1:Nfsize
                tmp += Enk[spin,μ,ik] * conj(Umnk[spin][ik][μ,pst]) * Umnk[spin][ik][μ,qst]
            end

            kRn = dot(MPI_kpts[ik], cell_list_ijk[cell])
            ex = cispi(-2*kRn)
            Sum += tmp*ex
        end

        HmnR[qst,pst,cell,spin] = Sum/Nkpt
    end

    MPI.Allreduce!(HmnR, MPI.SUM, comm)
end