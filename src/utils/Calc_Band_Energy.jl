function Calc_Band_Energy!(electron::CrystalBloch, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Spindeg = electron.Spindeg
    TotalZ = electron.TotalZ
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    spinsize = electron.spinsize
    Nfsize = electron.Nfsize
    Nkpt = electron.Nkpt
    Enk = electron.Enk
    FF = electron.FF
    AllNkpt = kpoints.AllNkpt
    All_kweight = kpoints.All_kweight

    ChemP = Calc_ChemP(spinsize, Spindeg, Nfsize, Nkpt, TotalZ, All_kweight, Beta, Enk)
    ChemP = MPI.Bcast(ChemP, 0, comm)
    electron.ChemP = ChemP


    Eele0 = 0.0
    @inbounds for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:Nfsize
            
        x = (Enk[μ,ik,spin]-ChemP)*Beta
        if x <= -max_x
            x = -max_x
         end

        if x >= max_x
            x = max_x
        end
            
        FermiF = 1/(1 + exp(x))
        FF[μ,ik,spin] = FermiF
        Eele0 += FermiF*Enk[μ,ik,spin]*All_kweight[ik]
    end
    Eele0 = Spindeg*Eele0/AllNkpt
    Eele0 = MPI.Bcast(Eele0, 0, comm)

    electron.FF = FF
    electron.Eele = Eele0
end