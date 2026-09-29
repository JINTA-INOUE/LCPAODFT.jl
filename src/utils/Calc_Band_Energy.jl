function Calc_Band_Energy!(SpinPol::String, Hks, DM, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_Hoffset = system_grid.MPI_Hoffset
    Total_NumOrbs = system_grid.Total_NumOrbs
    band_energy = 0.0
    spinsize = ifelse(SpinPol=="off", 1, 2)
    Spindeg = ifelse(SpinPol=="off", 2, 1)

    for spin = 1:spinsize
        local_energy = 0.0
        for pair in eachindex(MPI_atom)
            atom = MPI_atom[pair]
            jatom = MPI_natn[pair]
            offset = MPI_Hoffset[pair]
            block_size = Total_NumOrbs[atom]*Total_NumOrbs[jatom]
            for index = offset+1:offset+block_size
                local_energy += DM[spin][index]*Hks[spin][index]
            end
        end
        band_energy += MPI.Allreduce(local_energy, MPI.SUM, comm)
    end
    Eele = MPI.Bcast(Spindeg*band_energy, 0, comm)

    return Eele
end


function Calc_Band_Energy!(Hks, iHks, DM, iDM, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_Hoffset = system_grid.MPI_Hoffset
    Total_NumOrbs = system_grid.Total_NumOrbs

    local_energy = 0.0
    for pair in eachindex(MPI_atom)
        atom = MPI_atom[pair]
        Rn = MPI_FNAN[pair]
        jatom = MPI_natn[pair]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        offset = MPI_Hoffset[pair]

        iHks1 = iHks[1][atom][Rn]
        iHks2 = iHks[2][atom][Rn]
        iHks3 = iHks[3][atom][Rn]
        index = offset
        for ist = 1:NO0, jst = 1:NO1
            index += 1
            local_energy +=    DM[1][index]*Hks[1][index]
            local_energy -=   iDM[1][index]*iHks1[ist][jst]
            local_energy +=    DM[2][index]*Hks[2][index]
            local_energy -=   iDM[2][index]*iHks2[ist][jst]
            local_energy +=  2*DM[3][index]*Hks[3][index]
            local_energy -=  2*DM[4][index]*(Hks[4][index] + iHks3[ist][jst])
        end
    end
    band_energy = MPI.Allreduce(local_energy, MPI.SUM, comm)
    Eele = MPI.Bcast(band_energy, 0, comm)

    return Eele
end


#=
function Calc_Band_Energy!(electron::ClusterBloch)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Spindeg = electron.Spindeg
    TotalZ = electron.TotalZ
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    spinsize = electron.spinsize
    Nfsize = electron.Nfsize
    ChemP = electron.ChemP
    Enk = electron.Enk
    FF = electron.FF
    ChemP = MPI.Bcast(ChemP, 0, comm)
    electron.ChemP = ChemP


    Eele0 = 0.0
    @inbounds for spin = 1:spinsize, μ = 1:Nfsize
            
        x = (Enk[μ,spin]-ChemP)*Beta
        if x <= -max_x
            x = -max_x
         end

        if x >= max_x
            x = max_x
        end
            
        FermiF = 1/(1 + exp(x))
        FF[μ,spin] = FermiF
        Eele0 += FermiF*Enk[μ,spin]
    end
    Eele0 = Spindeg*Eele0
    Eele0 = MPI.Bcast(Eele0, 0, comm)

    electron.FF = FF
    electron.Eele = Eele0
end


function Calc_Band_Energy!(electron::CrystalBloch, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Spindeg = electron.Spindeg
    TotalZ = electron.TotalZ
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    spinsize = electron.spinsize
    Nfsize = electron.Nfsize
    Nkpt = electron.Nkpt
    ChemP = electron.ChemP
    Enk = electron.Enk
    FF = electron.FF
    AllNkpt = kpoints.AllNkpt
    All_kweight = kpoints.All_kweight

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
=#