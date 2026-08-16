@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Col!(cal_force::Bool, OLP, Hks, DM, EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPI_kweight = kpoints.MPI_kweight
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    spinsize = electron.spinsize
    Enk = electron.Enk
    FF = electron.FF

    for spin = 1:spinsize
        fill!(DM[spin], 0.0)
        fill!(EDM[spin], 0.0)
    end
    fill!(Enk, 0.0)

    Swork = zeros(ComplexF64, fsize, fsize)
    Hwork = zeros(ComplexF64, fsize, fsize)
    if cal_force
        Cwork = zeros(ComplexF64, fsize, fsize)
    end

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        HS_matrix!(Swork, Hwork, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[:,knum+ik,spin] .= eigvals!(Hermitian(Hwork), Hermitian(Swork))
    end

    MPI.Allreduce!(Enk, MPI.SUM, comm)
    Calc_Band_Energy!(electron, kpoints)

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        HS_matrix!(Swork, Hwork, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Cnk = eigen!(Hermitian(Hwork), Hermitian(Swork)).vectors
        if cal_force
            copyto!(Cwork, Cnk)
            Cnk = Cwork
        end
        @inbounds for μ = 1:fsize
            weights = sqrt(FF[μ,knum+ik,spin]*MPI_kweight[ik])
            @views rmul!(Cnk[:,μ], weights)
        end

        if cal_force
            copyto!(Swork, Cnk)
            @inbounds for μ = 1:fsize
                @views rmul!(Swork[:,μ], Enk[μ,knum+ik,spin])
            end
            mul!(Hwork, Swork, adjoint(Cnk))
            mul!(Swork, Cnk, adjoint(Cnk))
            _accumulate_collinear!(DM[spin], EDM[spin], Swork, Hwork, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik], AllNkpt)
        else
            mul!(Swork, Cnk, adjoint(Cnk))
            _accumulate_collinear!(DM[spin], Swork, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik], AllNkpt)
        end
    end

    for spin = 1:spinsize
        MPI.Allreduce!(DM[spin], MPI.SUM, comm)
        if cal_force
            MPI.Allreduce!(EDM[spin], MPI.SUM, comm)
        end
    end
end


@timeit timer "Crystal_DFT_NonCollinear" function Crystal_DFT_NonCol!(cal_force, OLP, Hks, iHks, DM, iDM, EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPI_kweight = kpoints.MPI_kweight
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    FF = electron.FF

    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)
    fill!(DM[3], 0.0)
    fill!(DM[4], 0.0)
    fill!(EDM[1], 0.0)
    fill!(EDM[2], 0.0)
    fill!(Enk, 0.0)

    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
       fill!(iDM[1][atom][Rn][ist], 0.0) 
       fill!(iDM[2][atom][Rn][ist], 0.0) 
    end

    tmpH = zeros(ComplexF64, fsize, fsize)
    Swork = zeros(ComplexF64, Nfsize, Nfsize)
    Hwork = zeros(ComplexF64, Nfsize, Nfsize)
    if cal_force
        Cwork = zeros(ComplexF64, Nfsize, Nfsize)
    end

    for ik = 1:MPI_Nkpt
        HS_matrix_NC!(Hwork, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        fill!(Swork, 0.0)
        @. @views Swork[1:fsize, 1:fsize] = tmpH
        @. @views Swork[fsize+1:end, fsize+1:end] = tmpH
        Enk[:,knum+ik,1] .= eigvals!(Hermitian(Hwork), Hermitian(Swork))
    end

    MPI.Allreduce!(Enk, MPI.SUM, comm)
    Calc_Band_Energy!(electron, kpoints)

    for ik = 1:MPI_Nkpt
        HS_matrix_NC!(Hwork, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        fill!(Swork, 0.0)
        @. @views Swork[1:fsize, 1:fsize] = tmpH
        @. @views Swork[fsize+1:end, fsize+1:end] = tmpH
        Cnk = eigen!(Hermitian(Hwork), Hermitian(Swork)).vectors
        if cal_force
            copyto!(Cwork, Cnk)
            Cnk = Cwork
        end
        @inbounds for μ = 1:Nfsize
            weights = sqrt(FF[μ,knum+ik,1]*MPI_kweight[ik])
            @views rmul!(Cnk[:,μ], weights)
        end
        if cal_force
            copyto!(Swork, Cnk)
            @inbounds for μ = 1:Nfsize
                @views rmul!(Swork[:,μ], Enk[μ,knum+ik,1])
            end
            mul!(Hwork, Swork, adjoint(Cnk))
            mul!(Swork, Cnk, adjoint(Cnk))
            _accumulate_noncollinear!(DM, iDM, EDM, Swork, Hwork, Natom, fsize, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik], AllNkpt)
        else
            mul!(Swork, Cnk, adjoint(Cnk))
            _accumulate_noncollinear!(DM, iDM, Swork, Natom, fsize, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik], AllNkpt)
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
    MPI.Allreduce!(DM[3], MPI.SUM, comm)
    MPI.Allreduce!(DM[4], MPI.SUM, comm)
    if cal_force
        MPI.Allreduce!(EDM[1], MPI.SUM, comm)
        MPI.Allreduce!(EDM[2], MPI.SUM, comm)
    end
end


@timeit timer "_accumulate_collinear!" function _accumulate_collinear!(DM, Swork, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts, AllNkpt)
    ka, kb, kc = kpts
    DMst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        Anum = MP[atom]
        NO0 = Total_NumOrbs[atom]
        cell_index = ncn[atom][Rn]+1
        jatom = natn[atom][Rn]
        Bnum = MP[jatom]
        NO1 = Total_NumOrbs[jatom]
        cell = atv_ijk[cell_index]
        kRn = ka*cell[1] + kb*cell[2] + kc*cell[3]
        ex = cispi(2*kRn)/AllNkpt
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DMelement = Swork[Bnum+jst,Anum+ist]
            DM[DMst] += real(DMelement*ex)
        end
    end
end


@timeit timer "_accumulate_collinear!" function _accumulate_collinear!(DM, EDM, Swork, Hwork, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts, AllNkpt)
    ka, kb, kc = kpts
    DMst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        Anum = MP[atom]
        NO0 = Total_NumOrbs[atom]
        cell_index = ncn[atom][Rn]+1
        jatom = natn[atom][Rn]
        Bnum = MP[jatom]
        NO1 = Total_NumOrbs[jatom]
        cell = atv_ijk[cell_index]
        kRn = ka*cell[1] + kb*cell[2] + kc*cell[3]
        ex = cispi(2*kRn)/AllNkpt
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DMelement = Swork[Bnum+jst,Anum+ist]
            EDMelement = Hwork[Bnum+jst,Anum+ist]
            DM[DMst] += real(DMelement*ex)
            EDM[DMst] += real(EDMelement*ex)
        end
    end
end


@timeit timer "_accumulate_noncollinear!" function _accumulate_noncollinear!(DM, iDM, Swork, Natom, fsize, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts, AllNkpt)
    ka, kb, kc = kpts
    DMst = 0
    DM1 = DM[1]
    DM2 = DM[2]
    DM3 = DM[3]
    DM4 = DM[4]
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        Anum = MP[atom]
        NO0 = Total_NumOrbs[atom]
        cell_index = ncn[atom][Rn]+1
        jatom = natn[atom][Rn]
        Bnum = MP[jatom]
        NO1 = Total_NumOrbs[jatom]
        cell = atv_ijk[cell_index]
        kRn = ka*cell[1] + kb*cell[2] + kc*cell[3]
        ex = cispi(2*kRn)/AllNkpt
        iDM1 = iDM[1][atom][Rn]
        iDM2 = iDM[2][atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DMelement_uu = Swork[Bnum+jst,Anum+ist]
            DMelement_dd = Swork[fsize+Bnum+jst,fsize+Anum+ist]
            DMelement_ud = Swork[fsize+Bnum+jst,Anum+ist]
            DM1[DMst] += real(DMelement_uu*ex)
            DM2[DMst] += real(DMelement_dd*ex)
            DM3[DMst] += real(DMelement_ud*ex)
            DM4[DMst] += imag(DMelement_ud*ex)
            iDM1[ist][jst] += imag(DMelement_uu*ex)
            iDM2[ist][jst] += imag(DMelement_dd*ex)
        end
    end
end


@timeit timer "_accumulate_noncollinear!" function _accumulate_noncollinear!(DM, iDM, EDM, Swork, Hwork, Natom, fsize, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts, AllNkpt)
    ka, kb, kc = kpts
    DMst = 0
    DM1 = DM[1]
    DM2 = DM[2]
    DM3 = DM[3]
    DM4 = DM[4]
    EDM1 = EDM[1]
    EDM2 = EDM[2]
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        Anum = MP[atom]
        NO0 = Total_NumOrbs[atom]
        cell_index = ncn[atom][Rn]+1
        jatom = natn[atom][Rn]
        Bnum = MP[jatom]
        NO1 = Total_NumOrbs[jatom]
        cell = atv_ijk[cell_index]
        kRn = ka*cell[1] + kb*cell[2] + kc*cell[3]
        ex = cispi(2*kRn)/AllNkpt
        iDM1 = iDM[1][atom][Rn]
        iDM2 = iDM[2][atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DMelement_uu = Swork[Bnum+jst,Anum+ist]
            DMelement_dd = Swork[fsize+Bnum+jst,fsize+Anum+ist]
            DMelement_ud = Swork[fsize+Bnum+jst,Anum+ist]
            EDMelement_uu = Hwork[Bnum+jst,Anum+ist]
            EDMelement_dd = Hwork[fsize+Bnum+jst,fsize+Anum+ist]
            DM1[DMst] += real(DMelement_uu*ex)
            DM2[DMst] += real(DMelement_dd*ex)
            DM3[DMst] += real(DMelement_ud*ex)
            DM4[DMst] += imag(DMelement_ud*ex)
            iDM1[ist][jst] += imag(DMelement_uu*ex)
            iDM2[ist][jst] += imag(DMelement_dd*ex)
            EDM1[DMst] += real(EDMelement_uu*ex)
            EDM2[DMst] += real(EDMelement_dd*ex)
        end
    end
end