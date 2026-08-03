@timeit timer "Set_Density_Grid" function Set_Density_Grid_nonpol!(ucell::UCell, Orbs_Grid, DM, Density_Grid)
   
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]

    DMsum = 0
    DM_atom = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs))
    ai_tempDGs = zeros(Float64, maximum(GridN_Atom))

    fill!(Density_Grid[1], 0.0)

    for loop = 1:MPI_size

        fill!(ai_tempDGs, 0.0)

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        
        DMst = 0
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DM_atom[jst,ist] = DM[1][DMnum+DMsum+DMst]
        end
        DMsum += DMst

        _Calc_Den2_nonpol!(ai_tempDGs, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], DM_atom)
        
        @inbounds for xyz = 1:GridN_Atom[atom]
            r = GridListAtom[atom][xyz]+1
            Density_Grid[1][r] += ai_tempDGs[xyz]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
end


@timeit timer "Set_Density_Grid" function Set_Density_Grid_pol!(ucell::UCell, Orbs_Grid, DM, Density_Grid)
   
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]

    DMsum = 0
    DM_atom = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 2)
    ai_tempDGs = zeros(Float64, maximum(GridN_Atom), 2)

    fill!(Density_Grid[1], 0.0)
    fill!(Density_Grid[2], 0.0)

    for loop = 1:MPI_size

        fill!(ai_tempDGs, 0.0)

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        DMst = 0
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DM_atom[jst,ist,1] = DM[1][DMnum+DMsum+DMst]
            DM_atom[jst,ist,2] = DM[2][DMnum+DMsum+DMst]
        end
        DMsum += DMst

        _Calc_Den2_pol!(ai_tempDGs, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], DM_atom)
        
        @inbounds for xyz = 1:GridN_Atom[atom]
            r = GridListAtom[atom][xyz]+1
            Density_Grid[1][r] += ai_tempDGs[xyz,1]
            Density_Grid[2][r] += ai_tempDGs[xyz,2]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[2], MPI.SUM, comm)
end


@timeit timer "Set_Density_Grid" function Set_Density_Grid_nc!(ucell::UCell, Orbs_Grid, DM, Density_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]

    DMsum = 0
    DM_atom = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 4)
    ai_tempDGs = zeros(Float64, maximum(GridN_Atom), 4)

    fill!(Density_Grid[1], 0.0)
    fill!(Density_Grid[2], 0.0)
    fill!(Density_Grid[3], 0.0)
    fill!(Density_Grid[4], 0.0)

    for loop = 1:MPI_size

        fill!(ai_tempDGs, 0.0)

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        
        DMst = 0
        @inbounds for ist = 1:NO0, jst = 1:NO1
            DMst += 1
            DM_atom[jst,ist,1] = DM[1][DMnum+DMsum+DMst]
            DM_atom[jst,ist,2] = DM[2][DMnum+DMsum+DMst]
            DM_atom[jst,ist,3] = DM[3][DMnum+DMsum+DMst]
            DM_atom[jst,ist,4] = DM[4][DMnum+DMsum+DMst]
        end
        DMsum += DMst

        _Calc_Den2_nc!(ai_tempDGs, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], DM_atom)
        
        @inbounds for xyz = 1:GridN_Atom[atom]
            r = GridListAtom[atom][xyz]+1
            Density_Grid[1][r] += ai_tempDGs[xyz,1]
            Density_Grid[2][r] += ai_tempDGs[xyz,2]
            Density_Grid[3][r] += ai_tempDGs[xyz,3]
            Density_Grid[4][r] += ai_tempDGs[xyz,4]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[2], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[3], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[4], MPI.SUM, comm)
end


function diagonalize_nc_density!(Density_Grid)

    Ngrid = length(Density_Grid[1])

    @inbounds for i = 1:Ngrid

        Re11 = Density_Grid[1][i]
        Re22 = Density_Grid[2][i]
        Re12 = Density_Grid[3][i]
        Im12 = Density_Grid[4][i]

        Nup, Ndown, theta, phi = EulerAngle_Spin(Re11, Re22, Re12, Im12)

        Density_Grid[1][i] = Nup
        Density_Grid[2][i] = Ndown
        Density_Grid[3][i] = theta
        Density_Grid[4][i] = phi 
    end
end


function _Calc_Den2_nonpol!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:2:NumOLG-1

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1
        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1
        
        Sum0 = 0.0
        Sum1 = 0.0
        for ist = 1:NO0
            temp0 = 0.0
            temp1 = 0.0
            @inbounds for jst = 1:NO1
                tmp = DM[jst,ist]
                temp0  += Orbs_Grid2[Nh0][jst]*tmp
                temp1  += Orbs_Grid2[Nh1][jst]*tmp
            end
            Sum0  += Orbs_Grid1[Nc0][ist]*temp0
            Sum1  += Orbs_Grid1[Nc1][ist]*temp1
        end
        ai_tempDGs[Nc0] += Sum0
        ai_tempDGs[Nc1] += Sum1
    end


    Nog1 = 2*div(NumOLG, 2)
    rem_NumOLG = rem(NumOLG, 2)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum = 0.0
        for ist = 1:NO0
            temp = 0.0
            @inbounds for jst = 1:NO1
                temp += Orbs_Grid2[Nh][jst]*DM[jst,ist]
            end
            Sum += Orbs_Grid1[Nc][ist]*temp
        end
        ai_tempDGs[Nc] += Sum
    end
end


function _Calc_Den2_pol!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:2:NumOLG-1

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1

        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1

        Sum0_up = 0.0
        Sum1_up = 0.0

        Sum0_dn = 0.0
        Sum1_dn = 0.0

        for ist = 1:NO0
            temp0_up = 0.0
            temp1_up = 0.0

            temp0_dn = 0.0
            temp1_dn = 0.0

            @inbounds for jst = 1:NO1
                orbs2_0  = Orbs_Grid2[Nh0][jst]
                orbs2_1  = Orbs_Grid2[Nh1][jst]

                DM_up = DM[jst,ist,1]
                DM_dn = DM[jst,ist,2]
                
                temp0_up += orbs2_0*DM_up
                temp1_up += orbs2_1*DM_up

                temp0_dn += orbs2_0*DM_dn
                temp1_dn += orbs2_1*DM_dn
            end

            orbs1_0  = Orbs_Grid1[Nc0][ist]
            orbs1_1  = Orbs_Grid1[Nc1][ist]

            Sum0_up += orbs1_0*temp0_up
            Sum1_up += orbs1_1*temp1_up

            Sum0_dn += orbs1_0*temp0_dn
            Sum1_dn += orbs1_1*temp1_dn
        end
        ai_tempDGs[Nc0,1] += Sum0_up
        ai_tempDGs[Nc1,1] += Sum1_up

        ai_tempDGs[Nc0,2] += Sum0_dn
        ai_tempDGs[Nc1,2] += Sum1_dn
    end


    Nog1 = 2*div(NumOLG, 2)
    rem_NumOLG = rem(NumOLG, 2)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum_up = 0.0
        Sum_dn = 0.0
        for ist = 1:NO0
            temp_up = 0.0
            temp_dn = 0.0
            @inbounds for jst = 1:NO1
                temp_up += Orbs_Grid2[Nh][jst]*DM[jst,ist,1]
                temp_dn += Orbs_Grid2[Nh][jst]*DM[jst,ist,2]
            end
            Sum_up += Orbs_Grid1[Nc][ist]*temp_up
            Sum_dn += Orbs_Grid1[Nc][ist]*temp_dn
        end
        ai_tempDGs[Nc,1] += Sum_up
        ai_tempDGs[Nc,2] += Sum_dn
    end
end


function _Calc_Den2_nc!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:2:NumOLG-1

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1

        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1
        
        Sum0_uu = 0.0
        Sum1_uu = 0.0

        Sum0_dd = 0.0
        Sum1_dd = 0.0

        Sum0_ud_r = 0.0
        Sum1_ud_r = 0.0

        Sum0_ud_i = 0.0
        Sum1_ud_i = 0.0

        for ist = 1:NO0
            temp0_uu = 0.0
            temp1_uu = 0.0

            temp0_dd = 0.0
            temp1_dd = 0.0

            temp0_ud_r = 0.0
            temp1_ud_r = 0.0

            temp0_ud_i = 0.0
            temp1_ud_i = 0.0

            @inbounds for jst = 1:NO1

                orbs2_0  = Orbs_Grid2[Nh0][jst]
                orbs2_1  = Orbs_Grid2[Nh1][jst]

                DM_uu = DM[jst,ist,1]
                DM_dd = DM[jst,ist,2]
                DM_ud_r = DM[jst,ist,3]
                DM_ud_i = DM[jst,ist,4]
                
                temp0_uu += orbs2_0*DM_uu
                temp1_uu += orbs2_1*DM_uu

                temp0_dd += orbs2_0*DM_dd
                temp1_dd += orbs2_1*DM_dd

                temp0_ud_r += orbs2_0*DM_ud_r
                temp1_ud_r += orbs2_1*DM_ud_r

                temp0_ud_i += orbs2_0*DM_ud_i
                temp1_ud_i += orbs2_1*DM_ud_i
            end

            orbs1_0  = Orbs_Grid1[Nc0][ist]
            orbs1_1  = Orbs_Grid1[Nc1][ist]

            Sum0_uu += orbs1_0*temp0_uu
            Sum1_uu += orbs1_1*temp1_uu

            Sum0_dd += orbs1_0*temp0_dd
            Sum1_dd += orbs1_1*temp1_dd

            Sum0_ud_r += orbs1_0*temp0_ud_r
            Sum1_ud_r += orbs1_1*temp1_ud_r

            Sum0_ud_i += orbs1_0*temp0_ud_i
            Sum1_ud_i += orbs1_1*temp1_ud_i
        end
        ai_tempDGs[Nc0,1] += Sum0_uu
        ai_tempDGs[Nc1,1] += Sum1_uu

        ai_tempDGs[Nc0,2] += Sum0_dd
        ai_tempDGs[Nc1,2] += Sum1_dd

        ai_tempDGs[Nc0,3] += Sum0_ud_r
        ai_tempDGs[Nc1,3] += Sum1_ud_r

        ai_tempDGs[Nc0,4] += Sum0_ud_i
        ai_tempDGs[Nc1,4] += Sum1_ud_i
    end


    Nog1 = 2*div(NumOLG, 2)
    rem_NumOLG = rem(NumOLG, 2)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum_uu = 0.0
        Sum_dd = 0.0
        Sum_ud_r = 0.0
        Sum_ud_i = 0.0
        for ist = 1:NO0
            temp_uu = 0.0
            temp_dd = 0.0
            temp_ud_r = 0.0
            temp_ud_i = 0.0
            @inbounds for jst = 1:NO1
                orbs2 = Orbs_Grid2[Nh][jst]
                temp_uu += orbs2*DM[jst,ist,1]
                temp_dd += orbs2*DM[jst,ist,2]
                temp_ud_r += orbs2*DM[jst,ist,3]
                temp_ud_i += orbs2*DM[jst,ist,4]
            end
            orbs1 = Orbs_Grid1[Nc][ist]
            Sum_uu += orbs1*temp_uu
            Sum_dd += orbs1*temp_dd
            Sum_ud_r += orbs1*temp_ud_r
            Sum_ud_i += orbs1*temp_ud_i
        end
        ai_tempDGs[Nc,1] += Sum_uu
        ai_tempDGs[Nc,2] += Sum_dd
        ai_tempDGs[Nc,3] += Sum_ud_r
        ai_tempDGs[Nc,4] += Sum_ud_i
    end
end