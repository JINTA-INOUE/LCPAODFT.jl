function Calc_partialDM_Crystal_Collinear!(parDM, material::LCPAO_model, kpoints::KPoints, Dis_Enegry, weight_type, Enk, Cnk)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nspin = material.Nspin
    SpinPol = material.SpinPol
    Natom = material.Natom
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol=="off", 1, 2)

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    Nloop = sum(FNAN.+1)
    OneD2atom = zeros(Int32, Nloop)
    OneD2FNAN = zeros(Int32, Nloop)
    OneD2natn = zeros(Int32, Nloop)

    counts = 1
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        OneD2atom[counts] = atom
        OneD2FNAN[counts] = Rn
        OneD2natn[counts] = natn[atom][Rn]

        counts += 1
    end

    myrange = split_evenly(1:Nloop, nprocs)
    MPI_size = length(myrange[myrank+1])
    MPI_atom = OneD2atom[myrange[myrank+1]]
    MPI_natn = OneD2natn[myrange[myrank+1]]

    myHsize = 0
    for loop = 1:MPI_size, _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:Total_NumOrbs[MPI_natn[loop]]
        myHsize += 1
    end
    Total_Hsize = MPI.Allreduce(myHsize, MPI.SUM, comm)


    DM_1D = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        DM_1D[spin] = zeros(Float64, Total_Hsize)
    end


    tmpCnk = zeros(ComplexF64, fsize)
    @inbounds for spin = 1:spinsize, ik = 1:MPI_Nkpt, μ = 1:fsize
        if weight_type == "fermi"
            FF = CWF_weight(Enk[spin][ik][μ], ChemP, Dis_Enegry)
        elseif weight_type == "poly"
            FF = CWF_weight2(Enk[spin][ik][μ], ChemP, Dis_Enegry)
        else
            error("please check weight_type")
        end
        F = sqrt(FF)
        @views mul!(tmpCnk, F, Cnk[spin][ik][:,μ])
        @. @views Cnk[spin][ik][:,μ] = tmpCnk
    end

    
    ctemp = zeros(ComplexF64, fsize, fsize)
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        @. ctemp = Cnk[spin][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            coskRn = cos(2*pi*kRn)/AllNkpt
            sinkRn = sin(2*pi*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                tmp_re = 0.0
                tmp_im = 0.0
                @inbounds for μ = 1:fsize
                    a = ctemp[Anum+ist,μ]
                    b = ctemp[Bnum+jst,μ]
                    tmp_re += real(a)*real(b) + imag(a)*imag(b)
                    tmp_im += imag(a)*real(b) - real(a)*imag(b)
                end
                DMst += 1
                DM_1D[spin][DMst] += tmp_re*coskRn + tmp_im*sinkRn
            end
        end
    end


    for spin = 1:Nspin
        MPI.Allreduce!(DM_1D[spin], MPI.SUM, comm)
    end


    for spin = 1:Nspin
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            DMst += 1
            parDM[spin][atom][Rn][ist][jst] = DM_1D[spin][DMst]
        end
    end
end


function Set_CWF_Guiding_MOs(cwf_setup::CWF_Setup_MO, kpoints::KPoints, Enk, Cnk)

    material = cwf_setup.material
    spinsize = cwf_setup.spinsize
    Natom = material.Natom
    SpinPol = material.SpinPol
    TCpyCell = material.TCpyCell
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol
    OLP = material.OLP
    Hks = material.Hks
    parDM = material.DM

    Num_CWF_Grouped_Atoms = cwf_setup.Num_CWF_Grouped_Atoms
    CWF_Grouped_Atoms_EachNum = cwf_setup.CWF_Grouped_Atoms_EachNum
    Num_CWF_MOs_Group = cwf_setup.Num_CWF_MOs_Group
    CWF_Grouped_Atoms = cwf_setup.CWF_Grouped_Atoms
    CWF_Total_NumOrbs = cwf_setup.CWF_Total_NumOrbs
    MP3 = cwf_setup.MP3
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type

    if SpinPol == "nc"
        error("not support Noncollinear case.")
    end



    Calc_partialDM_Crystal_Collinear!(parDM, material, kpoints, Dis_Energy, weight_type, Enk, Cnk)



    CpyCell = div(Int64(cbrt(TCpyCell+1))-1,2)
    RMI = Get_RMI(Natom, CpyCell, FNAN, natn, ncn)



    RMI0 = Vector{Vector{Vector{Int32}}}(undef, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        RMI0[gidx] = Vector{Vector{Int32}}(undef, CWF_Grouped_Atoms_EachNum[gidx])
        for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            RMI0[gidx][Lidx] = zeros(Int32, (FNAN[atom]+1)^2)
        end
    end

    for gidx = 1:Num_CWF_Grouped_Atoms, Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
        atom = CWF_Grouped_Atoms[gidx][Lidx]
        for i = 1:FNAN[atom]+1, j = 1:FNAN[atom]+1
            RMI0[gidx][Lidx][(i-1)*(FNAN[atom]+1)+j] = RMI[atom][i][j]
        end
    end


    CWF_Guiding_MOs = Vector{Vector{Vector{Float64}}}(undef, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        CWF_Guiding_MOs[gidx] = Vector{Vector{Float64}}(undef, Num_CWF_MOs_Group[gidx])
        for i = 1:Num_CWF_MOs_Group[gidx]
            CWF_Guiding_MOs[gidx][i] = zeros(Float64, CWF_Total_NumOrbs[gidx])
        end
    end


    # data = open("MOs_info.dat", "w")

    for gidx = 1:Num_CWF_Grouped_Atoms

        DMS = zeros(Float64, CWF_Total_NumOrbs[gidx], CWF_Total_NumOrbs[gidx])
        DMS2 = zeros(Float64, CWF_Total_NumOrbs[gidx], CWF_Total_NumOrbs[gidx])

        Smo2 = Vector{Vector{Vector{Vector{Float64}}}}(undef, CWF_Grouped_Atoms_EachNum[gidx])
        DMmo2 = Vector{Vector{Vector{Vector{Float64}}}}(undef, CWF_Grouped_Atoms_EachNum[gidx])
        Hmo2 = Vector{Vector{Vector{Vector{Float64}}}}(undef, CWF_Grouped_Atoms_EachNum[gidx])
        for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            Smo2[Lidx] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            DMmo2[Lidx] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            Hmo2[Lidx] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                Smo2[Lidx][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                DMmo2[Lidx][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                Hmo2[Lidx][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    Smo2[Lidx][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                    DMmo2[Lidx][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                    Hmo2[Lidx][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end

        for spin = 1:spinsize, Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            for Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
                Smo2[Lidx][Rn][ist][jst] = OLP[atom][Rn][ist][jst]
                DMmo2[Lidx][Rn][ist][jst] += parDM[spin][atom][Rn][ist][jst]/spinsize
                Hmo2[Lidx][Rn][ist][jst] += Hks[spin][atom][Rn][ist][jst]/spinsize
            end
        end
    

        for Lidx1 = 1:CWF_Grouped_Atoms_EachNum[gidx]
            atom1 = CWF_Grouped_Atoms[gidx][Lidx1]
            Anum = MP3[gidx][Lidx1]
            for Lidx2 = 1:CWF_Grouped_Atoms_EachNum[gidx]
                atom2 = CWF_Grouped_Atoms[gidx][Lidx2]
                Bnum = MP3[gidx][Lidx2]

                # find Rn3
                Rn3 = -1
                for Rn = 0:FNAN[atom1]
                    atom3 = natn[atom1][Rn+1]
                    RnC = ncn[atom1][Rn+1]+1
                    l1, l2, l3 = atv_ijk[RnC]
                    if atom2 == atom3 && l1==0 && l2==0 && l3==0
                        Rn3 = Rn
                    end
                end

                if Rn3 >= 0
                    for Rn = 1:FNAN[atom1]+1
                        atom3 = natn[atom1][Rn]
                        Rn2 = RMI0[gidx][Lidx1][Rn3*(FNAN[atom1]+1)+Rn]
                        if Rn2 >= 0
                            for ist = 1:Total_NumOrbs[atom1], jst = 1:Total_NumOrbs[atom2]
                                Sum = 0.0
                                for kst = 1:Total_NumOrbs[atom3]
                                    Sum += DMmo2[Lidx1][Rn][ist][kst]*Smo2[Lidx2][Rn2+1][jst][kst]
                                end
                                DMS[Anum+ist,Bnum+jst] += Sum
                            end
                        end
                    end
                end
            end
        end


        OLPgidx = zeros(Float64, CWF_Total_NumOrbs[gidx], CWF_Total_NumOrbs[gidx])
        Hgidx = zeros(Float64, CWF_Total_NumOrbs[gidx], CWF_Total_NumOrbs[gidx])
        WI = zeros(Float64, CWF_Total_NumOrbs[gidx])

        for ist = 1:CWF_Total_NumOrbs[gidx], jst = 1:CWF_Total_NumOrbs[gidx]
            Sum = 0.0
            for kst = 1:CWF_Total_NumOrbs[gidx]
                Sum += DMS[kst,ist]*DMS[kst,jst]
            end
            DMS2[ist,jst] = -Sum
        end

        Eval, Evec = eigen(Symmetric(DMS2, :L))
        for μ = 1:CWF_Total_NumOrbs[gidx]
            Eval[μ] = sqrt(abs(Eval[μ]))
        end

    

        for Lidx1 = 1:CWF_Grouped_Atoms_EachNum[gidx]
            GA_AN = CWF_Grouped_Atoms[gidx][Lidx1]
            Anum = MP3[gidx][Lidx1]

            for Rn = 1:FNAN[GA_AN]+1
                GC_AN = natn[GA_AN][Rn]
                RnC = ncn[GA_AN][Rn]+1
                
                if atv_ijk[RnC] == [0, 0, 0]
                    for Lidx2 = 1:CWF_Grouped_Atoms_EachNum[gidx]
                        GB_AN = CWF_Grouped_Atoms[gidx][Lidx2]
                        Bnum = MP3[gidx][Lidx2]

                        if GB_AN == GC_AN
                            for ist = 1:Total_NumOrbs[GA_AN], jst = 1:Total_NumOrbs[GB_AN]
                                OLPgidx[Anum+ist,Bnum+jst] = Smo2[Lidx1][Rn][ist][jst]
                                Hgidx[Anum+ist,Bnum+jst] = Hmo2[Lidx1][Rn][ist][jst]
                            end
                        end
                    end
                end
            end
        end


        for p = 1:CWF_Total_NumOrbs[gidx]
            Sum = 0.0
            for ist = 1:CWF_Total_NumOrbs[gidx], jst = 1:CWF_Total_NumOrbs[gidx]
                Sum += OLPgidx[ist,jst]*Evec[ist,p]*Evec[jst,p]
            end

            Sum = 1/sqrt(abs(Sum))
            for ist = 1:CWF_Total_NumOrbs[gidx]
                Evec[ist,p] *= Sum
            end

            Sum = 0.0
            for ist = 1:CWF_Total_NumOrbs[gidx], jst = 1:CWF_Total_NumOrbs[gidx]
                Sum += Hgidx[ist,jst]*Evec[ist,p]*Evec[jst,p]
            end
            WI[p] = Sum
        end

        #=
        @printf(data, "# Group index=%2d\n", gidx)
        @printf(data,"# 1st column: serial number\n")
        @printf(data,"# 2nd column: on-site energy (eV) relative to chemical potential\n")
        @printf(data,"# 3rd column: population\n")
        for p = 1:CWF_Total_NumOrbs[gidx]
            @printf(data,"%2d %18.12f %18.12f\n",p, (WI[p]-ChemP)*27.2113845, Eval[p])
        end=#


        for p = 1:Num_CWF_MOs_Group[gidx]
            for j = 1:CWF_Total_NumOrbs[gidx]
                CWF_Guiding_MOs[gidx][p][j] = Evec[j,p]
            end
        end
    end
    # close(data)


    return CWF_Guiding_MOs
end