function CWF2Wannier90(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, Enk, Cnk, Amnk)

    filename = cwf_setup.filename
    material = cwf_setup.material
    SpinPol = material.SpinPol
    ChemP = material.ChemP
    spinsize = cwf_setup.spinsize
    Dis_Energy = cwf_setup.Dis_Energy
    Nwann = cwf_setup.Ngsize
    Nkpt = kpoints.MPI_Nkpt
    

    MinN, MaxN = Find_MinN_MaxN(material, Dis_Energy, Enk)
    BANDNUM = MaxN - MinN + 1


    if SpinPol ∈ ("off", "on")
        println("<Generate_Amnk_Col>")
        Generate_Amnk_Col(filename, spinsize, MinN, MaxN, Nkpt, Nwann, Amnk)
    else
        println("<Generate_Amnk_NonCol>")
        Generate_Amnk_NonCol(filename, MinN, MaxN, Nkpt, Nwann, Amnk)
    end

    
    if SpinPol ∈ ("off", "on")
        println("<Generate_Mmnkb_Col>")
        Generate_Mmnkb_Col(filename, cwf_setup, kpoints, MinN, MaxN, Cnk)
    else
        println("<Generate_Mmnkb_NonCol>")
        Generate_Mmnkb_NonCol(filename, cwf_setup, kpoints, MinN, MaxN, Cnk)
    end
    

    if SpinPol ∈ ("off", "on")
        println("<Generate_eig_Col>")
        Generate_eig_Col(filename, ChemP, spinsize, MinN, MaxN, Nkpt, Enk)
    else
        println("<Generate_eig_NonCol>")
        Generate_eig_NonCol(filename, ChemP, MinN, MaxN, Nkpt, Enk)
    end


    Write_win(filename, BANDNUM, kpoints, cwf_setup)
end


function Generate_Amnk_Col(filename::String, spinsize, MinN, MaxN, Nkpt, Nwann, Amnk)
    
    BANDNUM = MaxN - MinN + 1
    for spin = 1:spinsize
        data = open(filename*"_$spin.amn", "w")
        @printf(data, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
        @printf(data, "%d %d %d\n", BANDNUM, Nkpt, Nwann)
        for ik = 1:Nkpt
            for m = 1:Nwann, μ = 1:BANDNUM
                @printf(data, "%5d %5d %5d %22.14f %22.14f\n", μ, m, ik, real(Amnk[spin][ik][μ+MinN-1][m]), imag(Amnk[spin][ik][μ+MinN-1][m]))
            end
        end
        close(data)
    end
end


function Generate_Amnk_NonCol(filename::String, MinN, MaxN, Nkpt, Nwann, Amnk)
    
    BANDNUM = MaxN - MinN + 1
    data = open(filename*".amn", "w")
    @printf(data, "2nd line: BANDNUM, KPTNUM, WANNUM, Nexts: band, wan, k, ReAmn, ImAmn\n")
    @printf(data, "%d %d %d\n", BANDNUM, Nkpt, Nwann)

    for ik = 1:Nkpt
        for m = 1:Nwann, μ = 1:BANDNUM
            @printf(data, "%5d %5d %5d %22.14f %22.14f\n", μ, m, ik, real(Amnk[1][ik][μ+MinN-1][m]), imag(Amnk[1][ik][μ+MinN-1][m]))
        end
    end
    close(data)
end


function Generate_Mmnkb_Col(filename::String, cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, MinN, MaxN, Cnk)

    material = cwf_setup.material
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    spinsize = cwf_setup.spinsize
    Nwann = cwf_setup.Ngsize

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = kpoints.MPI_Nkpt
    kpts = kpoints.MPI_kpts

    BANDNUM = MaxN - MinN + 1
    disentangle = ifelse(BANDNUM > Nwann, true, false)

    MAXSHELL = 30
    println("<Set_MLWF_kgrid>")
    _, tot_bvector, bvector, frac_bv_int, frac_bv, kplusb, _ = Set_MLWF_kgrid(Latvecs, kmesh, MAXSHELL)

    println("<Set_OLPexp>")
    OLPexp = Set_OLPexp(cwf_setup, tot_bvector, bvector)


    println("Calculating Mmnk...")

    k1 = zeros(Float64, 3)
    k2 = zeros(Float64, 3)
    dk = zeros(Float64, 3)
    Mmnkb = zeros(ComplexF64, fsize, fsize)
    for spin = 1:spinsize

        println("spin: $spin")

        data = open(filename*"_$spin.mmn", "w")
        @printf(data, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
        @printf(data, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

        for ik = 1:Nkpt
            @printf("\t%d/%d\n", ik, Nkpt)

            @. k1 = kpts[ik]
            i = div(ik-1, kmesh2*kmesh3)
            j = div(ik-1 - i*kmesh2*kmesh3, kmesh3)
            k = ik-1 - i*kmesh2*kmesh3 - j*kmesh3
            for ib = 1:tot_bvector

                kk = kplusb[ik][ib]+1
                @. k2 = kpts[ik] + frac_bv[ib]
                @. dk = frac_bv[ib]
                
                Calc_Mmnkb_Col!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], Cnk[spin][ik], Cnk[spin][kk])

                ii = ((i+frac_bv_int[ib][1])%kmesh1 + kmesh1)%kmesh1
                jj = ((j+frac_bv_int[ib][2])%kmesh2 + kmesh2)%kmesh2
                kk = ((k+frac_bv_int[ib][3])%kmesh3 + kmesh3)%kmesh3

                ikb = ii*kmesh2*kmesh3 + jj*kmesh3 + kk

                g1 = div((i+frac_bv_int[ib][1])-ii, kmesh1)
                g2 = div((j+frac_bv_int[ib][2])-jj, kmesh2)
                g3 = div((k+frac_bv_int[ib][3])-kk, kmesh3)
                
                @printf(data, "%d %d  %d %d %d\n", ik, ikb+1, g1, g2, g3)
                for ν = 1:BANDNUM, μ = 1:BANDNUM
                    @printf(data, " %22.14f %22.14f\n", real(Mmnkb[μ,ν]), imag(Mmnkb[μ,ν]))
                end
            end
        end
        close(data)
    end
end


function Generate_Mmnkb_NonCol(filename::String, cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, MinN, MaxN, Cnk)

    material = cwf_setup.material
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    Nwann = cwf_setup.Ngsize

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = kpoints.MPI_Nkpt
    kpts = kpoints.MPI_kpts

    BANDNUM = MaxN - MinN + 1
    disentangle = ifelse(BANDNUM > Nwann, true, false)

    MAXSHELL = 30
    println("<Set_MLWF_kgrid>")
    _, tot_bvector, bvector, frac_bv_int, frac_bv, kplusb, _ = Set_MLWF_kgrid(Latvecs, kmesh, MAXSHELL)

    println("<Set_OLPexp>")
    OLPexp = Set_OLPexp(cwf_setup, tot_bvector, bvector)


    println("Calculating Mmnk...")

    k1 = zeros(Float64, 3)
    k2 = zeros(Float64, 3)
    dk = zeros(Float64, 3)
    Mmnkb = zeros(ComplexF64, fsize, fsize)

    data = open(filename*".mmn", "w")
    @printf(data, "2nd line: BANDNUM, KPTNUM, NUMB, Nexts: KPTNUM x NUMBband elements block\n")
    @printf(data, "%d %d %d\n", BANDNUM, Nkpt, tot_bvector)

    for ik = 1:Nkpt
        @printf("\t%d/%d\n", ik, Nkpt)

        @. k1 = kpts[ik]
        i = div(ik-1, kmesh2*kmesh3)
        j = div(ik-1 - i*kmesh2*kmesh3, kmesh3)
        k = ik-1 - i*kmesh2*kmesh3 - j*kmesh3
        
        for ib = 1:tot_bvector

            kk = kplusb[ik][ib]+1
            @. k2 = kpts[ik] + frac_bv[ib]
            @. dk = frac_bv[ib]
                
            Calc_Mmnkb_NonCol!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp[ib], Cnk[1][ik], Cnk[1][kk])

            ii = ((i+frac_bv_int[ib][1])%kmesh1 + kmesh1)%kmesh1
            jj = ((j+frac_bv_int[ib][2])%kmesh2 + kmesh2)%kmesh2
            kk = ((k+frac_bv_int[ib][3])%kmesh3 + kmesh3)%kmesh3

            ikb = ii*kmesh2*kmesh3 + jj*kmesh3 + kk

            g1 = div((i+frac_bv_int[ib][1])-ii, kmesh1)
            g2 = div((j+frac_bv_int[ib][2])-jj, kmesh2)
            g3 = div((k+frac_bv_int[ib][3])-kk, kmesh3)
                
            @printf(data, "%d %d  %d %d %d\n", ik, ikb+1, g1, g2, g3)
            for ν = 1:BANDNUM, μ = 1:BANDNUM
                @printf(data, " %22.14f %22.14f\n", real(Mmnkb[μ,ν]), imag(Mmnkb[μ,ν]))
            end
        end
    end
    close(data)
end


function Generate_eig_Col(filename::String, ChemP, spinsize, MinN, MaxN, Nkpt, Enk)

    BANDNUM = MaxN - MinN + 1
    for spin = 1:spinsize
        data = open(filename*"_$spin.eig", "w")
        for ik = 1:Nkpt, μ = 1:BANDNUM
            @printf(data, "%5d %5d %22.14f\n", μ, ik, (Enk[spin,μ+MinN-1,ik]-ChemP)*27.2113845)
        end
        close(data)
    end
end


function Generate_eig_NonCol(filename::String, ChemP, MinN, MaxN, Nkpt, Enk)

    BANDNUM = MaxN - MinN + 1
    data = open(filename*".eig", "w")
    for ik = 1:Nkpt, μ = 1:BANDNUM
        @printf(data, "%5d %5d %22.14f\n", μ, ik, (Enk[1,μ+MinN-1,ik]-ChemP)*27.2113845)
    end
    close(data)
end


function Write_win(filename::String, BANDNUM, kpoints::KPoints, cwf_setup::Union{CWF_Setup,CWF_Setup_MO})

    material = cwf_setup.material
    SpinPol = material.SpinPol
    Natom = material.Natom
    Latvecs = material.Latvecs
    Recvecs = 2*pi*inv(Latvecs')
    Gxyz_AU = material.Gxyz
    Atoms_symbol = material.Atoms_symbol
    spinsize = cwf_setup.spinsize
    WANNUM = cwf_setup.Ngsize

    kmesh = kpoints.kmesh
    Nkpt = kpoints.Nkpt
    kpts = kpoints.MPI_kpts


    Gxyz_frac = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Gxyz_frac[atom] = zeros(Float64, 3)
    end

    for atom = 1:Natom
        Gxyz_frac[atom][1] = dot(Gxyz_AU[atom], Recvecs[1,:])*0.5/pi
        Gxyz_frac[atom][2] = dot(Gxyz_AU[atom], Recvecs[2,:])*0.5/pi
        Gxyz_frac[atom][3] = dot(Gxyz_AU[atom], Recvecs[3,:])*0.5/pi
                
        for i = 1:3
            tmp = floor(Int64, Gxyz_frac[atom][i])
            if Gxyz_frac[atom][i] > 1.0
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]-tmp)
            elseif Gxyz_frac[atom][i] < -1e-13
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]+abs(tmp)+1)
            end
        end
    end


    for spin = 1:spinsize
        if SpinPol ∈ ("off", "on")
            data = open(filename*"_$spin.win", "w")
        else
            data = open(filename*".win", "w")
        end

        @printf(data, "num_bands %d\n", BANDNUM)
        @printf(data, "num_wann %d\n", WANNUM)
        @printf(data, "\n")
        @printf(data, "!spinors T\n")
        @printf(data, "\n")
        @printf(data, "dis_num_iter 0\n")
        @printf(data, "num_iter 0\n")
        @printf(data, "\n")
        @printf(data, "write_rmn  = .false.\n")
        @printf(data, "write_r2mn = .false.\n")
        @printf(data, "write_hr = .false.\n")
        @printf(data, "write_u_matrices = .false.\n")
        @printf(data, "\n")
        @printf(data, "! kmesh_tol = 0.00001\n")
        @printf(data, "\n")
        @printf(data, "bands_plot F\n")
        @printf(data, "bands_num_points 200\n")
        @printf(data, "bands_plot_format gnuplot\n")
        @printf(data, "!begin kpoint_path\n")
        @printf(data, "!end kpoint_path\n")
        @printf(data, "\n")
        @printf(data, "begin unit_cell_cart\n")
        @printf(data, "Ang\n")
        for i = 1:3
            @printf(data, "%lf %lf %lf\n", Latvecs[i,1]*0.529177249, Latvecs[i,2]*0.529177249, Latvecs[i,3]*0.529177249)
        end
        @printf(data, "end unit_cell_cart\n")
        @printf(data, "\n")
        @printf(data, "begin atoms_frac\n")
        for atom = 1:Natom
            @printf(data, "%4s %18.14f %18.14f %18.14f\n", Atoms_symbol[atom], Gxyz_frac[atom][1], Gxyz_frac[atom][2], Gxyz_frac[atom][3])
        end
        @printf(data, "end atoms_frac\n")
        @printf(data, "\n")
        @printf(data, "mp_grid %d %d %d\n", kmesh[1], kmesh[2], kmesh[3])
        @printf(data, "\n")
        @printf(data, "begin_kpoints\n")
        for ik = 1:Nkpt
            @printf(data, "%12.8f%12.8f%12.8f\n", kpts[ik][1], kpts[ik][2], kpts[ik][3])
        end
        @printf(data, "end_kpoints\n")
        close(data)
    end
end


function Calc_Mmnkb_Col!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, Cnk1, Cnk2)
    
    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]

    fill!(Mmnkb, 0.0)
    
    band_offset = MinN - 1

    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        cell  = ncn[atom][Rn]+1
        l1, l2, l3 = atv_ijk[cell]
        NO0  = Total_NumOrbs[atom]
        NO1  = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        kRn1 = -2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn
        kRn2 =  2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn
        ex1 = cis(kRn2)
        ex2 = cis(kRn1)

        _OLPexp = OLPexp[atom][Rn]

        for jst = 1:NO1, ist = 1:NO0
            aidx = Anum + ist
            bidx = Bnum + jst

            sij = _OLPexp[ist][jst]
            α1 = 0.5 * ex1 * sij
            α2 = 0.5 * ex2 * sij
            
            for ν = 1:BANDNUM
                νidx = ν + band_offset

                c2b = Cnk2[bidx,νidx]
                c2a = Cnk2[aidx,νidx]

                @inbounds for μ = 1:BANDNUM
                    μidx = μ + band_offset
                    c1a = Cnk1[aidx, μidx]
                    c1b = Cnk1[bidx, μidx]

                    Mmnkb[μ,ν] += α1*conj(c1a)*c2b + α2*conj(c1b)*c2a
                end
            end
        end
    end
end


function Calc_Mmnkb_NonCol!(Mmnkb, MinN, BANDNUM, Natom, Gxyz, FNAN, natn, ncn, Total_NumOrbs, MP, atv_ijk, Recvecs, k1, k2, dk, OLPexp, Cnk1, Cnk2)
    
    fsize = sum(Total_NumOrbs)
    k1a, k1b, k1c = k1
    k2a, k2b, k2c = k2
    dkx = dk[1]*Recvecs[1,1] + dk[2]*Recvecs[2,1] + dk[3]*Recvecs[3,1]
    dky = dk[1]*Recvecs[1,2] + dk[2]*Recvecs[2,2] + dk[3]*Recvecs[3,2]
    dkz = dk[1]*Recvecs[1,3] + dk[2]*Recvecs[2,3] + dk[3]*Recvecs[3,3]

    fill!(Mmnkb, 0.0)

    band_offset = MinN - 1
    
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        cell  = ncn[atom][Rn]+1
        l1, l2, l3 = atv_ijk[cell]
        NO0  = Total_NumOrbs[atom]
        NO1  = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        dkRn = dkx*Gxyz[atom][1] + dky*Gxyz[atom][2] + dkz*Gxyz[atom][3]
        kRn1 = -2*pi*(k1a*l1 + k1b*l2 + k1c*l3) - dkRn
        kRn2 =  2*pi*(k2a*l1 + k2b*l2 + k2c*l3) - dkRn
        ex1 = cis(kRn2)
        ex2 = cis(kRn1)

        _OLPexp = OLPexp[atom][Rn]

        for jst = 1:NO1, ist = 1:NO0
            aidx1 = Anum + ist
            aidx2 = Anum + ist + fsize
            bidx1 = Bnum + jst
            bidx2 = Bnum + jst + fsize

            sij = _OLPexp[ist][jst]
            α1 = 0.5 * ex1 * sij
            α2 = 0.5 * ex2 * sij
            
            for ν = 1:BANDNUM
                νidx = ν + band_offset

                c2b1 = Cnk2[bidx1,νidx]
                c2b2 = Cnk2[bidx2,νidx]
                c2a1 = Cnk2[aidx1,νidx]
                c2a2 = Cnk2[aidx2,νidx]
                

                @inbounds for μ = 1:BANDNUM
                    μidx = μ + band_offset
                    c1a1 = Cnk1[aidx1, μidx]
                    c1a2 = Cnk1[aidx2, μidx]
                    c1b1 = Cnk1[bidx1, μidx]
                    c1b2 = Cnk1[bidx2, μidx]

                    Mmnkb[μ,ν] += α1*conj(c1a1)*c2b1 + α2*conj(c1b1)*c2a1
                    Mmnkb[μ,ν] += α1*conj(c1a2)*c2b2 + α2*conj(c1b2)*c2a2
                end
            end
        end
    end
end