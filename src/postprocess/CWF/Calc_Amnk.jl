@timeit timer "Calc_Amnk" function Calc_Amnk(MinN, BANDNUM, Enk, Cnk, cwf_setup::CWF_Setup, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    SpinPol = material.SpinPol
    spinsize = cwf_setup.spinsize
    Ngsize = cwf_setup.Ngsize
    MPI_Nkpt = kpoints.MPI_Nkpt


    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, BANDNUM)
            for μ = 1:BANDNUM
                Amnk[spin][ik][μ] = zeros(ComplexF64, Ngsize)
            end
        end
    end


    if SpinPol == "off"
        Calc_Amnk_nonpol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup, kpoints)
    elseif SpinPol == "on"
        Calc_Amnk_pol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup, kpoints)
    else
        Calc_Amnk_NonCol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup, kpoints)
    end


    return Amnk
end


function Calc_Amnk_nonpol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup::CWF_Setup, kpoints::KPoints)

    material = cwf_setup.material
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    OLP = material.OLP
    GNatom = cwf_setup.GNatom
    Guide_index = cwf_setup.Guide_index
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        for patom = 1:GNatom, proj in Guide_index[patom]    
            pst += 1
            Sum = ComplexF64(0.0, 0.0)
            for Rn = 1:FNAN[patom]+1
    
                jatom = natn[patom][Rn]
                cell = ncn[patom][Rn]+1
                kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
                ex = cispi(-2*kRn)
                Bnum = MP[jatom]
                tmp1 = ComplexF64(0.0, 0.0)
                for jβ = 1:Total_NumOrbs[jatom]
                    tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ+MinN-1])*OLP[patom][Rn][proj][jβ]
                end
                Sum += tmp1*ex
            end
    
            if weight_type == "fermi"
                weight = CWF_weight(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            else    # weight_type == "poly"
                weight = CWF_weight2(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            end
            Amnk[1][ik][μ][pst] = weight * Sum
        end
    end
end


function Calc_Amnk_pol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup::CWF_Setup, kpoints::KPoints)

    material = cwf_setup.material
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    OLP = material.OLP
    GNatom = cwf_setup.GNatom
    Guide_index = cwf_setup.Guide_index
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        for patom = 1:GNatom, proj in Guide_index[patom]
            pst += 1
            Sum1 = ComplexF64(0.0, 0.0)
            Sum2 = ComplexF64(0.0, 0.0)
            for Rn = 1:FNAN[patom]+1
                jatom = natn[patom][Rn]
                cell = ncn[patom][Rn]+1
                kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
                ex = cispi(-2*kRn)
                Bnum = MP[jatom]
                tmp1 = ComplexF64(0.0, 0.0)
                tmp2 = ComplexF64(0.0, 0.0)
                for jβ = 1:Total_NumOrbs[jatom]
                    tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ+MinN-1])*OLP[patom][Rn][proj][jβ]
                    tmp2 += conj(Cnk[2][ik][Bnum+jβ,μ+MinN-1])*OLP[patom][Rn][proj][jβ]
                end
                Sum1 += tmp1*ex
                Sum2 += tmp2*ex
            end
                
            if weight_type == "fermi"
                weight1 = CWF_weight(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
                weight2 = CWF_weight(Enk[2][ik][μ+MinN-1], ChemP, Dis_Energy)
            else    # weight_type == "poly"
                weight1 = CWF_weight2(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
                weight2 = CWF_weight2(Enk[2][ik][μ+MinN-1], ChemP, Dis_Energy)
            end

            Amnk[1][ik][μ][pst] = weight1 * Sum1
            Amnk[2][ik][μ][pst] = weight2 * Sum2
        end
    end
end


function Calc_Amnk_NonCol!(MinN, BANDNUM, Amnk, Enk, Cnk, cwf_setup::CWF_Setup, kpoints::KPoints)

    material = cwf_setup.material
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    ChemP = material.ChemP
    OLP = material.OLP

    GNatom = cwf_setup.GNatom
    gsize = cwf_setup.gsize
    Guide_index = cwf_setup.Guide_index
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        for patom = 1:GNatom, proj in Guide_index[patom]    
            pst += 1
            Sum1 = ComplexF64(0.0, 0.0)
            Sum2 = ComplexF64(0.0, 0.0)
            for Rn = 1:FNAN[patom]+1
    
                jatom = natn[patom][Rn]
                cell = ncn[patom][Rn]+1
                kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
                ex = cispi(-2*kRn)
                Bnum = MP[jatom]
                tmp1 = ComplexF64(0.0, 0.0)
                tmp2 = ComplexF64(0.0, 0.0)
                for jβ = 1:Total_NumOrbs[jatom]
                    tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ+MinN-1])*OLP[patom][Rn][proj][jβ]
                    tmp2 += conj(Cnk[1][ik][fsize+Bnum+jβ,μ+MinN-1])*OLP[patom][Rn][proj][jβ]
                end
                Sum1 += tmp1*ex
                Sum2 += tmp2*ex
            end
    
            if weight_type == "fermi"
                weight = CWF_weight(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            else    # weight_type == "poly"
                weight = CWF_weight2(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            end
            Amnk[1][ik][μ][pst] = weight * Sum1
            Amnk[1][ik][μ][gsize+pst] = weight * Sum2
        end
    end
end


@timeit timer "Calc_Amnk" function Calc_Amnk(MinN, BANDNUM, Enk, Cnk, cwf_setup::CWF_Setup_MO, kpoints::KPoints, CWF_Guiding_MOs)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    FNAN = material.FNAN
    natn = material.natn
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol
    OLP = material.OLP

    spinsize = cwf_setup.spinsize
    Num_CWF_Grouped_Atoms = cwf_setup.Num_CWF_Grouped_Atoms
    Num_CWF_MOs_Group = cwf_setup.Num_CWF_MOs_Group
    CWF_Grouped_Atoms_EachNum = cwf_setup.CWF_Grouped_Atoms_EachNum
    CWF_Grouped_Atoms = cwf_setup.CWF_Grouped_Atoms
    MP3 = cwf_setup.MP3
    Ngsize = cwf_setup.Ngsize
    MPI_Nkpt = kpoints.MPI_Nkpt


    OLPproj = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Num_CWF_Grouped_Atoms)
    for gidx = 1:Num_CWF_Grouped_Atoms
        OLPproj[gidx] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Num_CWF_MOs_Group[gidx])
        for p = 1:Num_CWF_MOs_Group[gidx]
            OLPproj[gidx][p] = Vector{Vector{Vector{Float64}}}(undef, CWF_Grouped_Atoms_EachNum[gidx])
            for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
                atom = CWF_Grouped_Atoms[gidx][Lidx]
                OLPproj[gidx][p][Lidx] = Vector{Vector{Float64}}(undef, FNAN[atom]+1)
                for Rn = 1:FNAN[atom]+1
                    OLPproj[gidx][p][Lidx][Rn] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end


    for gidx = 1:Num_CWF_Grouped_Atoms, p = 1:Num_CWF_MOs_Group[gidx]
        for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
            atom = CWF_Grouped_Atoms[gidx][Lidx]
            Anum = MP3[gidx][Lidx]

            for Rn = 1:FNAN[atom]+1, jst = 1:Total_NumOrbs[natn[atom][Rn]]
                Sum = 0.0
                for ist = 1:Total_NumOrbs[atom]
                    Sum += OLP[atom][Rn][ist][jst]*CWF_Guiding_MOs[gidx][p][Anum+ist]
                end
                OLPproj[gidx][p][Lidx][Rn][jst] = Sum
            end
        end
    end



    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, BANDNUM)
            for μ = 1:BANDNUM
                Amnk[spin][ik][μ] = zeros(ComplexF64, Ngsize)
            end
        end
    end


    if SpinPol == "off"
        Calc_Amnk_nonpol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup, kpoints)
    elseif SpinPol == "on"
        Calc_Amnk_pol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup, kpoints)
    else
        Calc_Amnk_NonCol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup, kpoints)
    end


    return Amnk
end


function Calc_Amnk_nonpol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup::CWF_Setup_MO, kpoints::KPoints)

    material = cwf_setup.material
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type
    Num_CWF_Grouped_Atoms = cwf_setup.Num_CWF_Grouped_Atoms
    Num_CWF_MOs_Group = cwf_setup.Num_CWF_MOs_Group
    CWF_Grouped_Atoms_EachNum = cwf_setup.CWF_Grouped_Atoms_EachNum
    CWF_Grouped_Atoms = cwf_setup.CWF_Grouped_Atoms

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        for gidx = 1:Num_CWF_Grouped_Atoms, p = 1:Num_CWF_MOs_Group[gidx]
            pst += 1
            Sum = ComplexF64(0.0, 0.0)
            for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
                atom = CWF_Grouped_Atoms[gidx][Lidx]
                for Rn = 1:FNAN[atom]+1
        
                    jatom = natn[atom][Rn]
                    cell = ncn[atom][Rn]+1
                    NO1 = Total_NumOrbs[jatom]
                    kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
                    ex = cispi(-2*kRn)
                    Bnum = MP[jatom]
                    tmp1 = 0.0 + im*0.0
                    for jβ = 1:NO1
                        tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ+MinN-1])*OLPproj[gidx][p][Lidx][Rn][jβ]
                    end
                    Sum += tmp1*ex
                end
            end
    
            if weight_type == "fermi"
                weight = CWF_weight(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            else    # weight_type == "poly"
                weight = CWF_weight2(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
            end

            Amnk[1][ik][μ][pst] = weight * Sum
        end
    end
end


function Calc_Amnk_pol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup::CWF_Setup_MO, kpoints::KPoints)

    material = cwf_setup.material
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    Dis_Energy = cwf_setup.Dis_Energy
    weight_type = cwf_setup.weight_type
    Num_CWF_Grouped_Atoms = cwf_setup.Num_CWF_Grouped_Atoms
    Num_CWF_MOs_Group = cwf_setup.Num_CWF_MOs_Group
    CWF_Grouped_Atoms_EachNum = cwf_setup.CWF_Grouped_Atoms_EachNum
    CWF_Grouped_Atoms = cwf_setup.CWF_Grouped_Atoms

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        for gidx = 1:Num_CWF_Grouped_Atoms, p = 1:Num_CWF_MOs_Group[gidx]
            pst += 1
            Sum1 = ComplexF64(0.0, 0.0)
            Sum2 = ComplexF64(0.0, 0.0)
            for Lidx = 1:CWF_Grouped_Atoms_EachNum[gidx]
                atom = CWF_Grouped_Atoms[gidx][Lidx]
                for Rn = 1:FNAN[atom]+1
                    jatom = natn[atom][Rn]
                    cell = ncn[atom][Rn]+1
                    NO1 = Total_NumOrbs[jatom]
                    kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
                    ex = cispi(-2*kRn)
                    Bnum = MP[jatom]
                    tmp1 = 0.0 + im*0.0
                    tmp2 = 0.0 + im*0.0
                    for jβ = 1:NO1
                        tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ+MinN-1])*OLPproj[gidx][p][Lidx][Rn][jβ]
                        tmp2 += conj(Cnk[2][ik][Bnum+jβ,μ+MinN-1])*OLPproj[gidx][p][Lidx][Rn][jβ]
                    end
                    Sum1 += tmp1*ex
                    Sum2 += tmp2*ex
                end
            end
    
            if weight_type == "fermi"
                weight1 = CWF_weight(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
                weight2 = CWF_weight(Enk[2][ik][μ+MinN-1], ChemP, Dis_Energy)
            else    # weight_type == "poly"
                weight1 = CWF_weight2(Enk[1][ik][μ+MinN-1], ChemP, Dis_Energy)
                weight2 = CWF_weight2(Enk[2][ik][μ+MinN-1], ChemP, Dis_Energy)
            end

            Amnk[1][ik][μ][pst] = weight1 * Sum1
            Amnk[2][ik][μ][pst] = weight2 * Sum2
        end
    end
end


function Calc_Amnk_NonCol!(MinN, BANDNUM, OLPproj, Amnk, Enk, Cnk, cwf_setup::CWF_Setup_MO, kpoints::KPoints)

end