@timeit timer "Calc_MLWF_Amnk" function Calc_MLWF_Amnk(BANDNUM, Nk, Cnk, mlwf_setup::MLWF_Setup, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    WANNUM = mlwf_setup.WANNUM
    MPI_Nkpt = kpoints.MPI_Nkpt

    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, BANDNUM)
            for μ = 1:BANDNUM
                Amnk[spin][ik][μ] = zeros(ComplexF64, WANNUM)
            end
        end
    end

    if SpinPol ∈ ("off", "on")
        Calc_MLWF_Amnk_Col!(BANDNUM, Nk, Amnk, Cnk, mlwf_setup, kpoints)
    else
        Calc_MLWF_Amnk_NonCol!(BANDNUM, Nk, Amnk, Cnk, mlwf_setup, kpoints)
    end


    return Amnk
end


function Calc_MLWF_Amnk_Col!(BANDNUM, Nk, Amnk, Cnk, mlwf_setup::MLWF_Setup, kpoints::KPoints)

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    OLP = material.OLP
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Guide_index = mlwf_setup.Guide_index
    GNatom = length(Guide_index)

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        _Cnk = Cnk[spin][ik]
        _Amnk = Amnk[spin][ik]
        for μ = 1:BANDNUM
            ka, kb, kc = MPI_kpts[ik]
            pst = 0
            midx = μ + Nk[spin,ik,1]
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
                        tmp1 += conj(_Cnk[Bnum+jβ,midx])*OLP[patom][Rn][proj][jβ]
                    end
                    Sum += tmp1*ex
                end
                _Amnk[μ][pst] = Sum
            end
        end
    end
end


function Calc_MLWF_Amnk_NonCol!(BANDNUM, Nk, Amnk, Cnk, mlwf_setup::MLWF_Setup, kpoints::KPoints)

    material = mlwf_setup.material
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    OLP = material.OLP
    Guide_index = mlwf_setup.Guide_index
    GNatom = length(Guide_index)
    WANNUM_half = 0
    for atom = 1:GNatom
        WANNUM_half += length(Guide_index[atom])
    end
    
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    for ik = 1:MPI_Nkpt, μ = 1:BANDNUM
        ka, kb, kc = MPI_kpts[ik]
        pst = 0
        midx = μ + Nk[1,ik,1]
        _Cnk = Cnk[1][ik]
        _Amnk = Amnk[1][ik]
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
                    tmp1 += conj(_Cnk[Bnum+jβ,midx])*OLP[patom][Rn][proj][jβ]
                    tmp2 += conj(_Cnk[fsize+Bnum+jβ,midx])*OLP[patom][Rn][proj][jβ]
                end
                Sum1 += tmp1*ex
                Sum2 += tmp2*ex
            end
        
            _Amnk[μ][pst] = Sum1
            _Amnk[μ][WANNUM_half+pst] =  Sum2
        end
    end
end