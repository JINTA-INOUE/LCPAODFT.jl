function Calc_Enk_Cnk_Dos(material::LCPAO_model, kpoints::KPoints, iemin, iemax)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol ∈ ("off", "on"), fsize, 2*fsize)
    spinsize = ifelse(SpinPol == "off", 1, 2)
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks
    ChemP = material.ChemP

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    neg = iemax-iemin+1

    Enk = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, spinsize)
    for spin = 1:spinsize
        Enk[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Enk[spin][ik] = Vector{Vector{Vector{Float64}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Enk[spin][ik][jk] = Vector{Vector{Float64}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Enk[spin][ik][jk][kk] = zeros(Float64, neg)
                end
            end
        end
    end
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, Nkpt)
        for ik = 1:Nkpt
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end


    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    tmpH = zeros(ComplexF64, fsize, fsize)
    Enk_tmp = zeros(Float64, Nfsize)


    if SpinPol ∈ ("off", "on")
        for spin = 1:spinsize
            kp = 0
            for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
                kp += 1
                HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp])

                Enk_tmp, Cnk[spin][kp] = eigen(Hermitian(H), Hermitian(S))
                for μ = 1:neg
                    Enk[spin][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP
                end
            end
        end
    elseif SpinPol == "nc"
        kp = 0
        for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            kp += 1
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk_tmp, Cnk[1][kp] = eigen(Hermitian(H), Hermitian(S))

            for μ = 1:neg
                Enk[1][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP
                Enk[2][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP      # copy for EVec
            end
        end
    end


    return Enk, Cnk
end


function Calc_Enk_Cnk_Dos(material::CWF_model, kpoints::KPoints, iemin, iemax)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Spindeg = ifelse(SpinPol=="off", 1, 2)
    Nwann = material.Ngsize
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    ChemP = material.ChemP
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = kpoints.AllNkpt
    kpts = kpoints.MPI_kpts
    neg = iemax-iemin+1


    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    Enk = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Spindeg)
    for spin = 1:Spindeg
        Enk[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Enk[spin][ik] = Vector{Vector{Vector{Float64}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Enk[spin][ik][jk] = Vector{Vector{Float64}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Enk[spin][ik][jk][kk] = zeros(Float64, neg)
                end
            end
        end
    end
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, Nkpt)
        for ik = 1:Nkpt
            Cnk[spin][ik] = zeros(ComplexF64, Nwann, Nwann)
        end
    end


    H = zeros(ComplexF64, Nwann, Nwann)
    for spin = 1:spinsize, k = 1:Nkpt
        k1, k2, k3 = kpts[k]
        fill!(H, 0.0)
        for cell = 1:NCell
            kRn = k1*cell_list_ijk[cell][1] + k2*cell_list_ijk[cell][2] + k3*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)

            for ist = 1:Nwann, jst = 1:Nwann
                H[ist,jst] += HmnR[jst,ist,cell,spin] * ex
            end
        end

        ik, jk, kk = kindex[k,:]
        Enk_tmp, Cnk[spin][k] = eigen(Hermitian(H))
        for μ = 1:neg
            Enk[spin][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP
        end
    end

    

    if SpinPol == "nc"
        for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3, μ = 1:neg
            Enk[2][ik][jk][kk][μ] = Enk[1][ik][jk][kk][μ]
        end
    end


    return Enk, Cnk
end