function Calc_Enk_Cnk(material::LCPAO_model, kpoints::KPoints, type::Integer)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    MPI_Nkpt = kpoints.MPI_Nkpt

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else
        spinsize = 2
    end

    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    else
        Nfsize = 2*fsize
    end

    if type == 1
        Enk = zeros(Float64, spinsize, Nfsize, MPI_Nkpt)
    elseif type == 2
        Enk = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, spinsize)
        for spin = 1:spinsize
            Enk[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh1)
            for ik = 1:kmesh1
                Enk[spin][ik] = Vector{Vector{Vector{Float64}}}(undef, kmesh2)
                for jk = 1:kmesh2
                    Enk[spin][ik][jk] = Vector{Vector{Float64}}(undef, kmesh3)
                    for kk = 1:kmesh3
                        Enk[spin][ik][jk][kk] = zeros(Float64, Nfsize)
                    end
                end
            end
        end
    else
        error("please check type.")
    end

    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end

    Calc_Enk_Cnk!(material, kpoints, Enk, Cnk)


    return Enk, Cnk
end


function Calc_Enk_Cnk!(
    material::LCPAO_model, 
    kpoints::KPoints, 
    Enk::Array{Float64,3}, 
    Cnk::Vector{Vector{Matrix{ComplexF64}}})

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    Hks = material.Hks
    iHks = material.iHks
    OLP = material.OLP

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else SpinPol == "on"
        spinsize = 2
    end


    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        @inbounds for spin = 1:spinsize, ik = 1:MPI_Nkpt
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            Enk[spin,:,ik], Cnk[spin][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:MPI_Nkpt
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk[1,:,ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    end
end


function Calc_Enk_Cnk!(
    material::LCPAO_model, 
    kpoints::KPoints, 
    Enk::Vector{Vector{Vector{Vector{Vector{Float64}}}}}, 
    Cnk::Vector{Vector{Matrix{ComplexF64}}})

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    Hks = material.Hks
    iHks = material.iHks
    OLP = material.OLP

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    MPI_kpts = kpoints.MPI_kpts
    

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else
        spinsize = 2
    end


    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        @inbounds for spin = 1:spinsize
            k = 0
            for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
                k += 1
                HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
                Enk[spin][ik][jk][kk], Cnk[spin][k] = eigen(Hermitian(H), Hermitian(S))
            end
        end
    else
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        k = 0
        @inbounds for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            k += 1
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk[1][ik][jk][kk], Cnk[1][k] = eigen(Hermitian(H), Hermitian(S))
        end
    end
end


function Calc_Enk_Cnk(material::CWF_model, kpoints::KPoints, type::Integer)

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    MPI_Nkpt = kpoints.MPI_Nkpt

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else SpinPol == "on"
        spinsize = 2
    end


    if type == 1
        Enk = zeros(Float64, spinsize, Nwann, MPI_Nkpt)
    elseif type == 2
        Enk = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, spinsize)
        for spin = 1:spinsize
            Enk[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh1)
            for ik = 1:kmesh1
                Enk[spin][ik] = Vector{Vector{Vector{Float64}}}(undef, kmesh2)
                for jk = 1:kmesh2
                    Enk[spin][ik][jk] = Vector{Vector{Float64}}(undef, kmesh3)
                    for kk = 1:kmesh3
                        Enk[spin][ik][jk][kk] = zeros(Float64, Nwann)
                    end
                end
            end
        end
    else
        error("please check type.")
    end

    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Cnk[spin][ik] = zeros(ComplexF64, Nwann, Nwann)
        end
    end

    Calc_Enk_Cnk!(material, kpoints, Enk, Cnk)


    return Enk, Cnk
end


function Calc_Enk_Cnk!(
    material::CWF_model,
    kpoints::KPoints, 
    Enk::Array{Float64,3}, 
    Cnk::Vector{Vector{Matrix{ComplexF64}}})

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    Nkpt = kpoints.MPI_Nkpt
    kpts = kpoints.MPI_kpts

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else
        spinsize = 2
    end
    

    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, ik = 1:Nkpt
        for cell = 1:NCell
            kRn = kpts[ik][1]*cell_list_ijk[cell][1] + kpts[ik][2]*cell_list_ijk[cell][2] + kpts[ik][3]*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)

            for ist = 1:Nwann, jst = 1:Nwann
                H[ist,jst] += HmnR[jst,ist,cell,spin] * ex
            end
        end

        Enk[1,:,ik], Cnk[spin][k] = eigen(Hermitian(H))
    end
end


function Calc_Enk_Cnk!(
    material::CWF_model,
    kpoints::KPoints, 
    Enk::Vector{Vector{Vector{Vector{Vector{Float64}}}}}, 
    Cnk::Vector{Vector{Matrix{ComplexF64}}})

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = kpoints.MPI_Nkpt
    kpts = kpoints.MPI_kpts

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else
        spinsize = 2
    end
    

    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
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
        Enk[spin][ik][jk][kk], Cnk[spin][k] = eigen(Hermitian(H))
    end
end