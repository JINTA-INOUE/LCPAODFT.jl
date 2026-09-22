@timeit timer "Calc_Enk_Cnk" function Calc_Enk_Cnk(material::LCPAO_model, kpoints::KPoints, type::Integer)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    MPI_Nkpt = kpoints.MPI_Nkpt
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)

    if type == 1
        Enk = zeros(Float64, Nfsize, MPI_Nkpt, spinsize)
    elseif type == 2
        Enk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
        for spin = 1:spinsize
            Enk[spin] = Vector{Vector{Float64}}(undef, MPI_Nkpt)
            for ik = 1:MPI_Nkpt
                Enk[spin][ik] = zeros(Float64, Nfsize)
            end
        end
    elseif type == 3
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

    if SpinPol == "off"
        S = zeros(ComplexF64, fsize, fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H = Cnk[1][ik]
            HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[:,ik,1] .= decomposition.values
        end
    elseif SpinPol == "on"
        S_base = zeros(ComplexF64, fsize, fsize)
        S_work = zeros(ComplexF64, fsize, fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H_up = Cnk[1][ik]
            HS_matrix!(S_base, H_up, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            copyto!(S_work, S_base)
            decomposition_up = eigen!(Hermitian(H_up), Hermitian(S_work))
            @views Enk[:,ik,1] .= decomposition_up.values

            H_dn = Cnk[2][ik]
            HS_matrix!(H_dn, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            copyto!(S_work, S_base)
            decomposition_dn = eigen!(Hermitian(H_dn), Hermitian(S_work))
            @views Enk[:,ik,2] .= decomposition_dn.values
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H = Cnk[1][ik]
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            fill!(S, 0.0)
            @. @views S[1:fsize, 1:fsize] = tmpH
            @. @views S[fsize+1:end, fsize+1:end] = tmpH
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[:,ik,1] .= decomposition.values
        end
    end
end


function Calc_Enk_Cnk!(
    material::LCPAO_model, 
    kpoints::KPoints, 
    Enk::Vector{Vector{Vector{Float64}}}, 
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

    if SpinPol == "off"
        S = zeros(ComplexF64, fsize, fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H = Cnk[1][ik]
            HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[1][ik] .= decomposition.values
        end
    elseif SpinPol == "on"
        S_base = zeros(ComplexF64, fsize, fsize)
        S_work = zeros(ComplexF64, fsize, fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H_up = Cnk[1][ik]
            HS_matrix!(S_base, H_up, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            copyto!(S_work, S_base)
            decomposition_up = eigen!(Hermitian(H_up), Hermitian(S_work))
            @views Enk[1][ik] .= decomposition_up.values

            H_dn = Cnk[2][ik]
            HS_matrix!(H_dn, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            copyto!(S_work, S_base)
            decomposition_dn = eigen!(Hermitian(H_dn), Hermitian(S_work))
            @views Enk[2][ik] .= decomposition_dn.values
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:MPI_Nkpt
            H = Cnk[1][ik]
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            fill!(S, 0.0)
            @. @views S[1:fsize, 1:fsize] = tmpH
            @. @views S[fsize+1:end, fsize+1:end] = tmpH
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[1][ik] .= decomposition.values
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
    

    if SpinPol == "off"
        S = zeros(ComplexF64, fsize, fsize)
        k = 0
        @inbounds for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            k += 1
            H = Cnk[spin][ik]
            HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[spin][ik] .= decomposition.values
        end
    elseif SpinPol == "on"
        S_base = zeros(ComplexF64, fsize, fsize)
        S_work = zeros(ComplexF64, fsize, fsize)
        k = 0
        @inbounds for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            k += 1
            H_up = Cnk[1][k]
            HS_matrix!(S_base, H_up, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            copyto!(S_work, S_base)
            decomposition_up = eigen!(Hermitian(H_up), Hermitian(S_work))
            @views Enk[1][k] .= decomposition_up.values

            H_dn = Cnk[2][k]
            HS_matrix!(H_dn, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            copyto!(S_work, S_base)
            decomposition_dn = eigen!(Hermitian(H_dn), Hermitian(S_work))
            @views Enk[2][k] .= decomposition_dn.values
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        k = 0
        @inbounds for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            k += 1
            H = Cnk[1][k]
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[k])
            fill!(S, 0.0)
            @. @views S[1:fsize, 1:fsize] = tmpH
            @. @views S[fsize+1:end, fsize+1:end] = tmpH
            decomposition = eigen!(Hermitian(H), Hermitian(S))
            @views Enk[1][ik] .= decomposition.values
        end
    end
end


function Calc_Enk_Cnk(material::CWF_model, kpoints::KPoints, type::Integer)

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    MPI_Nkpt = kpoints.MPI_Nkpt
    spinsize = ifelse(SpinPol=="on", 2, 1)

    if type == 1
        Enk = zeros(Float64, spinsize, Nwann, MPI_Nkpt)
    elseif type == 2
        Enk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
        for spin = 1:spinsize
            Enk[spin] = Vector{Vector{Float64}}(undef, MPI_Nkpt)
            for ik = 1:MPI_Nkpt
                Enk[spin][ik] = zeros(Float64, Nwann)
            end
        end
    elseif type == 3
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
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    spinsize = ifelse(SpinPol=="on", 2, 1)


    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        k1, k2, k3 = MPI_kpts[ik]
        fill!(H, 0.0)
        for cell = 1:NCell
            kRn = k1*cell_list_ijk[cell][1] + k2*cell_list_ijk[cell][2] + k3*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)

            for ist = 1:Nwann, jst = 1:Nwann
                H[ist,jst] += HmnR[jst,ist,cell,spin] * ex
            end
        end

        Enk[:,ik,1], Cnk[spin][ik] = eigen(Hermitian(H))
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
    spinsize = ifelse(SpinPol=="on", 2, 1)

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