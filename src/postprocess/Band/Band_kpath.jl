function Band_kpath(filepath::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}; Nk=50)

    filename, _ = splitext(basename(filepath))

    model = Load_JLD2(filepath)

    Natom = model.Natom
    Nspin = model.Nspin
    SpinPol = model.SpinPol
    Recvecs = model.Recvecs
    Total_NumOrbs = model.Total_NumOrbs
    MP = model.MP
    FNAN = model.FNAN
    natn = model.natn
    ncn = model.ncn
    atv_ijk = model.atv_ijk
    fsize = sum(Total_NumOrbs)

    Hks = model.Hks
    OLP = model.OLP
    iHks = model.iHks
    ChemP = model.ChemP
    
    spinsize = ifelse(SpinPol ∈ ("off","nc"), 1, 2)
    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
    end
    
    

    Nkpath = length(kpath)-1
    kpath_Nk = ones(Int64, Nkpath)*Nk


    kpath_start = Vector{Vector{Float64}}(undef, Nkpath)
    kpath_end = Vector{Vector{Float64}}(undef, Nkpath)
    for ik = 1:Nkpath
        kpath_start[ik] = kpath[ik]
        kpath_end[ik] = kpath[ik+1]
    end

    Enk = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nkpath)
    for ik = 1:Nkpath
        Enk[ik] = Vector{Vector{Vector{Float64}}}(undef, kpath_Nk[ik])
        for ipath = 1:kpath_Nk[ik]
            Enk[ik][ipath] = Vector{Vector{Float64}}(undef, spinsize)
            for spin = 1:spinsize
                Enk[ik][ipath][spin] = zeros(Float64, Nfsize)
            end
        end
    end




    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

            k1 = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
            k2 = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
            k3 = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

            for spin = 1:Nspin
                HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, [k1,k2,k3])
                HS_matrix!(H, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, [k1,k2,k3])

                Enk[ik][ipath][spin] = eigvals(Hermitian(H), Hermitian(S))
            end
        end
    else SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

            k1 = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
            k2 = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
            k3 = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, [k1,k2,k3])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, [k1,k2,k3])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH

            Enk[ik][ipath][1] = eigvals(Hermitian(H), Hermitian(S))
        end
    end

    Write_BANDDAT(spinsize, filename, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, ChemP, Recvecs)
    Write_GNUBAND(spinsize, filename, Nkpath, kpath, kname, Recvecs)
end