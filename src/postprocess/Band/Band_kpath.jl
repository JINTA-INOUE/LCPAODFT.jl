function Band_kpath_LCPAO(filepath::String, filename::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}, Nk)

    model = Load_LCPAODFT_model(filepath)

    Natom = model.Natom
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
    
    Enk = Vector{Vector{Vector{Float64}}}(undef, Nkpath)
    for ik = 1:Nkpath
        Enk[ik] = Vector{Vector{Float64}}(undef, kpath_Nk[ik])
        for ipath = 1:kpath_Nk[ik]
            Enk[ik][ipath] = zeros(Float64, Nfsize)
        end
    end


    kpts = zeros(Float64, 3)

    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for spin = 1:spinsize
            for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

                kpts[1] = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
                kpts[2] = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
                kpts[3] = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

                HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)

                Enk[ik][ipath] = eigvals(Hermitian(H), Hermitian(S))
            end

            Write_BANDDAT(filename, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, ChemP, Recvecs)
        end
    else SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

            kpts[1] = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
            kpts[2] = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
            kpts[3] = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH

            Enk[ik][ipath] = eigvals(Hermitian(H), Hermitian(S))
        end

        Write_BANDDAT(filename, 1, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, ChemP, Recvecs)
    end

    Write_GNUBAND(filename, spinsize, Nkpath, kpath, kname, Recvecs)
end


function Band_kpath_CWF(filepath::String, filename::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}, Nk)
    
    cwf_model = Load_CWF_model(filepath)
    spinsize = cwf_model.spinsize
    Recvecs = cwf_model.Recvecs
    NCell = cwf_model.NCell
    cell_list_ijk = cwf_model.cell_list_ijk
    Nwann = cwf_model.Ngsize
    HmnR = cwf_model.HmnR
    ChemP = cwf_model.ChemP

    Nkpath = length(kpath)-1
    kpath_Nk = ones(Int64, Nkpath)*Nk


    kpath_start = Vector{Vector{Float64}}(undef, Nkpath)
    kpath_end = Vector{Vector{Float64}}(undef, Nkpath)
    for ik = 1:Nkpath
        kpath_start[ik] = kpath[ik]
        kpath_end[ik] = kpath[ik+1]
    end

    
    Enk = Vector{Vector{Vector{Float64}}}(undef, Nkpath)
    for ik = 1:Nkpath
        Enk[ik] = Vector{Vector{Float64}}(undef, kpath_Nk[ik])
        for ipath = 1:kpath_Nk[ik]
            Enk[ik][ipath] = zeros(Float64, Nwann)
        end
    end


    H = zeros(ComplexF64, Nwann, Nwann)
    for spin = 1:spinsize, ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

        fill!(H, 0.0)

        k1 = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
        k2 = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
        k3 = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

        for cell = 1:NCell
            kRn = k1*cell_list_ijk[cell][1] + k2*cell_list_ijk[cell][2] + k3*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)

            for ist = 1:Nwann, jst = 1:Nwann
                H[jst,ist] += HmnR[jst,ist,cell,spin] * ex
            end
        end

        Enk[ik][ipath] = eigvals(Hermitian(H))
        
        Write_BANDDAT(filename, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nwann, Enk, ChemP, Recvecs)
    end


    Write_GNUBAND(filename, spinsize, Nkpath, kpath, kname, Recvecs)
end


function Band_kpath(filepath::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}; seedname=splitext(basename(filepath))[1], Nk=50)

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    if nprocs > 1
        error("please run serial.")
    end

    base_filepath = basename(filepath)
    file = split(base_filepath, ".")
    if file[end] ≠ "jld2"
        error("please check filepath.")
    end
    if length(file) == 2
        ext = "jld2"
    else
        ext = file[end-1]*"."*file[end]
    end
    if ext == "jld2"
        Band_kpath_LCPAO(filepath, seedname, kpath, kname, Nk)
    elseif ext == "CWF.jld2"
        Band_kpath_CWF(filepath, seedname, kpath, kname, Nk)
    else
        error("please check filepath.")
    end
end