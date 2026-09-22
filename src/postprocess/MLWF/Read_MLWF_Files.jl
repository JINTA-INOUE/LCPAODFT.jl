function Read_Mmnkb(filepath::String)

    data = open(filepath, "r")
    f = readline(data)
    f = readline(data)
    BANDNUM = parse(Int64, split(f)[1])
    Nkpt = parse(Int64, split(f)[2])
    tot_bvector = parse(Int64, split(f)[3])

    Mmnkb = Vector{Vector{Matrix{ComplexF64}}}(undef, Nkpt)
    for ik = 1:Nkpt
        Mmnkb[ik] = Vector{Matrix{ComplexF64}}(undef, tot_bvector)
        for ib = 1:tot_bvector
            Mmnkb[ik][ib] = zeros(ComplexF64, BANDNUM, BANDNUM)
        end
    end

    
    for ik = 1:Nkpt, ib = 1:tot_bvector
        f = readline(data)
        for m = 1:BANDNUM, n = 1:BANDNUM
            f = readline(data)
            re = parse(Float64, split(f)[1])
            im = parse(Float64, split(f)[2])
            Mmnkb[ik][ib][n,m] = ComplexF64(re, im)
        end
    end
    close(data)
    

    return Mmnkb
end


function Read_Amnk(filepath::String)

    data = open(filepath, "r")
    f = readline(data)
    f = readline(data)
    BANDNUM = parse(Int64, split(f)[1])
    Nkpt = parse(Int64, split(f)[2])
    WANNUM = parse(Int64, split(f)[3])
    spinsize = parse(Int64, split(f)[4])


    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, Nkpt)
        for ik = 1:Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, BANDNUM)
            for μ = 1:BANDNUM
                Amnk[spin][ik][μ] = Vector{ComplexF64}(undef, WANNUM)
            end
        end
    end

    for spin = 1:spinsize, ik = 1:Nkpt
        for n = 1:WANNUM, μ = 1:BANDNUM
            f = readline(data)
            re = parse(Float64, split(f)[4])
            im = parse(Float64, split(f)[5])
            Amnk[spin][ik][μ][n] = ComplexF64(re,im)
        end
    end
    close(data)
    

    return Amnk
end


function Read_eigen(filepath::String, spinsize, Nkpt, BANDNUM, ChemP)

    eV2Hartree = 27.2113845

    eigen = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        eigen[spin] = Vector{Vector{Float64}}(undef, Nkpt)
        for ik = 1:Nkpt
            eigen[spin][ik] = Vector{Float64}(undef, BANDNUM)
        end
    end

    data = open(filepath, "r")
    for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:BANDNUM
        f = readline(data)
        eigen[spin][ik][μ] = parse(Float64, split(f)[3])/eV2Hartree+ChemP
    end
    close(data)


    return eigen
end


function Read_chk(filename::AbstractString)

    data = FortranFile(filename)

    # strip and read line
    header_len = 33
    header = trimstring(read(data, FString{header_len}))

    # gfortran default integer is 4 bytes
    Tint = Int32
    BANDNUM = read(data, Tint)
    n_exclude_bands = read(data, Tint)
    exclude_bands = zeros(Int, n_exclude_bands)
    exclude_bands .= read(data, (Tint, n_exclude_bands))
    n_exclude_bands == 0 && read(data)
    Latvecs = read(data, (Float64, 3, 3))
    Recvecs = read(data, (Float64, 3, 3))
    Nkpt = read(data, Tint)
    kmesh = read(data, (Tint, 3))
    kpoints = zeros(Float64, 3*Nkpt)
    read(data, kpoints)
    kpts = Vector{Vector{Float64}}(undef, Nkpt)
    for ik = 1:Nkpt
        kpts[ik] = zeros(Float64, 3)
        for xyz = 1:3
            kpts[ik][xyz] = kpoints[3*(ik-1)+xyz]
        end
    end
    n_bvecs = read(data, Tint)
    WANNUM = read(data, Tint)
    checkpoint = trimstring(read(data, FString{20}))
    have_disentangled = parse_bool(read(data, Tint))

    if have_disentangled
        # omega_invariant
        ΩI = read(data, Float64)

        tmp = parse_bool.(read(data, (Tint, BANDNUM, Nkpt)))
        dis_bands = [tmp[:, i] for i in 1:Nkpt]
        n_dis = zeros(Int, Nkpt)
        n_dis .= read(data, (Tint, Nkpt))
        for ik in 1:Nkpt
            n_dis[ik] == count(dis_bands[ik]) || error("Inconsistent number of disentangled bands")
        end

        # u_matrix_opt
        U_tmp = zeros(ComplexF64, BANDNUM, WANNUM, Nkpt)
        read(data, U_tmp)
        Udis = [U_tmp[:, :, ik] for ik = 1:Nkpt]
    else
        ΩI = -1.0
        dis_bands = BitVector[]
        n_dis = Int[]
        Udis = Matrix{ComplexF64}[]
    end

    Uml = zeros(ComplexF64, WANNUM, WANNUM, Nkpt)
    Mmnkb = zeros(ComplexF64, WANNUM, WANNUM, n_bvecs, Nkpt)
    WannierCenter = zeros(Float64, 3, WANNUM)
    Omega_tot = zeros(Float64, WANNUM)
    read(data, Uml)
    read(data, Mmnkb)
    read(data, WannierCenter)
    read(data, Omega_tot)
    close(data)


    return w90_Chk(
        header, BANDNUM, n_exclude_bands, exclude_bands, Latvecs, Recvecs,
        Nkpt, kmesh, kpts, n_bvecs, WANNUM, checkpoint, have_disentangled,
        ΩI, dis_bands, n_dis, Udis,
        [Uml[:, :, ik] for ik = 1:Nkpt],
        [[Mmnkb[:, :, ib, ik] for ib = 1:n_bvecs] for ik = 1:Nkpt],
        WannierCenter, Omega_tot)
end