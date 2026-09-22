function UnfoldBand_kpath_LCPAO(filepath::String, filename::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}, Nk::Integer)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    # This is the pre-shared-memory implementation: every MPI process loads
    # and owns a complete LCPAO model independently.
    material = Load_LCPAODFT_model(filepath)
    myrank == 0 && Print_LCPAO_model(filepath, material)

    Natom = material.Natom
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "on", 2, 1)
    EVec_spinsize = ifelse(SpinPol == "off", 1, 2)
    Recvecs = material.Recvecs
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol == "nc", 2*fsize, fsize)
    Hks = material.Hks
    OLP = material.OLP
    iHks = material.iHks
    ChemP = material.ChemP

    Nkpath, kpath_Nk, kpath_start, kpath_end, all_kpts = _setup_band_kpath(kpath, Nk)
    MPI_krange, MPI_kpts = _split_band_kpoints(all_kpts, comm)
    MPI_Nkpt = length(MPI_kpts)

    # k point is the last dimension so that contiguous rank-owned ranges can be
    # gathered directly into the same global ordering on rank 0.
    local_Enk = zeros(Float64, Nfsize, spinsize, MPI_Nkpt)
    local_EVec = zeros(Float32, Nfsize, fsize, EVec_spinsize, MPI_Nkpt)


    if SpinPol ∈ ("off", "on")
        SD = zeros(Float32, fsize)
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for spin = 1:spinsize, ik = 1:MPI_Nkpt
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            local_Enk[:,spin,ik], Cnk = eigen(Hermitian(H), Hermitian(S))
            Calc_EVec_Col!(spin, ik, Natom, Total_NumOrbs, MP, S, SD, Cnk, local_EVec)
        end
    elseif SpinPol == "nc"
        Angle_spin = Get_Angle_spin(material)
        SD = zeros(Float32, 2*fsize)
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)

        for ik = 1:MPI_Nkpt
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
            fill!(S, 0.0)
            @views S[1:fsize,1:fsize] .= tmpH
            @views S[fsize+1:end,fsize+1:end] .= tmpH
            local_Enk[:,1,ik], Cnk = eigen(Hermitian(H), Hermitian(S))
            Calc_EVec_NonCol!(ik, Natom, Total_NumOrbs, MP, Angle_spin, tmpH, SD, Cnk, local_EVec)
        end
    else
        error("unsupported spin polarization: $SpinPol")
    end
    

    Enk = _gather_band_energies(local_Enk, MPI_krange, length(all_kpts), comm)
    EVec = _gather_evec(local_EVec, MPI_krange, length(all_kpts), comm)
    if myrank == 0
        for spin = 1:spinsize
            Write_UnfoldBANDDAT(filename, SpinPol, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, EVec, ChemP, Recvecs)
        end
        Write_GNUUNFOLDBAND(filename, spinsize, Nkpath, kpath, kname, Recvecs)
    end
    MPI.Barrier(comm)
end


function UnfoldBand_kpath_CWF(filepath::String, filename::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}, Nk::Integer)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    # This is the pre-shared-memory implementation: every MPI process loads
    # and owns a complete CWF model independently.
    material = Load_CWF_model(filepath)
    myrank == 0 && Print_CWF_model(filepath, material)

    SpinPol = material.SpinPol
    spinsize = material.spinsize
    EVec_spinsize = ifelse(SpinPol == "off", 1, 2)
    Recvecs = material.Recvecs
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    Nwann = material.Ngsize
    HmnR = material.HmnR
    ChemP = material.ChemP
    Nwann_EVec = ifelse(SpinPol=="nc", div(Nwann,2), Nwann)

    Nkpath, kpath_Nk, kpath_start, kpath_end, all_kpts = _setup_band_kpath(kpath, Nk)
    MPI_krange, MPI_kpts = _split_band_kpoints(all_kpts, comm)
    MPI_Nkpt = length(MPI_kpts)

    local_Enk = zeros(Float64, Nwann, spinsize, MPI_Nkpt)
    local_EVec = zeros(Float32, Nwann, Nwann_EVec, EVec_spinsize, MPI_Nkpt)
    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        fill!(H, 0.0)
        k1, k2, k3 = MPI_kpts[ik]

        for cell = 1:NCell
            kRn = k1*cell_list_ijk[cell][1] + k2*cell_list_ijk[cell][2] + k3*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)
            @inbounds for ist = 1:Nwann, jst = 1:Nwann
                H[jst,ist] += HmnR[jst,ist,cell,spin]*ex
            end
        end

        local_Enk[:,spin,ik], Cnk = eigen(Hermitian(H))
        Calc_EVec_CWF!(SpinPol, spin, ik, Nwann, Cnk, local_EVec)
    end


    Enk = _gather_band_energies(local_Enk, MPI_krange, length(all_kpts), comm)
    EVec = _gather_evec(local_EVec, MPI_krange, length(all_kpts), comm)
    if myrank == 0
        for spin = 1:spinsize
            Write_UnfoldBANDDAT(filename, SpinPol, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nwann, Enk, EVec, ChemP, Recvecs)
        end
        Write_GNUUNFOLDBAND(filename, spinsize, Nkpath, kpath, kname, Recvecs)
    end

    MPI.Barrier(comm)
end


function UnfoldBand_kpath(filepath::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}; seedname=splitext(basename(filepath))[1], Nk::Integer=50)

    Threads.nthreads() == 1 || error(
        "MPI-flat Band_kpath requires exactly one Julia thread per MPI process; " *
        "start Julia with --threads=1",
    )
    provided_thread_level = MPI.Init(; threadlevel=:single)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    BLAS.set_num_threads(1)

    length(kname) == length(kpath) || error("kname and kpath must contain the same number of entries")

    if myrank == 0
        println("<Band MPI-flat configuration>")
        println("\t$nprocs MPI processes x 1 Julia thread")
        println("\t$(BLAS.get_num_threads()) BLAS thread per process")
        println("\tMPI thread level: $provided_thread_level")
    end

    model = select_model(filepath)
    if model == 1
        UnfoldBand_kpath_LCPAO(filepath, seedname, kpath, kname, Nk)
    elseif model == 2
        UnfoldBand_kpath_CWF(filepath, seedname, kpath, kname, Nk)
    else
        error("unsupported model file: $filepath")
    end
end
