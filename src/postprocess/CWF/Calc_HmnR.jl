@timeit timer "Calc_HmnR!" function Calc_HmnR!(HmnR, spinsize, MinN, MaxN, NCell, Ngsize, cell_list_ijk, Umnk, Enk, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    Nkpt = kpoints.Nkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    BANDNUM = MaxN - MinN + 1

    Hk = zeros(ComplexF64, Ngsize, Ngsize, MPI_Nkpt)
    phase = zeros(ComplexF64, MPI_Nkpt)

    for spin = 1:spinsize
        for pst = 1:Ngsize, qst = 1:Ngsize, ik = 1:MPI_Nkpt
            tmp = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:BANDNUM
                tmp += Enk[spin][ik][μ+MinN-1] * conj(Umnk[spin][ik][μ,pst]) * Umnk[spin][ik][μ,qst]
            end
            Hk[qst,pst,ik] = tmp
        end

        for cell = 1:NCell
            @inbounds for ik = 1:MPI_Nkpt
                kRn = dot(MPI_kpts[ik], cell_list_ijk[cell])
                phase[ik] = cispi(-2*kRn)
            end
            for pst = 1:Ngsize, qst = 1:Ngsize
                Sum = ComplexF64(0.0, 0.0)
                @inbounds for ik = 1:MPI_Nkpt
                    Sum += Hk[qst,pst,ik]*phase[ik]
                end
                HmnR[qst,pst,cell,spin] = Sum/Nkpt
            end
        end
    end


    MPI.Allreduce!(HmnR, MPI.SUM, comm)
end


const HMNR_BLOCK_TARGET_BYTES = 64 * 1024^2


function _HmnR_block_cells(Ngsize, NCell, MPI_Nkpt=0; target_bytes=HMNR_BLOCK_TARGET_BYTES)

    # Each cell in a block needs one Ngsize x Ngsize result matrix and one
    # phase per rank-local k point.
    bytes_per_cell = sizeof(ComplexF64) * (Int(Ngsize) * Int(Ngsize) + Int(MPI_Nkpt))
    return clamp(div(target_bytes, max(bytes_per_cell, 1)), 1, Int(NCell))
end


function _HmnR_fft_block_cells(Ngsize, NCell, Nkpt; target_bytes=HMNR_BLOCK_TARGET_BYTES)

    # Rank 0 keeps the complete H(k) grid for the in-place FFT.  Use the
    # remaining workspace for the contiguous output block.
    matrix_bytes = sizeof(ComplexF64) * Int128(Ngsize) * Int128(Ngsize)
    fixed_bytes = matrix_bytes * Int128(Nkpt)
    available_bytes = max(Int128(target_bytes) - fixed_bytes, matrix_bytes)
    return clamp(Int(div(available_bytes, matrix_bytes)), 1, Int(NCell))
end


function _HmnR_memory_GiB(number_of_elements)
    return Float64(number_of_elements * sizeof(ComplexF64)) / 1024^3
end


function _Write_HmnR_block!(HmnR_dataset, HmnR_block, first_cell, last_cell, spin, NCell, Ngsize)

    cells_in_block = last_cell - first_cell + 1
    block_view = @view HmnR_block[:, :, 1:cells_in_block]

    # JLD2.ArrayDataset currently implements a range assignment through its
    # scalar setindex! fallback.  HmnR is an uncompressed contiguous
    # ComplexF64 dataset, so write the already-contiguous cell block in one IO
    # operation.  Keep the indexed fallback for other JLD2 implementations.
    if hasproperty(HmnR_dataset, :f) && hasproperty(HmnR_dataset, :data_address)
        matrix_elements = Int64(Ngsize) * Int64(Ngsize)
        element_offset = ((Int64(spin) - 1) * Int64(NCell) + Int64(first_cell) - 1) * matrix_elements
        byte_offset = element_offset * sizeof(ComplexF64)
        expected_bytes = length(block_view) * sizeof(ComplexF64)

        io = HmnR_dataset.f.io
        seek(io, HmnR_dataset.data_address + byte_offset)
        written_bytes = write(io, block_view)
        written_bytes == expected_bytes || error("incomplete HmnR block write: wrote $written_bytes of $expected_bytes bytes",)
    else
        HmnR_dataset[:, :, first_cell:last_cell, spin] = block_view
    end

    return nothing
end


function _HmnR_fft_compatible(kpoints::KPoints)
    return Int(kpoints.Nkpt) == prod(Int.(kpoints.kmesh)) && lowercase(kpoints.KP_flag) == "gcenter" && iszero(kpoints.Shift_K_Point)
end


@timeit timer "Calc_HmnR_stream!" function Calc_HmnR_stream!(
    HmnR_dataset, spinsize,
    MinN, MaxN, NCell, Ngsize,
    cell_list_ijk, Umnk, Enk, kpoints::KPoints; 
    block_cells=nothing,)

    if _HmnR_fft_compatible(kpoints)
        resolved_block_cells = isnothing(block_cells) ? _HmnR_fft_block_cells(Ngsize, NCell, kpoints.Nkpt) : clamp(Int(block_cells), 1, Int(NCell))
        return _Calc_HmnR_stream_fft!(
            HmnR_dataset, spinsize,
            MinN, MaxN,
            NCell, Ngsize,
            cell_list_ijk, Umnk, Enk, kpoints,
            resolved_block_cells)
    end

    resolved_block_cells = isnothing(block_cells) ? _HmnR_block_cells(Ngsize, NCell, kpoints.MPI_Nkpt) : clamp(Int(block_cells), 1, Int(NCell))
    return _Calc_HmnR_stream_direct!(
        HmnR_dataset,
        spinsize,
        MinN,
        MaxN,
        NCell,
        Ngsize,
        cell_list_ijk,
        Umnk,
        Enk,
        kpoints,
        resolved_block_cells,
    )
end


function _Calc_HmnR_stream_fft!(
    HmnR_dataset,
    spinsize,
    MinN,
    MaxN,
    NCell,
    Ngsize,
    cell_list_ijk,
    Umnk,
    Enk,
    kpoints::KPoints,
    block_cells,
)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = Int(kpoints.Nkpt)
    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    BANDNUM = MaxN - MinN + 1
    kmesh1, kmesh2, kmesh3 = Int.(kpoints.kmesh)
    matrix_elements = Int(Ngsize) * Int(Ngsize)

    # KPoints distributes contiguous ranges of the global k-point ordering.
    # On rank 0, its local range is already the first part of Hk_global and
    # can therefore participate in Gatherv! in place.
    Hk_global = myrank == 0 ? zeros(ComplexF64, Ngsize, Ngsize, Nkpt) : nothing
    Hk_local = myrank == 0 ? @view(Hk_global[:, :, 1:MPI_Nkpt]) : zeros(ComplexF64, Ngsize, Ngsize, MPI_Nkpt)
    weighted_conj_U = zeros(ComplexF64, BANDNUM, Ngsize)
    receive_counts = [matrix_elements * length(krange) for krange in kpoints.MPI_krange]

    Hk_grid = myrank == 0 ? reshape(Hk_global, Ngsize, Ngsize, kmesh3, kmesh2, kmesh1) : nothing
    fft_plan = myrank == 0 ? plan_fft!(Hk_grid, (3, 4, 5); flags=FFTW.ESTIMATE) : nothing
    HmnR_block = myrank == 0 ? zeros(ComplexF64, Ngsize, Ngsize, block_cells) : nothing

    if myrank == 0
        dense_GiB = _HmnR_memory_GiB(Int128(Ngsize) * Ngsize * NCell * spinsize)
        work_GiB = _HmnR_memory_GiB(
            Int128(Ngsize) * Ngsize * Nkpt +
            Int128(BANDNUM) * Ngsize +
            Int128(Ngsize) * Ngsize * block_cells,
        )
        @printf("\tHmnR dense allocation avoided: %.3f GiB per MPI rank\n", dense_GiB)
        @printf(
            "\tHmnR FFT workspace: %.3f GiB on rank 0 (%d cells/block, %d ranks)\n",
            work_GiB,
            block_cells,
            nprocs,
        )
        @printf("\tHmnR transform: %d x %d x %d FFT\n", kmesh1, kmesh2, kmesh3)
        flush(stdout)
    end

    normalization = inv(Float64(Nkpt))
    for spin = 1:spinsize
        for ik = 1:MPI_Nkpt
            U = Umnk[spin][ik]
            @inbounds for pst = 1:Ngsize, μ = 1:BANDNUM
                weighted_conj_U[μ,pst] =
                    Enk[spin][ik][μ+MinN-1] * conj(U[μ,pst])
            end

            # Preserve the historical storage convention
            # Hk[q,p] = sum_mu E_mu * conj(U[mu,p]) * U[mu,q].
            mul!(@view(Hk_local[:, :, ik]), transpose(U), weighted_conj_U)
        end

        if myrank == 0
            MPI.Gatherv!(
                MPI.IN_PLACE,
                MPI.VBuffer(Hk_global, receive_counts),
                comm;
                root=0,
            )
        else
            MPI.Gatherv!(Hk_local, nothing, comm; root=0)
        end

        if myrank == 0
            @printf("\tHmnR FFT spin %d/%d\n", spin, spinsize)
            flush(stdout)
            fft_plan * Hk_grid

            for first_cell = 1:block_cells:NCell
                last_cell = min(first_cell + block_cells - 1, NCell)
                cells_in_block = last_cell - first_cell + 1

                @inbounds for block_cell = 1:cells_in_block
                    cell = first_cell + block_cell - 1
                    R1, R2, R3 = cell_list_ijk[cell]
                    i1 = mod(Int(R1), kmesh1) + 1
                    i2 = mod(Int(R2), kmesh2) + 1
                    i3 = mod(Int(R3), kmesh3) + 1
                    source = @view Hk_grid[:, :, i3, i2, i1]
                    destination = @view HmnR_block[:, :, block_cell]
                    @. destination = normalization * source
                end

                _Write_HmnR_block!(
                    HmnR_dataset,
                    HmnR_block,
                    first_cell,
                    last_cell,
                    spin,
                    NCell,
                    Ngsize,
                )
            end
        end
    end

    return nothing
end


function _Calc_HmnR_stream_direct!(
    HmnR_dataset,
    spinsize,
    MinN,
    MaxN,
    NCell,
    Ngsize,
    cell_list_ijk,
    Umnk,
    Enk,
    kpoints::KPoints,
    block_cells,
)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nkpt = kpoints.Nkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    BANDNUM = MaxN - MinN + 1

    Hk = zeros(ComplexF64, Ngsize, Ngsize, MPI_Nkpt)
    weighted_conj_U = zeros(ComplexF64, BANDNUM, Ngsize)
    phase_block = zeros(ComplexF64, MPI_Nkpt, block_cells)
    HmnR_block = zeros(ComplexF64, Ngsize, Ngsize, block_cells)
    Hk_2D = reshape(Hk, Int(Ngsize) * Int(Ngsize), MPI_Nkpt)
    HmnR_block_2D = reshape(HmnR_block, Int(Ngsize) * Int(Ngsize), block_cells)

    if myrank == 0
        dense_GiB = _HmnR_memory_GiB(Int128(Ngsize) * Ngsize * NCell * spinsize)
        work_GiB = _HmnR_memory_GiB(
            Int128(Ngsize) * Ngsize * MPI_Nkpt +
            Int128(BANDNUM) * Ngsize +
            Int128(block_cells) * (Int128(Ngsize) * Ngsize + MPI_Nkpt),
        )
        @printf("\tHmnR dense allocation avoided: %.3f GiB per MPI rank\n", dense_GiB)
        @printf(
            "\tHmnR streaming workspace: %.3f GiB per MPI rank (%d cells/block, %d ranks)\n",
            work_GiB,
            block_cells,
            nprocs,
        )
    end

    for spin = 1:spinsize
        for ik = 1:MPI_Nkpt
            U = Umnk[spin][ik]
            @inbounds for pst = 1:Ngsize, μ = 1:BANDNUM
                weighted_conj_U[μ,pst] = Enk[spin][ik][μ+MinN-1] * conj(U[μ,pst])
            end

            # Preserve the historical storage convention
            # Hk[q,p] = sum_mu E_mu * conj(U[mu,p]) * U[mu,q].
            mul!(@view(Hk[:, :, ik]), transpose(U), weighted_conj_U)
        end

        for first_cell = 1:block_cells:NCell
            last_cell = min(first_cell + block_cells - 1, NCell)
            cells_in_block = last_cell - first_cell + 1

            for block_cell = 1:cells_in_block
                cell = first_cell + block_cell - 1
                @inbounds for ik = 1:MPI_Nkpt
                    kRn = dot(MPI_kpts[ik], cell_list_ijk[cell])
                    phase_block[ik,block_cell] = cispi(-2*kRn)
                end
            end

            HmnR_view = @view HmnR_block_2D[:, 1:cells_in_block]
            if MPI_Nkpt == 0
                fill!(HmnR_view, 0.0)
            else
                mul!(HmnR_view,Hk_2D,@view(phase_block[:, 1:cells_in_block]),inv(Float64(Nkpt)),0.0)
            end

            # Only rank 0 needs the completed block.  Using Reduce! instead of
            # Allreduce! avoids replicating the full real-space Hamiltonian on
            # every MPI process.
            MPI.Reduce!(HmnR_block, MPI.SUM, comm; root=0)

            if myrank == 0
                _Write_HmnR_block!(HmnR_dataset, HmnR_block, first_cell, last_cell, spin, NCell, Ngsize)
            end
        end
    end

    return nothing
end
