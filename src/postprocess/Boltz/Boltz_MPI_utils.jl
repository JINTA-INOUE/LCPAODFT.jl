"""Return the contiguous work range owned by a zero-based MPI rank."""
function _boltz_work_range(nwork::Integer, nworkers::Integer, worker::Integer)
    count, remainder = divrem(Int(nwork), Int(nworkers))
    first_index = Int(worker) * count + min(Int(worker), remainder) + 1
    local_count = count + (worker < remainder ? 1 : 0)
    return first_index:(first_index + local_count - 1)
end


"""MPI-distributed Gamma-centred regular k mesh for the Boltzmann solver."""
struct BoltzKPoints
    AllNkpt::Int
    MPI_Nkpt::Int
    kmesh::NTuple{3,Int}
    MPI_kpts::Vector{Vector{Float64}}
    MPI_krange::Vector{UnitRange{Int}}
    MPkpts::Vector{Int}
end


"""Build only the k-point coordinates owned by the calling MPI rank."""
function _boltz_kpoints(kmesh)
    kmesh_tuple = Tuple(Int.(kmesh))
    length(kmesh_tuple) == 3 || error("kmesh must contain three dimensions")
    all(>(0), kmesh_tuple) || error("all kmesh dimensions must be positive")

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    Nkpt = prod(kmesh_tuple)
    nprocs <= Nkpt || error(
        "the number of MPI processes ($nprocs) must not exceed " *
        "the number of k points ($Nkpt)",
    )

    ranges = [_boltz_work_range(Nkpt, nprocs, rank)
              for rank = 0:nprocs - 1]
    local_range = ranges[myrank + 1]
    local_kpoints = Vector{Vector{Float64}}(undef, length(local_range))
    kmesh1, kmesh2, kmesh3 = kmesh_tuple

    for (local_k, global_k) in enumerate(local_range)
        ik, jk, kk = _boltz_kindices(global_k, kmesh2, kmesh3)
        local_kpoints[local_k] = [
            (ik - 1) / kmesh1,
            (jk - 1) / kmesh2,
            (kk - 1) / kmesh3,
        ]
    end

    offsets = [first(range) - 1 for range in ranges]
    return BoltzKPoints(Nkpt, length(local_range), kmesh_tuple,
                        local_kpoints, ranges, offsets)
end


"""Convert a linear k-point index to its three regular-grid indices."""
@inline function _boltz_kindices(kp::Integer, kmesh2::Integer,
                                 kmesh3::Integer)
    k0 = Int(kp) - 1
    ik = div(k0, Int(kmesh2) * Int(kmesh3)) + 1
    remainder = rem(k0, Int(kmesh2) * Int(kmesh3))
    jk = div(remainder, Int(kmesh3)) + 1
    kk = rem(remainder, Int(kmesh3)) + 1
    return ik, jk, kk
end


"""Convert three regular-grid indices to a linear k-point index."""
@inline function _boltz_kindex(ik::Integer, jk::Integer, kk::Integer,
                               kmesh2::Integer, kmesh3::Integer)
    return (Int(ik) - 1) * Int(kmesh2) * Int(kmesh3) +
           (Int(jk) - 1) * Int(kmesh3) + Int(kk)
end


"""Write the eight periodic vertex indices of one regular-grid cell."""
function _boltz_cell_vertices!(vertices, kp::Integer, kmesh)
    kmesh1, kmesh2, kmesh3 = Int.(kmesh)
    i, j, k = _boltz_kindices(kp, kmesh2, kmesh3)
    i -= 1
    j -= 1
    k -= 1

    @inbounds for i_in = 0:1, j_in = 0:1, k_in = 0:1
        ii = mod(i + i_in, kmesh1) + 1
        jj = mod(j + j_in, kmesh2) + 1
        kk = mod(k + k_in, kmesh3) + 1
        vertex = 4 * i_in + 2 * j_in + k_in + 1
        vertices[vertex] = _boltz_kindex(ii, jj, kk, kmesh2, kmesh3)
    end
    return vertices
end


"""Communication and indexing metadata for the k-point halo."""
struct BoltzHaloPlan
    cell_range::UnitRange{Int}
    local_first::Int
    nlocal::Int
    remote_kpoints::Vector{Int}
    vertex_positions::Matrix{Int}
    request_send_counts::Vector{Cint}
    request_recv_counts::Vector{Cint}
    incoming_requests::Vector{Int}
end


"""A rank-local k-point array plus only the remotely required halo blocks."""
struct BoltzHaloData{T,N}
    owned::Array{T,N}
    halo::Array{T,N}
end


@inline function _boltz_halo_get(data::BoltzHaloData{T,3}, i::Integer,
                                 s::Integer, k::Integer) where {T}
    nlocal = size(data.owned, 3)
    return k <= nlocal ? data.owned[i, s, k] :
           data.halo[i, s, k - nlocal]
end


@inline function _boltz_halo_get(data::BoltzHaloData{T,4}, i::Integer,
                                 j::Integer, s::Integer,
                                 k::Integer) where {T}
    nlocal = size(data.owned, 4)
    return k <= nlocal ? data.owned[i, j, s, k] :
           data.halo[i, j, s, k - nlocal]
end


function _boltz_halo_plan(kpoints::BoltzKPoints, kmesh)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    Nkpt = prod(Int.(kmesh))
    Nkpt == kpoints.AllNkpt || error("k-point mesh and ownership disagree")

    cell_range = kpoints.MPI_krange[myrank + 1]
    nlocal = length(cell_range)
    local_first = first(cell_range)
    local_last = last(cell_range)
    owner_ends = [last(range) for range in kpoints.MPI_krange]
    request_sets = [Set{Int}() for _ = 1:nprocs]
    cell_vertices = zeros(Int, 8)

    for kp in cell_range
        _boltz_cell_vertices!(cell_vertices, kp, kmesh)
        for global_k in cell_vertices
            if !(local_first <= global_k <= local_last)
                owner = searchsortedfirst(owner_ends, global_k)
                owner <= nprocs || error("no owner for k-point $global_k")
                push!(request_sets[owner], global_k)
            end
        end
    end

    requests_by_owner = [sort!(collect(indices)) for indices in request_sets]
    request_send_counts = Cint[length(indices) for indices in requests_by_owner]
    remote_kpoints = Int[]
    for indices in requests_by_owner
        append!(remote_kpoints, indices)
    end

    request_recv_counts = Vector{Cint}(undef, nprocs)
    MPI.Alltoall!(MPI.UBuffer(request_send_counts, 1),
                  MPI.UBuffer(request_recv_counts, 1), comm)
    incoming_count = sum(Int(count) for count in request_recv_counts)
    incoming_requests = Vector{Int}(undef, incoming_count)
    MPI.Alltoallv!(MPI.VBuffer(remote_kpoints, request_send_counts),
                   MPI.VBuffer(incoming_requests, request_recv_counts), comm)

    remote_positions = Dict{Int,Int}(
        global_k => nlocal + position
        for (position, global_k) in enumerate(remote_kpoints)
    )
    vertex_positions = Matrix{Int}(undef, 8, nlocal)
    for (cell_position, kp) in enumerate(cell_range)
        _boltz_cell_vertices!(cell_vertices, kp, kmesh)
        for vertex = 1:8
            global_k = cell_vertices[vertex]
            vertex_positions[vertex, cell_position] =
                local_first <= global_k <= local_last ?
                global_k - local_first + 1 : remote_positions[global_k]
        end
    end

    return BoltzHaloPlan(cell_range, local_first, nlocal, remote_kpoints,
                         vertex_positions, request_send_counts,
                         request_recv_counts, incoming_requests)
end


function _boltz_scaled_counts(counts, multiplier::Integer)
    scaled = Vector{Cint}(undef, length(counts))
    for i in eachindex(counts)
        count = Int(counts[i]) * Int(multiplier)
        count <= typemax(Cint) || error("MPI halo message exceeds Cint count")
        scaled[i] = Cint(count)
    end
    return scaled
end


"""Exchange only remotely required k-point blocks."""
function _boltz_exchange_halo(local_data::Array{T,N},
                              plan::BoltzHaloPlan) where {T,N}
    size(local_data, N) == plan.nlocal ||
        error("local k-point array does not match the halo plan")

    blocksize = prod(size(local_data)[1:(N - 1)])
    halo_dims = ntuple(
        dimension -> dimension == N ? length(plan.remote_kpoints) :
                     size(local_data, dimension),
        N,
    )
    halo_data = similar(local_data, halo_dims)

    send_values = Vector{T}(undef, blocksize * length(plan.incoming_requests))
    for (request, global_k) in enumerate(plan.incoming_requests)
        local_k = global_k - plan.local_first + 1
        1 <= local_k <= plan.nlocal ||
            error("received a halo request for non-local k-point $global_k")
        copyto!(send_values, (request - 1) * blocksize + 1,
                local_data, (local_k - 1) * blocksize + 1, blocksize)
    end

    send_counts = _boltz_scaled_counts(plan.request_recv_counts, blocksize)
    receive_counts = _boltz_scaled_counts(plan.request_send_counts, blocksize)
    MPI.Alltoallv!(MPI.VBuffer(send_values, send_counts),
                   MPI.VBuffer(halo_data, receive_counts), MPI.COMM_WORLD)
    return BoltzHaloData(local_data, halo_data)
end


function _boltz_report_distributed_memory(plan::BoltzHaloPlan, Nkpt, arrays...)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    nhalo = length(plan.remote_kpoints)
    payload_bytes = sum(
        sizeof(array.owned) + sizeof(array.halo) for array in arrays)
    replicated_bytes = sum(
        prod(size(array.owned)[1:(ndims(array.owned) - 1)]) *
        Int(Nkpt) * sizeof(eltype(array.owned))
        for array in arrays
    )

    owned_min = MPI.Reduce(plan.nlocal, MPI.MIN, comm; root=0)
    owned_max = MPI.Reduce(plan.nlocal, MPI.MAX, comm; root=0)
    halo_min = MPI.Reduce(nhalo, MPI.MIN, comm; root=0)
    halo_max = MPI.Reduce(nhalo, MPI.MAX, comm; root=0)
    stored_total = MPI.Reduce(plan.nlocal + nhalo, MPI.SUM, comm; root=0)
    payload_total = MPI.Reduce(payload_bytes, MPI.SUM, comm; root=0)
    payload_max = MPI.Reduce(payload_bytes, MPI.MAX, comm; root=0)

    if myrank == 0
        replicated_total = replicated_bytes * nprocs
        println("<Boltz distributed k-point memory>")
        println("\towned k points per rank: $owned_min .. $owned_max")
        println("\thalo k points per rank: $halo_min .. $halo_max")
        println("\taggregate stored k-point slots: $stored_total " *
                "(full replication: $(Int(Nkpt) * nprocs))")
        println("\taggregate Enk/Vnk/EVec payload: " *
                "$(round(payload_total / 2.0^20; digits=3)) MiB " *
                "(full replication: " *
                "$(round(replicated_total / 2.0^20; digits=3)) MiB)")
        println("\tmaximum Enk/Vnk/EVec payload per rank: " *
                "$(round(payload_max / 2.0^20; digits=3)) MiB")
    end
    return nothing
end
