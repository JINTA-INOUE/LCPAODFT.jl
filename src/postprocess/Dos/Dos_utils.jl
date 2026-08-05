"""Return the contiguous work range owned by a zero-based worker id."""
function _dos_work_range(nwork::Integer, nworkers::Integer, worker::Integer)
    count, remainder = divrem(Int(nwork), Int(nworkers))
    first_index = Int(worker) * count + min(Int(worker), remainder) + 1
    local_count = count + (worker < remainder ? 1 : 0)
    return first_index:(first_index + local_count - 1)
end


"""MPI-distributed regular k mesh used only by the DOS implementation."""
struct DosKPoints
    AllNkpt::Int
    MPI_Nkpt::Int
    MPI_kpts::Vector{Vector{Float64}}
    MPI_krange::Vector{UnitRange{Int}}
    MPkpts::Vector{Int}
end


"""Build only the k-point coordinates owned by the calling MPI rank."""
function _dos_kpoints(kmesh)
    kmesh1, kmesh2, kmesh3 = Int.(kmesh)
    all(dimension -> dimension > 0, (kmesh1, kmesh2, kmesh3)) ||
        error("all kmesh dimensions must be positive")

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    Nkpt = kmesh1 * kmesh2 * kmesh3
    ranges = [_dos_work_range(Nkpt, nprocs, rank)
              for rank = 0:nprocs - 1]
    local_range = ranges[myrank + 1]
    local_kpoints = Vector{Vector{Float64}}(undef, length(local_range))

    for (local_k, global_k) in enumerate(local_range)
        ik, jk, kk = _dos_kindices(global_k, kmesh2, kmesh3)
        local_kpoints[local_k] = [
            (ik - 1) / kmesh1,
            (jk - 1) / kmesh2,
            (kk - 1) / kmesh3,
        ]
    end

    offsets = [first(range) - 1 for range in ranges]
    return DosKPoints(Nkpt, length(local_range), local_kpoints,
                      ranges, offsets)
end


"""Convert the linear DOS k-point index to its three grid indices."""
@inline function _dos_kindices(kp::Integer, kmesh2::Integer, kmesh3::Integer)
    k0 = Int(kp) - 1
    ik = div(k0, Int(kmesh2) * Int(kmesh3)) + 1
    remainder = rem(k0, Int(kmesh2) * Int(kmesh3))
    jk = div(remainder, Int(kmesh3)) + 1
    kk = rem(remainder, Int(kmesh3)) + 1
    return ik, jk, kk
end


"""Convert three DOS grid indices to the linear k-point index."""
@inline function _dos_kindex(ik::Integer, jk::Integer, kk::Integer,
                            kmesh2::Integer, kmesh3::Integer)
    return (Int(ik) - 1) * Int(kmesh2) * Int(kmesh3) +
           (Int(jk) - 1) * Int(kmesh3) + Int(kk)
end


"""Write the eight periodic vertex k-point indices of one regular-grid cell."""
function _dos_cell_vertices!(vertices, kp::Integer, kmesh)
    kmesh1, kmesh2, kmesh3 = Int.(kmesh)
    i, j, k = _dos_kindices(kp, kmesh2, kmesh3)
    i -= 1
    j -= 1
    k -= 1

    @inbounds for i_in = 0:1, j_in = 0:1, k_in = 0:1
        ii = mod(i + i_in, kmesh1) + 1
        jj = mod(j + j_in, kmesh2) + 1
        kk = mod(k + k_in, kmesh3) + 1
        vertex = 4 * i_in + 2 * j_in + k_in + 1
        vertices[vertex] = _dos_kindex(ii, jj, kk, kmesh2, kmesh3)
    end
    return vertices
end


"""Communication and indexing metadata for the DOS k-point halo."""
struct DosHaloPlan
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
struct DosHaloData{T,N}
    owned::Array{T,N}
    halo::Array{T,N}
end


@inline function _dos_halo_get(data::DosHaloData{T,3}, i::Integer,
                               j::Integer, k::Integer) where {T}
    nlocal = size(data.owned, 3)
    return k <= nlocal ? data.owned[i, j, k] :
           data.halo[i, j, k - nlocal]
end


@inline function _dos_halo_get(data::DosHaloData{T,4}, i::Integer,
                               j::Integer, s::Integer,
                               k::Integer) where {T}
    nlocal = size(data.owned, 4)
    return k <= nlocal ? data.owned[i, j, s, k] :
           data.halo[i, j, s, k - nlocal]
end


function _dos_halo_plan(kpoints::DosKPoints, kmesh)
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    Nkpt = prod(Int.(kmesh))
    Nkpt == kpoints.AllNkpt || error("k-point mesh and ownership disagree")

    cell_range = kpoints.MPI_krange[myrank + 1]
    nlocal = length(cell_range)
    local_first = first(cell_range)
    local_last = local_first + nlocal - 1
    owner_ends = [last(range) for range in kpoints.MPI_krange]
    request_sets = [Set{Int}() for _ = 1:nprocs]
    cell_vertices = zeros(Int, 8)

    for kp in cell_range
        _dos_cell_vertices!(cell_vertices, kp, kmesh)
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
        _dos_cell_vertices!(cell_vertices, kp, kmesh)
        for vertex = 1:8
            global_k = cell_vertices[vertex]
            vertex_positions[vertex, cell_position] =
                local_first <= global_k <= local_last ?
                global_k - local_first + 1 : remote_positions[global_k]
        end
    end

    return DosHaloPlan(cell_range, local_first, nlocal, remote_kpoints,
                       vertex_positions, request_send_counts,
                       request_recv_counts, incoming_requests)
end


function _dos_scaled_counts(counts, multiplier::Integer)
    scaled = Vector{Cint}(undef, length(counts))
    for i in eachindex(counts)
        count = Int(counts[i]) * Int(multiplier)
        count <= typemax(Cint) || error("MPI halo message exceeds Cint count")
        scaled[i] = Cint(count)
    end
    return scaled
end


"""Exchange only remotely required blocks while retaining rank-local storage."""
function _dos_exchange_halo(local_data::Array{T,N}, plan::DosHaloPlan) where {T,N}
    size(local_data, N) == plan.nlocal ||
        error("local k-point array does not match the halo plan")

    blocksize = prod(size(local_data)[1:(N - 1)])
    halo_dims = ntuple(
        dimension -> dimension == N ?
            length(plan.remote_kpoints) :
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

    send_counts = _dos_scaled_counts(plan.request_recv_counts, blocksize)
    receive_counts = _dos_scaled_counts(plan.request_send_counts, blocksize)
    MPI.Alltoallv!(MPI.VBuffer(send_values, send_counts),
                   MPI.VBuffer(halo_data, receive_counts),
                   MPI.COMM_WORLD)
    return DosHaloData(local_data, halo_data)
end


function _dos_report_distributed_memory(plan::DosHaloPlan, Nkpt, arrays...)
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
        println("<DOS distributed k-point memory>")
        println("\towned k points per rank: $owned_min .. $owned_max")
        println("\thalo k points per rank: $halo_min .. $halo_max")
        println("\taggregate stored k-point slots: $stored_total " *
                "(full replication: $(Int(Nkpt) * nprocs))")
        println("\taggregate Enk/EVec payload: " *
                "$(round(payload_total / 2.0^20; digits=3)) MiB " *
                "(full replication: " *
                "$(round(replicated_total / 2.0^20; digits=3)) MiB)")
        println("\tmaximum Enk/EVec payload per rank: " *
                "$(round(payload_max / 2.0^20; digits=3)) MiB")
    end
    return nothing
end


function Get_Angle_spin(material::LCPAO_model)

    SpinPol = material.SpinPol
    if SpinPol ≠ "nc"
        error("not support $SpinPol case")
    end

    Nspin = 4
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    Total_NumOrbs = material.Total_NumOrbs
    OLP = material.OLP
    DM = material.DM


    MulP = zeros(Float64, 4)
    Angle_Spin = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Angle_Spin[atom] = zeros(Float64, 2)
    end

    for atom = 1:Natom

        fill!(MulP, 0.0)
        for Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            for spin = 1:Nspin, ist = 1:Total_NumOrbs[atom]
                tmp0 = 0.0
                for jst = 1:Total_NumOrbs[jatom]
                    tmp0 += DM[spin][atom][Rn][ist][jst] * OLP[atom][Rn][ist][jst]
                end

                if spin == 4
                    MulP[spin] -= tmp0
                else
                    MulP[spin] += tmp0
                end
            end
        end


        Nup, Ndown, theta, phi = EulerAngle_Spin(MulP[1], MulP[2], MulP[3], MulP[4])

        MulP[1] = Nup
        MulP[2] = Ndown
        MulP[3] = theta
        MulP[4] = phi

        Angle_Spin[atom][1] = MulP[3]
        Angle_Spin[atom][2] = MulP[4]
    end


    return Angle_Spin
end


function Calc_Band_size(material::LCPAO_model, Dos_Erange)

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


    # Find iemin, iemax at Γ point
    # iemin : minimal band index
    # iemax : maximum band index
    iemin = 1
    iemax = 1
    n1min = 1
    kpts_zeros = zeros(Float64, 3)
    EΓ = zeros(Float64, Nfsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    tmpH = zeros(ComplexF64, fsize, fsize)

    for spin = 1:spinsize

        if SpinPol ∈ ("off", "on")
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        elseif SpinPol == "nc"
            HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        end

        iemin0 = 1
        n1min = ifelse(n1min<Nfsize, Nfsize, n1min)
        
        for μ = 1:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[1])
                iemin0 = μ - 1
                break
            end
        end

        iemin0 = ifelse(iemin0<1, 1, iemin0)
        iemax0 = Nfsize

        for μ = iemin0:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[2])
                iemax0 = μ
                break
            end
        end

        iemax0 = ifelse(iemax0>Nfsize, Nfsize, iemax0)
        iemin = ifelse(iemin>iemin0, iemin0, iemin)
        iemax = ifelse(iemax<iemax0, iemax0, iemax)
    end

    if SpinPol ∈ ("off", "on")
        iemin -= max(div(fsize, 20), 10)
        iemax += max(div(fsize, 20), 10)
    elseif SpinPol == "nc"
        iemin -= max(div(fsize, 10), 10)
        iemax += max(div(fsize, 10), 10)
    end

    iemin = ifelse(iemin<1, 1, iemin)
    iemax = ifelse(iemax>Nfsize, Nfsize, iemax)
    iemax = ifelse(iemax>n1min, n1min, iemax)

    return iemin, iemax
end


function Calc_Band_size(material::CWF_model, Dos_Erange)

    SpinPol = material.SpinPol
    Nwann = material.Ngsize
    NCell = material.NCell
    HmnR = material.HmnR
    ChemP = material.ChemP
    spinsize = ifelse(SpinPol == "on", 2, 1)    


    # Find iemin, iemax at Γ point
    # iemin : minimal band index
    # iemax : maximum band index
    iemin = 1
    iemax = 1
    n1min = 1
    EΓ = zeros(Float64, Nwann)
    H = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize
        for cell = 1:NCell, ist = 1:Nwann, jst = 1:Nwann
            H[ist,jst] += HmnR[jst,ist,cell,spin]
        end

        EΓ = eigvals(Hermitian(H))

        iemin0 = 1
        n1min = ifelse(n1min<Nwann, Nwann, n1min)
        
        for μ = 1:Nwann
            if EΓ[μ] > (ChemP + Dos_Erange[1])
                iemin0 = μ - 1
                break
            end
        end

        iemin0 = ifelse(iemin0<1, 1, iemin0)
        iemax0 = Nwann

        for μ = iemin0:Nwann
            if EΓ[μ] > (ChemP + Dos_Erange[2])
                iemax0 = μ
                break
            end
        end

        iemax0 = ifelse(iemax0>Nwann, Nwann, iemax0)
        iemin = ifelse(iemin>iemin0, iemin0, iemin)
        iemax = ifelse(iemax<iemax0, iemax0, iemax)
    end

    if SpinPol ∈ ("off", "on")
        iemin -= max(div(Nwann, 20), 10)
        iemax += max(div(Nwann, 20), 10)
    elseif SpinPol == "nc"
        iemin -= max(div(div(Nwann,2), 10), 10)
        iemax += max(div(div(Nwann,2), 10), 10)
    end

    iemin = ifelse(iemin<1, 1, iemin)
    iemax = ifelse(iemax>Nwann, Nwann, iemax)
    iemax = ifelse(iemax>n1min, n1min, iemax)

    return iemin, iemax
end


function Calc_PDos_Atom_proj(filename::String, material::LCPAO_model, Dos_Erange, DosE, Dos)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    spinsize = ifelse(SpinPol == "off", 1, 2)
    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]
    Dos_N = length(DosE)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)


    for atom = 1:Natom
        fill!(DosSum, 0.0)
        for spin = 1:spinsize
            for q = 1:Dos_N
                s1 = 0.0
                s2 = 0.0
                for ie = 2:2:q-1
                    DosSum[ie,spin] = 0.0
                    for ist = 1:Total_NumOrbs[atom]
                        DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                    end
                    s1 += DosSum[ie,spin]
                end

                for ie = 3:2:q-1
                    DosSum[ie,spin] = 0.0
                    for ist = 1:Total_NumOrbs[atom]
                        DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                    end
                    s2 += DosSum[ie,spin]
                end
                ssum[q,spin] = (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
            end
        end

        println("Write $(filename).PDOS.Tetrahedron.atom$(atom)")
        Dos_data = open("$(filename).PDOS.Tetrahedron.atom$(atom)", "w")

        for ie = 1:Dos_N
            fill!(DSum, 0.0)
            for ist = 1:Total_NumOrbs[atom], spin = 1:spinsize
                DSum[spin] += Dos[spin][atom][ist][ie]
            end

            if SpinPol ∈ ("on", "nc")
                @printf(Dos_data, "%lf %lf %lf %lf %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
            else
                @printf(Dos_data, "%lf %lf %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
            end
        end
        close(Dos_data)
    end
end


function Calc_PDos_Orbital_proj(filename::String, material::LCPAO_model, Spe_Num_Relation, Spe_Num_Basis, Dos_Erange, DosE, Dos)
    
    Natom = material.Natom
    SpinPol = material.SpinPol
    atom2spe = material.atom2spe
    Total_NumOrbs = material.Total_NumOrbs
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Lname = ["s", "p", "d", "f"]
    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]
    Dos_N = length(DosE)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)

    

    for atom = 1:Natom, L = 0:3
        spe = atom2spe[atom]
        if Spe_Num_Basis[spe][L+1] > 0
            for M = 0:2*L
                LM = 100*L + M + 1

                fill!(DosSum, 0.0)
                for spin = 1:spinsize
                    for q = 1:Dos_N
                        s1 = 0.0
                        s2 = 0.0
                        for ie = 2:2:q-1
                            DosSum[ie,spin] = 0.0
                            for ist = 1:Total_NumOrbs[atom]
                                if LM == Spe_Num_Relation[spe][ist]
                                    DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                                end
                            end
                            s1 += DosSum[ie,spin]
                        end
                        for ie = 3:2:q-1
                            DosSum[ie,spin] = 0.0
                            for ist = 1:Total_NumOrbs[atom]
                                if LM == Spe_Num_Relation[spe][ist]
                                    DosSum[ie,spin] += Dos[spin][atom][ist][ie]
                                end
                            end
                            s2 += DosSum[ie,spin]
                        end
                        ssum[q,spin] = (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
                    end
                end

                println("Write $(filename).PDOS.Tetrahedron.atom$(atom).$(Lname[L+1])$(M+1)")
                PDos_data = open("$(filename).PDOS.Tetrahedron.atom$(atom).$(Lname[L+1])$(M+1)", "w")
                for ie = 1:Dos_N
                    fill!(DSum, 0.0)
                    for ist = 1:Total_NumOrbs[atom]
                        if LM == Spe_Num_Relation[spe][ist]
                            for spin = 1:spinsize
                                DSum[spin] += Dos[spin][atom][ist][ie]
                            end
                        end
                    end

                    if SpinPol ∈ ("on", "nc")
                        @printf(PDos_data, "%lf  %lf  %lf  %lf  %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
                    else
                        @printf(PDos_data, "%lf  %lf  %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
                    end
                end
                close(PDos_data)
            end
        end
    end
end


function Calc_PDos_Orbital_proj(filename::String, material::CWF_model, Dos_Erange, DosE, Dos)
    
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    gsize = material.gsize
    DosEmin, DosEmax = Dos_Erange
    Dos_N = length(DosE)
    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree

    DSum = zeros(Float64, spinsize)
    DosSum = zeros(Float64, Dos_N, spinsize)
    ssum = zeros(Float64, Dos_N, spinsize)

    

    for ist = 1:gsize
        for spin = 1:spinsize, q = 1:Dos_N
            s1 = 0.0
            s2 = 0.0
            for ie = 2:2:q-1
                s1 += Dos[spin][ist][ie]
            end
                
            for ie = 3:2:q-1
                s2 += Dos[spin][ist][ie]
            end
            ssum[q,spin] =
                (DosSum[begin,spin] + 4*s1 + 2*s2 + DosSum[q,spin])*h/3
        end


        println("Write $(filename).PDOS.Tetrahedron.orb$(ist)")
        PDos_data = open("$(filename).PDOS.Tetrahedron.orb$(ist)", "w")
        for ie = 1:Dos_N
            for spin = 1:spinsize
                DSum[spin] = Dos[spin][ist][ie]
            end

            if SpinPol ∈ ("on", "nc")
                @printf(PDos_data, "%lf  %lf  %lf  %lf  %lf\n", DosE[ie]*eV2Hartree, DSum[1], -DSum[2], ssum[ie,1], ssum[ie,2])
            else
                @printf(PDos_data, "%lf  %lf  %lf\n", DosE[ie]*eV2Hartree, 2*DSum[1], 2*ssum[ie,1])
            end
        end
        close(PDos_data)
    end
end
