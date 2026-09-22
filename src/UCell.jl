struct System_Grid
    CpyCell::Int32
    Natom::Int32
    atom2spe::Vector{Int32}
    Latvecs::Matrix{Float64}
    atv::Vector{Vector{Float64}}
    atv_ijk::Vector{Vector{Int32}}
    Gxyz::Vector{Vector{Float64}}
    GridVol::Float64
    Grid_Origin::Vector{Float64}
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    ncn::Vector{Vector{Int32}}
    Dis::Vector{Vector{Float64}}
    RMI::Vector{Vector{Vector{Int32}}}
    Total_Hsize::Int32
    MPI_Hsize::Vector{Int32}
    MPHks::Vector{Int32}
    Nloop::Int32
    MPI_size::Int32
    MPI_atom::Vector{Int32}
    MPI_FNAN::Vector{Int32}
    MPI_natn::Vector{Int32}
    MPI_ncn::Vector{Int32}
    MPI_Hoffset::Vector{Int64}
    Atoms_Cut1::Vector{Float64}
    Total_NumOrbs::Vector{Int32}
    MP::Vector{Int32}
    Ngrid::Tuple{Int32,Int32,Int32}
end


struct UCell
    system_grid::System_Grid
    GridN_Atom::Vector{Int32}
    GridListAtom::Vector{Vector{Int32}}
    CellListAtom::Vector{Vector{Int32}} 
    MPI_NumOLG::Vector{Int32}
    MPI_GListTAtoms1::Vector{Vector{Int32}}
    MPI_GListTAtoms2::Vector{Vector{Int32}}
    density_scratch::Matrix{Float64}
    density_matrix_scratch::Array{Float64,3}
    density_orbital_scratch::NTuple{3,Matrix{Float64}}
    hamiltonian_orbital_scratch::NTuple{3,Matrix{Float64}}
    hamiltonian_product_scratch::Matrix{Float64}
    Ngrid::Tuple{Int32,Int32,Int32}
end



"""
Assign atom-pair work to MPI ranks.

Pairs sharing the same two atoms are considered together. Several deterministic
orders are tried; the partition with the smallest atom-grid memory footprint is
selected while respecting work and pair-count balance limits.
"""
function _split_weighted_vertex_cut(pair_atoms::AbstractVector{<:Integer},
                                    pair_neighbors::AbstractVector{<:Integer},
                                    pair_cells::AbstractVector{<:Integer},
                                    weights::AbstractVector{<:Real},
                                    atom_memory::AbstractVector{<:Real},
                                    nparts::Integer;
                                    max_load_ratio::Real=1.05,
                                    max_pair_count_ratio::Real=1.40,
                                    balance_weight::Real=0.05)
    nitems = length(weights)
    nitems >= nparts || throw(ArgumentError(
        "number of work items must be at least the number of MPI ranks"))
    length(pair_atoms) == nitems == length(pair_neighbors) == length(pair_cells) ||
        throw(ArgumentError("atom-pair work arrays must have the same length"))

    total_work = sum(Float64, weights)
    target_work = total_work/nparts
    load_limit = max_load_ratio*target_work
    target_pair_count = nitems/nparts
    pair_count_limit = max(ceil(Int, max_pair_count_ratio*target_pair_count), 1)
    memory_scale = max(sum(Float64, atom_memory), 1.0)

    # Periodic images and the two directions of the same physical atom pair
    # form one affinity group.  Records remain individually assignable so a
    # large group cannot force an excessive load imbalance.
    groups = Dict{Tuple{Int32,Int32},Vector{Int}}()
    for item in eachindex(weights)
        atom = Int32(pair_atoms[item])
        jatom = Int32(pair_neighbors[item])
        key = (min(atom, jatom), max(atom, jatom))
        push!(get!(groups, key, Int[]), item)
    end
    all_group_keys = collect(keys(groups))
    group_work = Dict(key => sum(Float64(weights[item]) for item in groups[key])
                      for key in all_group_keys)
    for items in values(groups)
        sort!(items; by=item -> (-Float64(weights[item]),
                                  pair_atoms[item], pair_neighbors[item],
                                  pair_cells[item]))
    end

    function pair_priority(key, trial)
        value = UInt64(key[1]) | (UInt64(key[2]) << 32)
        value += UInt64(trial)*0x9e3779b97f4a7c15
        value = xor(value, value >> 30)*0xbf58476d1ce4e5b9
        value = xor(value, value >> 27)*0x94d049bb133111eb
        return xor(value, value >> 31)
    end

    candidate_orders = Vector{Vector{Tuple{Int32,Int32}}}()
    for key_function in (
        key -> (key[1], key[2]),
        key -> (key[2], key[1]),
        key -> (key[1] + key[2], key[1], key[2]),
        key -> (key[2] - key[1], key[1], key[2]),
        key -> (-group_work[key], key[1], key[2]),
        key -> (-group_work[key], -key[1], -key[2]),
    )
        push!(candidate_orders, sort(copy(all_group_keys); by=key_function))
    end
    for trial = 1:8
        push!(candidate_orders, sort(copy(all_group_keys);
            by=key -> (-group_work[key], pair_priority(key, trial))))
    end
    for trial = 9:32
        push!(candidate_orders, sort(copy(all_group_keys);
            by=key -> pair_priority(key, trial)))
    end

    function partition_for_order(group_keys)
        buckets = [Int[] for _ = 1:nparts]
        loads = zeros(Float64, nparts)
        present = falses(nparts, length(atom_memory))

        for key in group_keys, item in groups[key]
            atom = Int(pair_atoms[item])
            jatom = Int(pair_neighbors[item])
            work = Float64(weights[item])

            feasible = Int[]
            for rank = 1:nparts
                if loads[rank] + work <= load_limit + eps(load_limit) &&
                   length(buckets[rank]) < pair_count_limit
                    push!(feasible, rank)
                end
            end
            isempty(feasible) && append!(feasible, 1:nparts)

            best_rank = feasible[1]
            best_score = (Inf, Inf, typemax(Int))
            for rank in feasible
                added_memory = present[rank, atom] ? 0.0 : atom_memory[atom]
                if atom != jatom && !present[rank, jatom]
                    added_memory += atom_memory[jatom]
                end
                projected_load = (loads[rank] + work)/target_work
                score = (added_memory/memory_scale +
                         balance_weight*projected_load^2,
                         projected_load, rank)
                if score < best_score
                    best_rank = rank
                    best_score = score
                end
            end

            push!(buckets[best_rank], item)
            loads[best_rank] += work
            present[best_rank, atom] = true
            present[best_rank, jatom] = true
        end

        # Empty ranks are undesirable even for unusually indivisible
        # workloads. Moving one light record is sufficient.
        for empty_rank in findall(isempty, buckets)
            donor = argmax(length.(buckets))
            length(buckets[donor]) > 1 || error("unable to give every MPI rank work")
            donor_weights = [Float64(weights[item]) for item in buckets[donor]]
            donor_position = argmin(donor_weights)
            item = splice!(buckets[donor], donor_position)
            push!(buckets[empty_rank], item)
            loads[donor] -= Float64(weights[item])
            loads[empty_rank] += Float64(weights[item])
        end

        fill!(present, false)
        for rank = 1:nparts, item in buckets[rank]
            present[rank, pair_atoms[item]] = true
            present[rank, pair_neighbors[item]] = true
        end
        rank_memory = [sum(atom_memory[atom] for atom in axes(present, 2)
                           if present[rank, atom]) for rank = 1:nparts]
        max_load = maximum(loads)/target_work
        max_pair_count = maximum(length, buckets)/target_pair_count
        work_overload = max(max_load - max_load_ratio, 0.0)
        count_overload = max(max_pair_count - max_pair_count_ratio, 0.0)
        total_overload = work_overload + count_overload
        partition_score = (total_overload > sqrt(eps(Float64)) ? 1 : 0,
                           total_overload, sum(rank_memory), maximum(rank_memory),
                           max_load, max_pair_count)
        return buckets, loads, partition_score
    end

    best_buckets = Vector{Vector{Int}}()
    best_loads = Float64[]
    best_partition_score = nothing
    for group_order in candidate_orders
        buckets, loads, partition_score = partition_for_order(group_order)
        if isnothing(best_partition_score) || partition_score < best_partition_score
            best_buckets = buckets
            best_loads = loads
            best_partition_score = partition_score
        end
    end

    # A deterministic local order keeps equal atom pairs adjacent during all
    # matrix and real-space grid integrations.
    for bucket in best_buckets
        sort!(bucket; by=item -> (
            min(pair_atoms[item], pair_neighbors[item]),
            max(pair_atoms[item], pair_neighbors[item]),
            pair_atoms[item], pair_neighbors[item], pair_cells[item]))
    end

    return best_buckets, best_loads
end


function split_system_grid(Natom, FNAN, natn, ncn, Total_NumOrbs;
                           GridN_Atom=nothing, CpyCell=nothing,
                           GridListAtom=nothing, CellListAtom=nothing,
                           atv_ijk=nothing, max_load_ratio=1.05,
                           max_pair_count_ratio=1.40)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    pair_count = sum(FNAN .+ 1)
    pair_atoms = zeros(Int32, pair_count)
    pair_neighbor_indices = zeros(Int32, pair_count)
    pair_neighbors = zeros(Int32, pair_count)
    pair_cells = zeros(Int32, pair_count)
    pair_hamiltonian_offsets = zeros(Int64, pair_count)
    pair_work = zeros(Float64, pair_count)

    pair = 1
    hamiltonian_offset = Int64(0)
    for atom = 1:Natom, neighbor_index = 1:FNAN[atom]+1
        neighbor = natn[atom][neighbor_index]
        pair_atoms[pair] = atom
        pair_neighbor_indices[pair] = neighbor_index
        pair_neighbors[pair] = neighbor
        pair_cells[pair] = ncn[atom][neighbor_index]
        pair_hamiltonian_offsets[pair] = hamiltonian_offset

        orbital_work = Total_NumOrbs[atom]*Total_NumOrbs[neighbor]
        grid_work = isnothing(GridN_Atom) ? 1 :
            min(GridN_Atom[atom], GridN_Atom[neighbor])
        pair_work[pair] = max(orbital_work*grid_work, 1)

        hamiltonian_offset += orbital_work
        pair += 1
    end

    overlap_inputs_available = !isnothing(CpyCell) &&
        !isnothing(GridN_Atom) && !isnothing(GridListAtom) &&
        !isnothing(CellListAtom) && !isnothing(atv_ijk)
    if nprocs > 1 && overlap_inputs_available
        overlap_counts = zeros(Int32, pair_count)
        if myrank == 0
            Count_AtomOverlap_Grid!(overlap_counts, CpyCell, pair_atoms,
                pair_neighbors, pair_cells, GridN_Atom, GridListAtom,
                CellListAtom, atv_ijk)
        end
        MPI.Bcast!(overlap_counts, 0, comm)
        for pair in eachindex(pair_work)
            orbital_work = Total_NumOrbs[pair_atoms[pair]]*
                           Total_NumOrbs[pair_neighbors[pair]]
            pair_work[pair] = max(
                orbital_work*max(overlap_counts[pair], 1), 1)
        end
    end

    atom_memory = isnothing(GridN_Atom) ? ones(Float64, Natom) :
        Float64.(Total_NumOrbs).*Float64.(GridN_Atom)
    buckets, _ = _split_weighted_vertex_cut(pair_atoms, pair_neighbors,
        pair_cells, pair_work, atom_memory, nprocs;
        max_load_ratio, max_pair_count_ratio)
    local_pairs = buckets[myrank+1]
    MPI_size = length(local_pairs)

    MPI_atom = pair_atoms[local_pairs]
    MPI_FNAN = pair_neighbor_indices[local_pairs]
    MPI_natn = pair_neighbors[local_pairs]
    MPI_ncn = pair_cells[local_pairs]
    MPI_Hoffset = pair_hamiltonian_offsets[local_pairs]

    MPI_Hsize = zeros(Int32, nprocs)
    MPHks = zeros(Int32, nprocs)

    myHsize = sum(Total_NumOrbs[MPI_atom[pair]]*
                  Total_NumOrbs[MPI_natn[pair]]
                  for pair in eachindex(MPI_atom))
    MPI_Hsize[myrank+1] = myHsize
    Total_Hsize = MPI.Allreduce(myHsize, MPI.SUM, comm)
    MPI.Allreduce!(MPI_Hsize, MPI.SUM, comm)

    rank_offset = 0
    for rank = 1:nprocs
        MPHks[rank] = rank_offset
        rank_offset += MPI_Hsize[rank]
    end

    return Total_Hsize, MPI_Hsize, MPHks, pair_count, MPI_size,
           MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Hoffset
end


"""Assemble a rank-local atom-pair vector in canonical `(atom, Rn)` order."""
function assemble_canonical!(global_data::AbstractVector, local_data::AbstractVector,
                             system_grid::System_Grid, comm=MPI.COMM_WORLD)
    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_Hoffset = system_grid.MPI_Hoffset
    Total_NumOrbs = system_grid.Total_NumOrbs

    fill!(global_data, zero(eltype(global_data)))
    local_offset = 0
    @inbounds for loop in eachindex(MPI_atom)
        block_size = Total_NumOrbs[MPI_atom[loop]]*Total_NumOrbs[MPI_natn[loop]]
        copyto!(global_data, MPI_Hoffset[loop] + 1,
                local_data, local_offset + 1, block_size)
        local_offset += block_size
    end
    local_offset == length(local_data) || error("rank-local Hamiltonian size mismatch")
    MPI.Allreduce!(global_data, MPI.SUM, comm)
    return global_data
end


"""Extract canonical atom-pair blocks into the rank-local calculation order."""
function extract_canonical!(local_data::AbstractVector, global_data::AbstractVector,
                            system_grid::System_Grid)
    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_Hoffset = system_grid.MPI_Hoffset
    Total_NumOrbs = system_grid.Total_NumOrbs

    local_offset = 0
    @inbounds for loop in eachindex(MPI_atom)
        block_size = Total_NumOrbs[MPI_atom[loop]]*Total_NumOrbs[MPI_natn[loop]]
        copyto!(local_data, local_offset + 1,
                global_data, MPI_Hoffset[loop] + 1, block_size)
        local_offset += block_size
    end
    local_offset == length(local_data) || error("rank-local Hamiltonian size mismatch")
    return local_data
end


# Set UCell for postprocess
@timeit timer "UCell" function UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs; verbosity=1)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    myrank == 0 && println("<UCell>  Setup Grid ...")

    CpyCell = Int64(0.5*(cbrt(TCpyCell+1)-1))
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)
    Generation_ATV_ijk!(CpyCell, atv_ijk)

    FNAN, natn, ncn, Dis = Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell)

    
    RMI = Get_RMI(Natom, CpyCell, FNAN, natn, ncn)
    if myrank == 0
        system = Check_system(FNAN, ncn, atv_ijk)
        if verbosity >= 1
            println("<Check_System> The system is $system.")
        end
    end
    
    MP = zeros(Int32, Natom+1)
    Sum = 0
    for atom = 1:Natom
        MP[atom+1] = Sum + Total_NumOrbs[atom]
        Sum += Total_NumOrbs[atom]
    end
    GridVol = abs(det(Latvecs))/prod(Ngrid)



    GridN_Atom, GridListAtom, CellListAtom = Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin)

    # split element for MPI
    # one dimensionalization Natom, FNAN, natn, ncn
    Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Hoffset = split_system_grid(
        Natom, FNAN, natn, ncn, Total_NumOrbs;
        GridN_Atom, CpyCell, GridListAtom, CellListAtom, atv_ijk)
    


    system_grid = System_Grid(CpyCell, Natom, atom2spe, Latvecs,
                              atv, atv_ijk, Gxyz, GridVol, Grid_Origin,
                              FNAN, natn, ncn, Dis, RMI, 
                              Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, 
                              MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Hoffset,
                              Atoms_Cut1, Total_NumOrbs, MP, Ngrid)
    

    MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2 = Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)
    density_scratch = zeros(Float64, maximum(GridN_Atom), Nspin)
    max_orbitals = maximum(Total_NumOrbs)
    density_matrix_scratch = zeros(Float64, max_orbitals, max_orbitals, Nspin)
    density_orbital_scratch = ntuple(_ -> zeros(Float64, max_orbitals, density_block_size), 3)
    hamiltonian_orbital_scratch = ntuple(_ -> zeros(Float64, max_orbitals, ham_block_size), 3)
    hamiltonian_product_scratch = zeros(Float64, max_orbitals, max_orbitals)


    return UCell(system_grid, 
                 GridN_Atom, GridListAtom, CellListAtom,
                 MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2,
                 density_scratch, density_matrix_scratch, density_orbital_scratch,
                 hamiltonian_orbital_scratch,
                 hamiltonian_product_scratch,
                 Ngrid)
end


@timeit timer "UCell" function UCell(Nspin, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs=nothing, verbosity=1)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    myrank == 0 && println("<UCell>  Setup Grid ...")


    CpyCell, FNAN, natn, ncn, Dis = Get_FNAN(Latvecs, Natom, Gxyz, Atoms_Cut1)

    if myrank == 0
        for atom = 1:Natom
            @printf("\tCpyCell = %d  atom = %3d  FNAN = %3d\n", CpyCell, atom, FNAN[atom])
        end
    end
    MPI.Barrier(comm)

    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)
    Generation_ATV_ijk!(CpyCell, atv_ijk)

    RMI = Get_RMI(Natom, CpyCell, FNAN, natn, ncn)
    if myrank == 0
        system = Check_system(FNAN, ncn, atv_ijk)
        if verbosity >= 1
            println("<Check_System> The system is $system.")
        end
    end
    

    if !isnothing(Total_NumOrbs)
        MP = zeros(Int32, Natom+1)
        Sum = 0
        for atom = 1:Natom
            MP[atom+1] = Sum + Total_NumOrbs[atom]
            Sum += Total_NumOrbs[atom]
        end
    else
        Total_NumOrbs = zeros(Int32, Natom)
        MP = zeros(Int32, Natom+1)
    end

    GridVol = abs(det(Latvecs))/prod(Ngrid)



    GridN_Atom, GridListAtom, CellListAtom = Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin)

    # split element for MPI
    # one dimensionalization Natom, FNAN, natn, ncn
    Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Hoffset = split_system_grid(
        Natom, FNAN, natn, ncn, Total_NumOrbs;
        GridN_Atom, CpyCell, GridListAtom, CellListAtom, atv_ijk)
    


    system_grid = System_Grid(CpyCell, Natom, atom2spe, Latvecs,
                              atv, atv_ijk, Gxyz, GridVol, Grid_Origin,
                              FNAN, natn, ncn, Dis, RMI, 
                              Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, 
                              MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Hoffset,
                              Atoms_Cut1, Total_NumOrbs, MP, Ngrid)
    

    MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2 = Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)
    density_scratch = zeros(Float64, maximum(GridN_Atom), Nspin)
    max_orbitals = maximum(Total_NumOrbs)
    density_matrix_scratch = zeros(Float64, max_orbitals, max_orbitals, Nspin)
    density_orbital_scratch = ntuple(_ -> zeros(Float64, max_orbitals, density_block_size), 3)
    hamiltonian_orbital_scratch = ntuple(_ -> zeros(Float64, max_orbitals, ham_block_size), 3)
    hamiltonian_product_scratch = zeros(Float64, max_orbitals, max_orbitals)


    return UCell(system_grid, 
                 GridN_Atom, GridListAtom, CellListAtom,
                 MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2,
                 density_scratch, density_matrix_scratch, density_orbital_scratch,
                 hamiltonian_orbital_scratch,
                 hamiltonian_product_scratch,
                 Ngrid)
end


function Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin)
    
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int32, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    Generation_ATV!(CpyCell, Latvecs, atv, atv_ijk, ratv)


    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    # reciprocal grid 
    gRecvecs = 2*pi*inv(gLatvecs')



    # if system is non-periodic when CellListAtom nouse
    GridN_Atom = zeros(Int32, Natom)
    GridListAtom = Vector{Vector{Int32}}(undef, Natom)
    CellListAtom = Vector{Vector{Int32}}(undef, Natom)


    kindex = [[2,3],[3,1],[1,2]]
    nmin = zeros(Int32, 3)
    nmax = zeros(Int32, 3)
    Cxyz = zeros(Float64, 3)
    NOC = zeros(Int32, 4)

    for atom = 1:Natom

        Gx, Gy, Gz = Gxyz[atom]
        rcut = Atoms_Cut1[atom] + 0.5

        for k = 1:3
            i, j = kindex[k]
            vx = Latvecs[i,2]*Latvecs[j,3] - Latvecs[i,3]*Latvecs[j,2]
            vy = Latvecs[i,3]*Latvecs[j,1] - Latvecs[i,1]*Latvecs[j,3]
            vz = Latvecs[i,1]*Latvecs[j,2] - Latvecs[i,2]*Latvecs[j,1]
            coef = inv(sqrt(vx*vx + vy*vy + vz*vz))
            vx = vx*coef
            vy = vy*coef
            vz = vz*coef

            Cx = Gx + rcut*vx - Grid_Origin[1]
            Cy = Gy + rcut*vy - Grid_Origin[2]
            Cz = Gz + rcut*vz - Grid_Origin[3]
            nmax[k] = trunc(Int32, (Cx*gRecvecs[k,1] + Cy*gRecvecs[k,2] + Cz*gRecvecs[k,3])*0.5/pi)

            Cx = Gx - rcut*vx - Grid_Origin[1]
            Cy = Gy - rcut*vy - Grid_Origin[2]
            Cz = Gz - rcut*vz - Grid_Origin[3]
            nmin[k] = trunc(Int32, (Cx*gRecvecs[k,1] + Cy*gRecvecs[k,2] + Cz*gRecvecs[k,3])*0.5/pi)

            if nmax[k] < nmin[k]
                nmin[k], nmax[k] = nmax[k], nmin[k]
            end
        end


        Np = floor(Int32, prod(nmax.-nmin.+1)*3/2)

        Nct = 0
        rcut = Atoms_Cut1[atom]

        tmp_GridListAtom = zeros(Int32, Np)
        tmp_CellListAtom = zeros(Int32, Np)

        for n1 = nmin[1]:nmax[1], n2 = nmin[2]:nmax[2], n3 = nmin[3]:nmax[3]

            Find_CGrids!(NOC, Cxyz, CpyCell, Ngrid, n1, n2, n3, atv, ratv, gLatvecs, Grid_Origin)

            Rn = NOC[1]
            l1 = NOC[2]
            l2 = NOC[3]
            l3 = NOC[4]
            N = l1*Ngrid2*Ngrid3 + l2*Ngrid3 + l3

            dx = Cxyz[1] - Gx
            dy = Cxyz[2] - Gy
            dz = Cxyz[3] - Gz
            R = sqrt(dx^2 + dy^2 + dz^2)

            if R <= rcut
                Nct = Nct + 1
                tmp_GridListAtom[Nct] = N
                tmp_CellListAtom[Nct] = Rn
            end
        end

        GridN_Atom[atom] = Nct
        GridListAtom[atom] = Vector{Int32}(undef, Nct)
        CellListAtom[atom] = Vector{Int32}(undef, Nct)
        for Nc = 1:Nct
            GridListAtom[atom][Nc] = tmp_GridListAtom[Nc]
            CellListAtom[atom][Nc] = tmp_CellListAtom[Nc]
        end
    end

    for atom = 1:Natom
        Grid_sort = sortperm(GridListAtom[atom])
        @. GridListAtom[atom] = GridListAtom[atom][Grid_sort]
        @. CellListAtom[atom] = CellListAtom[atom][Grid_sort]
    end

    
    return GridN_Atom, GridListAtom, CellListAtom
end



function Count_AtomOverlap_Grid!(overlap_counts, CpyCell, pair_atoms,
                                 pair_neighbors, pair_cells, GridN_Atom,
                                 GridListAtom, CellListAtom, atv_ijk)
    ratv = zeros(Int64, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    Generation_RATV!(CpyCell, ratv)

    fill!(overlap_counts, 0)
    for loop in eachindex(pair_atoms)
        atom = pair_atoms[loop]
        jatom = pair_neighbors[loop]
        cell = pair_cells[loop]
        l1, l2, l3 = atv_ijk[cell+1]

        overlap_count = 0
        Nc = 0
        for Nh = 1:GridN_Atom[jatom]
            GNh = GridListAtom[jatom][Nh]
            GRh = CellListAtom[jatom][Nh]
            ll1, ll2, ll3 = atv_ijk[GRh+1]
            lll1 = l1 + ll1
            lll2 = l2 + ll2
            lll3 = l3 + ll3

            if GridListAtom[atom][1] <= GNh
                if GNh == 0
                    Nc = 0
                else
                    while Nc != 0 && GNh <= GridListAtom[atom][Nc+1]
                        Nc = max(Nc - 10, 0)
                    end
                end

                if abs(lll1) <= CpyCell && abs(lll2) <= CpyCell &&
                   abs(lll3) <= CpyCell
                    GRh1 = ratv[lll1+CpyCell+1,
                                 lll2+CpyCell+1,
                                 lll3+CpyCell+1]
                    found = false
                    while !found && Nc < GridN_Atom[atom]
                        GNc = GridListAtom[atom][Nc+1]
                        GRc = CellListAtom[atom][Nc+1]
                        if GNc == GNh && GRc == GRh1
                            overlap_count += 1
                            found = true
                        elseif GNh < GNc
                            found = true
                        end
                        Nc += 1
                    end
                    Nc = max(Nc - 1, 0)
                end
            end
        end
        overlap_counts[loop] = overlap_count
    end

    return overlap_counts
end


"""
    atomoverlap_grid = AtomOverlap_Grid(...)

Create an instance of `AtomOverlap_Grid`.

Mandatory arguments:

- `system_grid`: an instance of `System_Grid`
- `atoms_grid`: an instance of `Atoms_Grid`
"""
function Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)
    
    ratv = zeros(Int64, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    Generation_RATV!(CpyCell, ratv)

    # Find overlap grids between two orbitals
    # GListTAtoms0, GListTAtoms1, GListTAtoms2
    MPI_NumOLG = zeros(Int32, MPI_size)
    MPI_GListTAtoms1 = Vector{Vector{Int32}}(undef, MPI_size)
    MPI_GListTAtoms2 = Vector{Vector{Int32}}(undef, MPI_size)
    
    TAtoms1 = zeros(Int32, maximum(GridN_Atom))
    TAtoms2 = zeros(Int32, maximum(GridN_Atom))

    for loop = 1:MPI_size

        fill!(TAtoms1, 0.0)
        fill!(TAtoms2, 0.0)

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        cell = MPI_ncn[loop]
        l1, l2, l3 = atv_ijk[cell+1]
        
        Nog = -1
        Nc = 0

        for Nh = 1:GridN_Atom[jatom]

            GNh = GridListAtom[jatom][Nh]
            GRh = CellListAtom[jatom][Nh]
            ll1, ll2, ll3 = atv_ijk[GRh+1]
            lll1 = l1 + ll1
            lll2 = l2 + ll2
            lll3 = l3 + ll3
                    
            if GridListAtom[atom][1] <= GNh
                if GNh == 0
                    Nc = 0
                else
                    while GNh <= GridListAtom[atom][Nc+1] && Nc ≠ 0
                        Nc = Nc - 10
                        if Nc < 0
                            Nc = 0
                        end
                    end
                end

                # find whether there is the overlapping or not
                if abs(lll1)<=CpyCell && abs(lll2)<=CpyCell && abs(lll3)<=CpyCell
                            
                    GRh1 = ratv[lll1+CpyCell+1,lll2+CpyCell+1,lll3+CpyCell+1]

                    po = 0

                    while po == 0 && Nc<GridN_Atom[atom]

                        GNc = GridListAtom[atom][Nc+1]
                        GRc = CellListAtom[atom][Nc+1]

                        if GNc==GNh && GRc==GRh1
                            Nog += 1
                            TAtoms1[Nog+1] = Nc
                            TAtoms2[Nog+1] = Nh-1
                            po = 1
                        elseif GNh < GNc
                            po = 1
                        end
                        Nc += 1
                    end
                    Nc -= 1
                    if Nc < 0
                         Nc = 0
                    end
                end
            end

            MPI_NumOLG[loop] = Nog + 1
        end

        MPI_GListTAtoms1[loop] = Vector{Int32}(undef, MPI_NumOLG[loop])
        MPI_GListTAtoms2[loop] = Vector{Int32}(undef, MPI_NumOLG[loop])

        for Nog = 1:MPI_NumOLG[loop]
            MPI_GListTAtoms1[loop][Nog] = TAtoms1[Nog]
            MPI_GListTAtoms2[loop][Nog] = TAtoms2[Nog]
        end
    end
    

    return MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2
end


# get the FNAN, natn, ncn, Dis
function Get_FNAN(Latvecs, Natom, Gxyz, Atoms_Cut1)
    
    po = 0
    CpyCell = 0
    TFNAN = 0
    count = 1

    countmax = 7

    # calc CpyCell
    while po == 0 && count < countmax

        CpyCell = CpyCell + 1
        atv, TCpyCell = Set_Periodic(Latvecs, CpyCell)

        TFNAN_temp = TFNAN

        FNAN = Estimate_Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell)
        TFNAN = sum(FNAN)

        if TFNAN == TFNAN_temp
            po = 1
        end

        count += 1
    end

    if count == countmax || po == 0
        println("count = $(count), po = $(po)")
        error("please check Get_FNAN")
    end

    atv, TCpyCell = Set_Periodic(Latvecs, CpyCell)
    FNAN, natn, ncn, Dis = Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell)
    
    return CpyCell, FNAN, natn, ncn, Dis
end


function Set_Periodic(Latvecs, CpyCell)
    
    TCpyCell = (2*CpyCell + 1)^3 - 1

    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)

    return atv, TCpyCell
end


function Estimate_Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell)
    
    FNAN = zeros(Int32, Natom)

    for atom = 1:Natom
        Gx1, Gy1, Gz1 = Gxyz[atom]
        rcutA = Atoms_Cut1[atom]
        FNAN[atom] = 0

        for jatom = 1:Natom
            Gx2, Gy2, Gz2 = Gxyz[jatom]
            rcutB = Atoms_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell
                if atom == jatom && iszero(Rn)
                    continue
                else
                    dx = abs(Gx1 - Gx2 - atv[Rn+1][1])
                    dy = abs(Gy1 - Gy2 - atv[Rn+1][2])
                    dz = abs(Gz1 - Gz2 - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut
                        r = sqrt(dx^2 + dy^2 + dz^2)

                        if r <= rcut
                            FNAN[atom] = FNAN[atom] + 1
                        end
                    end
                end
            end
        end
    end

    return FNAN
end



function Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell)

    Max_FNAN = maximum(Estimate_Trn_System(Natom, Gxyz, Atoms_Cut1, atv, TCpyCell))
    
    FNAN = zeros(Int32, Natom)

    natn_atom = zeros(Int32, Max_FNAN+1)
    ncn_atom = zeros(Int32, Max_FNAN+1)
    Dis_atom = zeros(Float64, Max_FNAN+1)

    natn = Vector{Vector{Int32}}(undef, Natom)
    ncn = Vector{Vector{Int32}}(undef, Natom)
    Dis = Vector{Vector{Float64}}(undef, Natom)


    for atom = 1:Natom
        Gx1, Gy1, Gz1 = Gxyz[atom]
        FNAN[atom] = 0
        rcutA = Atoms_Cut1[atom]

        for jatom = 1:Natom
            Gx2, Gy2, Gz2 = Gxyz[jatom]
            rcutB = Atoms_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell
                if atom == jatom && iszero(Rn)
                    natn_atom[1] = atom
                    ncn_atom[1] = 0
                    Dis_atom[1] = 0.0
                else
                    dx = abs(Gx1 - Gx2 - atv[Rn+1][1])
                    dy = abs(Gy1 - Gy2 - atv[Rn+1][2])
                    dz = abs(Gz1 - Gz2 - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut
                        r = sqrt(dx^2 + dy^2 + dz^2)
                        if r <= rcut
                            FNAN[atom] = FNAN[atom] + 1
                            natn_atom[FNAN[atom]+1] = jatom
                            ncn_atom[FNAN[atom]+1] = Rn
                            Dis_atom[FNAN[atom]+1] = r
                        end
                    end
                end
            end
        end

        natn[atom] = zeros(Int32, FNAN[atom]+1)
        ncn[atom] = zeros(Int32, FNAN[atom]+1)
        Dis[atom] = zeros(Float64, FNAN[atom]+1)
        for k = 1:FNAN[atom]+1
            natn[atom][k] = natn_atom[k]
            ncn[atom][k] = ncn_atom[k]
            Dis[atom][k] = Dis_atom[k]
        end
    end


    return FNAN, natn, ncn, Dis
end



function Get_RMI(Natom, CpyCell, FNAN, natn, ncn)

    RMI = Vector{Vector{Vector{Int32}}}(undef, Natom)
    for atom = 1:Natom
        RMI[atom] = Vector{Vector{Int32}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            RMI[atom][Rn] = zeros(Int32, FNAN[atom]+1)
        end
    end


    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int64, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)

    Generation_ATV_ijk!(CpyCell, atv_ijk)
    Generation_RATV!(CpyCell, ratv)



    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        ig = natn[atom][Rn]
        Rni = ncn[atom][Rn]+1
        for Rm = 1:FNAN[atom]+1

            jg = natn[atom][Rm]
            Rnj = ncn[atom][Rm]+1
            l1 = atv_ijk[Rnj][1] - atv_ijk[Rni][1]
            l2 = atv_ijk[Rnj][2] - atv_ijk[Rni][2]
            l3 = atv_ijk[Rnj][3] - atv_ijk[Rni][3]
            m1 = ifelse(l1<0, -l1, l1)
            m2 = ifelse(l2<0, -l2, l2)
            m3 = ifelse(l3<0, -l3, l3)  

            if m1 <= CpyCell && m2 <= CpyCell && m3 <= CpyCell  
                cell = ratv[l1+CpyCell+1,l2+CpyCell+1,l3+CpyCell+1] 
                k = 0
                po = 0
                RMI[atom][Rn][Rm] = -1  

                while po == 0 && k <= FNAN[ig]
                    if natn[ig][k+1]==jg && ncn[ig][k+1] == cell
                        RMI[atom][Rn][Rm] = k
                        po = 1
                    end
                    k = k + 1
                end
            else
                RMI[atom][Rn][Rm] = -1
            end
        end
    end
    

    return RMI
end


function Calc_Dis(Natom, FNAN, Gxyz, atv, Atoms_Cut1, CpyCell)

    TCpyCell = (2*CpyCell+1)^3-1
    Max_FNAN = maximum(FNAN)
    Dis_atom = zeros(Float64, Max_FNAN+1)


    Dis = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Dis[atom] = zeros(Float64, FNAN[atom]+1)
    end

    for atom = 1:Natom
        tmp = 0
        Gx1, Gy1, Gz1 = Gxyz[atom]
        rcutA = Atoms_Cut1[atom]
        for jatom = 1:Natom
            Gx2, Gy2, Gz2 = Gxyz[jatom]
            rcutB = Atoms_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell
                if atom == jatom && iszero(Rn)
                    Dis_atom[1] = 0.0
                else
                    dx = abs(Gx1 - Gx2 - atv[Rn+1][1])
                    dy = abs(Gy1 - Gy2 - atv[Rn+1][2])
                    dz = abs(Gz1 - Gz2 - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut
                        r = sqrt(dx^2 + dy^2 + dz^2)
                        if r <= rcut
                            tmp = tmp + 1
                            Dis_atom[tmp+1] = r
                        end
                    end
                end
            end
        end

        for Rn = 1:FNAN[atom]+1
            Dis[atom][Rn] = Dis_atom[Rn]
        end
    end


    return Dis
end


function Check_system(FNAN, ncn, atv_ijk)

    Natom = length(FNAN)
    po = zeros(Int64, 3)

    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        cell = ncn[atom][Rn]
        if cell ≠ 0
            if atv_ijk[cell+1][1] ≠ 0
                po[1] = 1
            elseif atv_ijk[cell+1][2] ≠ 0
                po[2] = 1
            elseif atv_ijk[cell+1][3] ≠ 0
                po[3] = 1
            end
        end
    end

    num = sum(po)

    if num == 0
        system = "molecule"
    elseif num == 1
        system = "chain"
    elseif num == 2
        system = "slab"
    elseif num == 3
        system = "bluk"
    else
        error("please check Check_System")
    end

    return system
end
