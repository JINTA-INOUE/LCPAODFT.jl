struct MPIDensityPlan
    comm::MPI.Comm
    Ngrid::NTuple{3,Int}
    global_indices::Vector{Int32}
    offsets::Vector{Int32}
    B::Array{Float64,3}
end


# DFT_Setup stores raw Spglib operations and a Cartesian grid origin.
# Convert them to the low-level convention x′ = transpose(S)*x - f.
function MPIDensityPlan(Latvecs, Grid_Origin, Ngrid, symmetry; global_time_reversal::Bool=false)
    
    comm = MPI.COMM_WORLD
    SpinPol = symmetry.SpinPol
    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end

    lattice = transpose(Latvecs)
    Grid_Origin_frac = lattice \ Grid_Origin
    rotations = [Matrix(transpose(W)) for W in symmetry.rotations]
    translations = [-t for t in symmetry.translations]
    cartesian_rotations = Nspin == 4 ? [lattice*W/lattice for W in symmetry.rotations] : nothing
    
    return MPIDensityPlan(
        Nspin, Ngrid, rotations, translations, Grid_Origin_frac;
        time_reversals=Bool.(symmetry.time_reversals), cartesian_rotations, global_time_reversal)
end


function MPIDensityPlan(
    Nspin::Integer, Ngrid::NTuple{3,<:Integer}, rotations, translations, Grid_Origin_frac;
    time_reversals=falses(length(rotations)), cartesian_rotations=nothing,
    global_time_reversal::Bool = false)

    comm = MPI.COMM_WORLD
    Ngrid = Int.(Ngrid)
    nops = length(rotations)
    maps,shifts,spin=Matrix{Int}[],Vector{Int}[],Matrix{Float64}[]
    for s = 1:nops
        S, f = rotations[s], translations[s]
        M = [S[i,j]*Ngrid[j]/Ngrid[i] for j=1:3, i=1:3]
        b = (transpose(S)*Grid_Origin_frac - f - Grid_Origin_frac).*collect(Ngrid)
        D = Matrix{Float64}(I, Nspin, Nspin)
        if Nspin == 2
            time_reversals[s] && (D = D[:,[2,1]])
            global_time_reversal && fill!(D, 0.5)
        elseif Nspin == 4
            C = cartesian_rotations[s]
            sign = (time_reversals[s] ? -1.0 : 1.0)*(det(C)<0 ? -1.0 : 1.0)
            D[2:4,2:4] .= sign.*C
            global_time_reversal && (D[2:4,2:4].=0.0)
        end
        push!(maps, round.(Int, M))
        push!(shifts, round.(Int, b))
        push!(spin, D)
    end

    # The smallest global index represents an orbit. Each candidate belongs
    # to one rank; reject immediately if a smaller image exists. This
    # avoids a global visited array and building a full plan on every rank.
    NN = prod(Ngrid)
    rank, nranks = MPI.Comm_rank(comm), MPI.Comm_size(comm)
    indices, offsets = Int[], Int[]
    B = Array{Float64}(undef, Nspin, Nspin, 0)
    targets, counts = zeros(Int, nops), zeros(Int, nops)
    block = zeros(Nspin, Nspin, nops)
    n1,n2,n3 = Ngrid

    # Count first, then allocate exactly once. No oversized growth buffers or
    # temporary copy of the large B array, even during construction.
    for pass = 1:2
        npoints, norbits = 0,0
        for seed = rank+1:nranks:NN
            q = seed-1
            z = mod(q,n3)
            q = div(q,n3)
            y = mod(q,n2)
            x = div(q,n2)
            nlocal = 0
            representative = true
            if pass == 2
                fill!(counts, 0)
                fill!(block, 0.0)
            end

            for s = 1:nops
                M, b = maps[s], shifts[s]
                a = mod(M[1,1]*x + M[1,2]*y + M[1,3]*z + b[1], n1)
                d = mod(M[2,1]*x + M[2,2]*y + M[2,3]*z + b[2], n2)
                e = mod(M[3,1]*x + M[3,2]*y + M[3,3]*z + b[3], n3)
                target = (a*n2 + d)*n3 + e + 1
                if target < seed
                    representative = false
                    break
                end
                k = 1
                while k<=nlocal && targets[k] ≠ target
                    k += 1
                end
                if k > nlocal
                    nlocal += 1
                    targets[k] = target
                end
                if pass == 2
                    counts[k] += 1
                    for j = 1:Nspin,i = 1:Nspin
                        block[i,j,k] += spin[s][i,j]
                    end
                end
            end
            representative || continue
            norbits += 1
            if pass == 2
                offsets[norbits] = npoints + 1
                for k = 1:nlocal
                    indices[npoints+k] = targets[k]
                    for j = 1:Nspin, i = 1:Nspin
                        B[i,j,npoints+k] = block[i,j,k]/counts[k]
                    end
                end
            end
            npoints += nlocal
        end
        if pass == 1
            indices = Vector{Int}(undef, npoints)
            offsets = Vector{Int}(undef, norbits+1)
            B = Array{Float64}(undef, Nspin, Nspin, npoints)
        else
            offsets[end] = npoints+1
        end
    end

    return MPIDensityPlan(comm, Ngrid, indices, offsets, B)
end


# A single vector is the full grid in the original ADensity_Grid order.
# Use a scalar (Nspin=1) plan and call on every rank with the same input.
function symmetrize_Grid!(plan::MPIDensityPlan, Density_Grid::Vector{Float64})
    indices = plan.global_indices
    @inbounds for orbit = 1:length(plan.offsets)-1
        lo, hi = plan.offsets[orbit], plan.offsets[orbit+1]-1
        value = 0.0
        for k = lo:hi
            value += Density_Grid[indices[k]]
        end
        value /= hi-lo+1
        for k = lo:hi
            Density_Grid[indices[k]] = value
        end
    end

    # Share completed orbits without another full-grid density buffer.
    comm = plan.comm
    nranks = MPI.Comm_size(comm)
    if nranks > 1
        rank = MPI.Comm_rank(comm)
        counts = MPI.Allgather(length(indices), comm)
        ids = similar(indices, maximum(counts))
        values = Vector{Float64}(undef, maximum(counts))
        for owner = 0:nranks-1
            n = counts[owner+1]
            if rank == owner
                for k = 1:n
                    ids[k] = indices[k]
                    values[k] = Density_Grid[indices[k]]
                end
            end
            MPI.Bcast!(view(ids, 1:n), owner, comm)
            MPI.Bcast!(view(values, 1:n), owner, comm)
            for k = 1:n
                Density_Grid[ids[k]] = values[k]
            end
        end
    end
end


# Full-grid channels: scatter, symmetrize, gather into the supplied arrays,
# then broadcast the result. Call collectively on every rank; root is the input.
# Set local_grid=true only for already distributed channels from scatter_density.
function symmetrize_Grid!(plan::MPIDensityPlan, Density_Grid::Vector{Vector{Float64}};
    root::Int=0, local_grid::Bool=false)
    local_density = local_grid ? Density_Grid : scatter_density(plan, Density_Grid; root)
    Nspin = size(plan.B, 1)
    if Nspin == 1
        _symmetrize_mpi_channels!(plan, (local_density[1],))
    elseif Nspin == 2
        _symmetrize_mpi_channels!(plan, (local_density[1], local_density[2]))
    else
        _symmetrize_mpi_channels!(plan, (local_density[1], local_density[2], local_density[3], local_density[4]))
    end
    if !local_grid
        gather_density(plan, local_density; root, destination=Density_Grid)
        for density in Density_Grid
            MPI.Bcast!(density, root, plan.comm)
        end
    end
    return Density_Grid
end


function _symmetrize_mpi_channels!(plan::MPIDensityPlan, Density_Grid::NTuple{Nspin,Vector{Float64}}) where {Nspin}
    B = plan.B
    @inbounds for orbit = 1:length(plan.offsets)-1
        lo, hi = plan.offsets[orbit], plan.offsets[orbit+1]-1
        v = ntuple(Val(Nspin)) do spin2
            value = 0.0
            for k = lo:hi, spin1 = 1:Nspin
                value += B[spin1,spin2,k]*Density_Grid[spin1][k]
            end
            value/(hi-lo+1)
        end
        for k = lo:hi, spin1 = 1:Nspin
            value = 0.0
            for spin2 = 1:Nspin
                value += B[spin1,spin2,k]*v[spin2]
            end
            Density_Grid[spin1][k] = value
        end
    end
    return Density_Grid
end


function symmetrize_noncollinear_matrix!(plan::MPIDensityPlan, channels;
    root::Int=0, local_grid::Bool=false)
    a,b,c,d = channels
    for i in eachindex(a)
        a[i],b[i],c[i],d[i] = a[i]+b[i],2c[i],-2d[i],a[i]-b[i]
    end
    symmetrize_Grid!(plan, channels; root, local_grid)
    for i in eachindex(a)
        a[i],b[i],c[i],d[i] = (a[i]+d[i])/2,(a[i]-d[i])/2,b[i]/2,-c[i]/2
    end
end

# Transfers are collective: all ranks call in the same order on an otherwise
# idle communicator. Root streams one rank at a time using one index buffer
# and one single-channel buffer. Never gather all indices or pack N*nc values.
function _transfer_mpi_density(plan::MPIDensityPlan,channels,root,scatter; destination=nothing)

    comm = plan.comm
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    Nspin = size(plan.B, 1)
    nlocal = length(plan.global_indices)
    counts = MPI.Allgather(nlocal,comm)
    localdata = scatter ? [Vector{Float64}(undef,nlocal) for _=1:Nspin] : channels
    globaldata = if rank == root
        scatter ? channels : (destination === nothing ?
            [Vector{Float64}(undef,prod(plan.Ngrid)) for _=1:Nspin] : destination)
    else
        nothing
    end

    if rank == root
        indices = similar(plan.global_indices, maximum(counts))
        buffer = Vector{Float64}(undef, maximum(counts))
        for peer = 0:nranks-1
            if peer == root
                for spin = 1:Nspin, k = 1:nlocal
                    g = plan.global_indices[k]
                    if scatter
                        localdata[spin][k] = globaldata[spin][g]
                    else
                        globaldata[spin][g] = localdata[spin][k]
                    end
                end
            else
                ids = view(indices, 1:counts[peer+1])
                buf = view(buffer, 1:counts[peer+1])
                MPI.Recv!(ids,comm; source=peer, tag=1701)
                for spin = 1:Nspin
                    if scatter
                        for k in eachindex(ids);buf[k]=globaldata[spin][ids[k]];end
                        MPI.Send(buf,comm; dest=peer, tag=1702)
                    else
                        MPI.Recv!(buf,comm; source=peer, tag=1702)
                        for k in eachindex(ids);globaldata[spin][ids[k]]=buf[k];end
                    end
                end
            end
        end
    else
        MPI.Send(plan.global_indices, comm; dest=root, tag=1701)
        for spin = 1:Nspin
            if scatter
                MPI.Recv!(localdata[spin], comm; source=root, tag=1702)
            else
                MPI.Send(localdata[spin], comm; dest=root, tag=1702)
            end
        end
    end
    scatter ? localdata : globaldata
end

scatter_density(plan::MPIDensityPlan, channels=nothing; root::Int=0) = _transfer_mpi_density(plan,channels,root,true)
gather_density(plan::MPIDensityPlan, channels; root::Int=0, destination=nothing) = _transfer_mpi_density(plan,channels,root,false; destination)

