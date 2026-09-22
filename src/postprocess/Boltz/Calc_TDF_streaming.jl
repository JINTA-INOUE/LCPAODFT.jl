"""Reduced orbital transport moments used when a full decomposed TDF is not requested."""
struct BoltzDecompMoments
    L0::Array{Float64,5}       # orbital, spin, component, temperature, mu
    L1::Array{Float64,5}       # orbital, spin, component, temperature, mu
    total_L0::Array{Float64,4} # spin, component, temperature, mu
end


function _boltz_tdf_energy_grid(boltz_setup::Boltz_Setup)
    TDF_Emin, TDF_Emax = boltz_setup.TDF_Erange
    TDF_EneNum = floor(Int,
        (TDF_Emax - TDF_Emin) / boltz_setup.TDF_dE)
    TDF_EneNum >= 2 ||
        error("the TDF energy grid must contain at least two points")

    TDF_Energy = Vector{Float64}(undef, TDF_EneNum)
    for ie = 1:TDF_EneNum
        TDF_Energy[ie] = TDF_Emin +
            (TDF_Emax - TDF_Emin) * (ie - 1) / (TDF_EneNum - 1)
    end
    return TDF_Energy
end


function _boltz_decomp_kernels(boltz_setup::Boltz_Setup, TDF_Energy)
    Nenergy = length(TDF_Energy)
    Ntemp = length(boltz_setup.Temp)
    Nmu = length(boltz_setup.muE)
    kernel0 = Array{Float64}(undef, Nenergy, Ntemp, Nmu)
    kernel1 = similar(kernel0)
    fermi_T_dE = Vector{Float64}(undef, Nenergy)
    dTDF = (boltz_setup.TDF_Erange[2] - boltz_setup.TDF_Erange[1]) /
           (Nenergy - 1)

    for itemp = 1:Ntemp, imu = 1:Nmu
        temperature = boltz_setup.Temp[itemp]
        mu = boltz_setup.muE[imu]
        kBT = temperature * k_B_SI / elem_charge_SI
        Set_dfdE!(fermi_T_dE, mu, TDF_Energy, Nenergy, kBT)
        @inbounds for ie = 1:Nenergy
            kernel0[ie, itemp, imu] = fermi_T_dE[ie] * dTDF
            kernel1[ie, itemp, imu] =
                fermi_T_dE[ie] * (TDF_Energy[ie] - mu) *
                dTDF / temperature
        end
    end
    return kernel0, kernel1
end


function _boltz_accumulate_tdf_block!(
    local_TDF, boltz_setup::Boltz_Setup, TDF_Energy,
    Enk, Vnk, vertex_positions)

    material = boltz_setup.material
    spinsize = material.SpinPol == "on" ? 2 : 1
    Nstate = Set_Nstate(material)
    plane_type = boltz_setup.plane_type
    ncomponents = plane_type ? 3 : 6
    velocity_pairs = plane_type ?
        ((1, 1), (1, 2), (2, 2)) :
        ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))
    TDF_Emin, TDF_Emax = boltz_setup.TDF_Erange
    TDF_EneNum = length(TDF_Energy)

    cell_e = Vector{Float64}(undef, 8)
    cell_a = Matrix{Float64}(undef, 8, ncomponents)
    tetra_e = Vector{Float64}(undef, 4)
    tetra_vertex = Vector{Int}(undef, 4)
    tetra_id = ((1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
                (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8))

    for spin = 1:spinsize, band = 1:Nstate
        for cell_position = axes(vertex_positions, 2)
            @inbounds for vertex = 1:8
                k_position = vertex_positions[vertex, cell_position]
                cell_e[vertex] = Enk[band, spin, k_position]
                for component = 1:ncomponents
                    xyz1, xyz2 = velocity_pairs[component]
                    cell_a[vertex, component] =
                        Vnk[band, xyz1, spin, k_position] *
                        Vnk[band, xyz2, spin, k_position]
                end
            end

            for vertices in tetra_id
                @inbounds for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                    tetra_vertex[ic] = vertices[ic]
                end
                OrderE!(tetra_e, tetra_vertex, 4)

                first_energy = trunc(Int,
                    (tetra_e[1] - TDF_Emin) /
                    (TDF_Emax - TDF_Emin) * (TDF_EneNum - 1) - 1)
                last_energy = trunc(Int,
                    (tetra_e[4] - TDF_Emin) /
                    (TDF_Emax - TDF_Emin) * (TDF_EneNum - 1) + 1)
                first_energy = max(first_energy, 0)
                last_energy = min(last_energy, TDF_EneNum - 1)

                if 0 <= first_energy < TDF_EneNum &&
                   0 <= last_energy < TDF_EneNum &&
                   first_energy <= last_energy
                    for ie = first_energy:last_energy
                        w1, w2, w3, w4 = ATM_Spectrum_weights(
                            tetra_e, TDF_Energy[ie + 1])
                        v1, v2, v3, v4 = tetra_vertex
                        @inbounds for component = 1:ncomponents
                            local_TDF[ie + 1, spin, component] +=
                                w1 * cell_a[v1, component] +
                                w2 * cell_a[v2, component] +
                                w3 * cell_a[v3, component] +
                                w4 * cell_a[v4, component]
                        end
                    end
                end
            end
        end
    end
    return nothing
end


function _boltz_accumulate_tdf_decomp_block!(
    local_TDF_decomp, boltz_setup::Boltz_Setup, TDF_Energy,
    Enk, EVec, Vnk, vertex_positions)

    material = boltz_setup.material
    spinsize = material.SpinPol == "on" ? 2 : 1
    Nwann = Int(material.Ngsize)
    plane_type = boltz_setup.plane_type
    ncomponents = plane_type ? 3 : 6
    velocity_pairs = plane_type ?
        ((1, 1), (1, 2), (2, 2)) :
        ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))
    TDF_Emin, TDF_Emax = boltz_setup.TDF_Erange
    TDF_EneNum = length(TDF_Energy)

    cell_e = Vector{Float64}(undef, 8)
    cell_a = Matrix{Float64}(undef, 8, ncomponents)
    tetra_e = Vector{Float64}(undef, 4)
    tetra_vertex = Vector{Int}(undef, 4)
    weighted_vertex = Matrix{Float64}(undef, 4, Nwann)
    tetra_id = ((1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
                (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8))

    for spin = 1:spinsize, band = 1:Nwann
        for cell_position = axes(vertex_positions, 2)
            @inbounds for vertex = 1:8
                k_position = vertex_positions[vertex, cell_position]
                cell_e[vertex] = Enk[band, spin, k_position]
                for component = 1:ncomponents
                    xyz1, xyz2 = velocity_pairs[component]
                    cell_a[vertex, component] =
                        Vnk[band, xyz1, spin, k_position] *
                        Vnk[band, xyz2, spin, k_position]
                end
            end

            for vertices in tetra_id
                @inbounds for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                    tetra_vertex[ic] = vertices[ic]
                end
                OrderE!(tetra_e, tetra_vertex, 4)

                first_energy = floor(Int,
                    (tetra_e[1] - TDF_Emin) /
                    (TDF_Emax - TDF_Emin) * (TDF_EneNum - 1) - 1)
                last_energy = floor(Int,
                    (tetra_e[4] - TDF_Emin) /
                    (TDF_Emax - TDF_Emin) * (TDF_EneNum - 1) + 1)
                first_energy = max(first_energy, 0)
                last_energy = min(last_energy, TDF_EneNum - 1)

                if 0 < first_energy < TDF_EneNum &&
                   0 <= last_energy < TDF_EneNum &&
                   first_energy <= last_energy
                    v1, v2, v3, v4 = tetra_vertex
                    p1 = vertex_positions[v1, cell_position]
                    p2 = vertex_positions[v2, cell_position]
                    p3 = vertex_positions[v3, cell_position]
                    p4 = vertex_positions[v4, cell_position]
                    for ie = first_energy:last_energy
                        w1, w2, w3, w4 = ATM_Spectrum_weights(
                            tetra_e, TDF_Energy[ie + 1])
                        @inbounds for orbital = 1:Nwann
                            weighted_vertex[1, orbital] = w1 *
                                EVec[orbital, band, spin, p1]
                            weighted_vertex[2, orbital] = w2 *
                                EVec[orbital, band, spin, p2]
                            weighted_vertex[3, orbital] = w3 *
                                EVec[orbital, band, spin, p3]
                            weighted_vertex[4, orbital] = w4 *
                                EVec[orbital, band, spin, p4]
                        end
                        @inbounds for component = 1:ncomponents,
                                      orbital = 1:Nwann
                            local_TDF_decomp[orbital][
                                ie, spin, component] +=
                                weighted_vertex[1, orbital] *
                                    cell_a[v1, component] +
                                weighted_vertex[2, orbital] *
                                    cell_a[v2, component] +
                                weighted_vertex[3, orbital] *
                                    cell_a[v3, component] +
                                weighted_vertex[4, orbital] *
                                    cell_a[v4, component]
                        end
                    end
                end
            end
        end
    end
    return nothing
end


"""
Accumulate the ordinary TDF and decomposed transport moments in one cell
traversal.  The tetrahedron energy ordering and interpolation weights are
shared by all tensor components and orbitals.
"""
function _boltz_accumulate_tdf_and_decomp_moments_block!(
    local_TDF, local_L0, local_L1, kernel0, kernel1,
    boltz_setup::Boltz_Setup, TDF_Energy,
    Enk, EVec, Vnk, vertex_positions)

    material = boltz_setup.material
    spinsize = material.SpinPol == "on" ? 2 : 1
    Nwann = Int(material.Ngsize)
    ncomponents = boltz_setup.plane_type ? 3 : 6
    velocity_pairs = boltz_setup.plane_type ?
        ((1, 1), (1, 2), (2, 2)) :
        ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))
    TDF_Emin, TDF_Emax = boltz_setup.TDF_Erange
    TDF_EneNum = length(TDF_Energy)
    energy_scale = (TDF_EneNum - 1) / (TDF_Emax - TDF_Emin)
    Ntemp = size(kernel0, 2)
    Nmu = size(kernel0, 3)

    cell_e = Vector{Float64}(undef, 8)
    cell_a = Matrix{Float64}(undef, 8, ncomponents)
    tetra_e = Vector{Float64}(undef, 4)
    tetra_vertex = Vector{Int}(undef, 4)
    weighted_vertex = Matrix{Float64}(undef, 4, Nwann)
    tetra_id = ((1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
                (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8))

    for spin = 1:spinsize, band = 1:Nwann
        for cell_position = axes(vertex_positions, 2)
            @inbounds for vertex = 1:8
                k_position = vertex_positions[vertex, cell_position]
                cell_e[vertex] = Enk[band, spin, k_position]
                for component = 1:ncomponents
                    xyz1, xyz2 = velocity_pairs[component]
                    cell_a[vertex, component] =
                        Vnk[band, xyz1, spin, k_position] *
                        Vnk[band, xyz2, spin, k_position]
                end
            end

            for vertices in tetra_id
                @inbounds for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                    tetra_vertex[ic] = vertices[ic]
                end
                OrderE!(tetra_e, tetra_vertex, 4)

                ordinary_first = trunc(Int,
                    (tetra_e[1] - TDF_Emin) * energy_scale - 1)
                ordinary_last = trunc(Int,
                    (tetra_e[4] - TDF_Emin) * energy_scale + 1)
                ordinary_first = max(ordinary_first, 0)
                ordinary_last = min(ordinary_last, TDF_EneNum - 1)
                ordinary_valid =
                    0 <= ordinary_first < TDF_EneNum &&
                    0 <= ordinary_last < TDF_EneNum &&
                    ordinary_first <= ordinary_last

                # Preserve the established decomposed-TDF indexing.  Its
                # zero-based energy index is paired with kernel index `ie`.
                decomp_first = floor(Int,
                    (tetra_e[1] - TDF_Emin) * energy_scale - 1)
                decomp_last = floor(Int,
                    (tetra_e[4] - TDF_Emin) * energy_scale + 1)
                decomp_first = max(decomp_first, 0)
                decomp_last = min(decomp_last, TDF_EneNum - 1)
                decomp_valid =
                    0 < decomp_first < TDF_EneNum &&
                    0 <= decomp_last < TDF_EneNum &&
                    decomp_first <= decomp_last

                (ordinary_valid || decomp_valid) || continue
                first_energy = ordinary_valid && decomp_valid ?
                    min(ordinary_first, decomp_first) :
                    (ordinary_valid ? ordinary_first : decomp_first)
                last_energy = ordinary_valid && decomp_valid ?
                    max(ordinary_last, decomp_last) :
                    (ordinary_valid ? ordinary_last : decomp_last)

                v1, v2, v3, v4 = tetra_vertex
                p1 = vertex_positions[v1, cell_position]
                p2 = vertex_positions[v2, cell_position]
                p3 = vertex_positions[v3, cell_position]
                p4 = vertex_positions[v4, cell_position]
                for ie = first_energy:last_energy
                    w1, w2, w3, w4 = ATM_Spectrum_weights(
                        tetra_e, TDF_Energy[ie + 1])

                    if ordinary_valid && ordinary_first <= ie <= ordinary_last
                        @inbounds for component = 1:ncomponents
                            local_TDF[ie + 1, spin, component] +=
                                w1 * cell_a[v1, component] +
                                w2 * cell_a[v2, component] +
                                w3 * cell_a[v3, component] +
                                w4 * cell_a[v4, component]
                        end
                    end

                    if decomp_valid && decomp_first <= ie <= decomp_last
                        @inbounds for orbital = 1:Nwann
                            weighted_vertex[1, orbital] = w1 *
                                EVec[orbital, band, spin, p1]
                            weighted_vertex[2, orbital] = w2 *
                                EVec[orbital, band, spin, p2]
                            weighted_vertex[3, orbital] = w3 *
                                EVec[orbital, band, spin, p3]
                            weighted_vertex[4, orbital] = w4 *
                                EVec[orbital, band, spin, p4]
                        end
                        @inbounds for component = 1:ncomponents,
                                      orbital = 1:Nwann
                            spectrum =
                                weighted_vertex[1, orbital] *
                                    cell_a[v1, component] +
                                weighted_vertex[2, orbital] *
                                    cell_a[v2, component] +
                                weighted_vertex[3, orbital] *
                                    cell_a[v3, component] +
                                weighted_vertex[4, orbital] *
                                    cell_a[v4, component]
                            for itemp = 1:Ntemp, imu = 1:Nmu
                                local_L0[orbital, spin, component,
                                         itemp, imu] +=
                                    spectrum * kernel0[ie, itemp, imu]
                                local_L1[orbital, spin, component,
                                         itemp, imu] +=
                                    spectrum * kernel1[ie, itemp, imu]
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end


"""
Generate the ordinary and optional orbital-decomposed TDF in bounded cell
blocks.  Only the unique vertex k points of the current block are resident.
"""
@timeit timer "Calc_TDF_streaming" function Calc_TDF_streaming(
    boltz_setup::Boltz_Setup, kpoints::BoltzKPoints;
    max_vertices_override::Union{Nothing,Integer}=nothing)

    material = boltz_setup.material
    SpinPol = material.SpinPol
    spinsize = SpinPol == "on" ? 2 : 1
    Nwann = Int(material.Ngsize)
    ncomponents = boltz_setup.plane_type ? 3 : 6
    TDF_Energy = _boltz_tdf_energy_grid(boltz_setup)
    TDF_EneNum = length(TDF_Energy)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    local_TDF = zeros(Float64, TDF_EneNum, spinsize, ncomponents)
    keep_decomp_tdf = boltz_setup.decomp && boltz_setup.Write_TDF
    use_decomp_moments = boltz_setup.decomp && !boltz_setup.Write_TDF
    local_TDF_decomp = keep_decomp_tdf ?
        [zeros(Float64, TDF_EneNum, spinsize, ncomponents)
         for _ = 1:Nwann] : nothing
    Ntemp = length(boltz_setup.Temp)
    Nmu = length(boltz_setup.muE)
    local_L0 = use_decomp_moments ?
        zeros(Float64, Nwann, spinsize, ncomponents, Ntemp, Nmu) : nothing
    local_L1 = use_decomp_moments ? similar(local_L0) : nothing
    if use_decomp_moments
        fill!(local_L1, 0.0)
        kernel0, kernel1 = _boltz_decomp_kernels(
            boltz_setup, TDF_Energy)
    end

    cell_range = kpoints.MPI_krange[myrank + 1]
    max_vertices = isnothing(max_vertices_override) ?
        _boltz_max_block_vertices(boltz_setup) : Int(max_vertices_override)
    max_vertices >= 8 || error("max_vertices_override must be at least 8")
    last_cell = last(cell_range)

    # Determine the exact workspace capacity without retaining all block
    # metadata.  Enk/Vnk/EVec and the Fourier work matrices are then reused
    # for every block, avoiding a full allocation and collection cycle per
    # block while preserving the configured memory bound.
    scan_cell = first(cell_range)
    workspace_capacity = 0
    while scan_cell <= last_cell
        scan_block, scan_cell = _boltz_next_cell_block(
            scan_cell, last_cell, boltz_setup.kmesh, max_vertices)
        workspace_capacity = max(
            workspace_capacity, length(scan_block.global_kpoints))
    end
    eigensystem_workspace = _boltz_eigensystem_workspace(
        material, workspace_capacity, Val(boltz_setup.decomp))

    next_cell = first(cell_range)
    nblocks = 0
    largest_block_vertices = 0

    while next_cell <= last_cell
        block, next_cell = _boltz_next_cell_block(
            next_cell, last_cell, boltz_setup.kmesh, max_vertices)
        coordinates = _boltz_kpoint_coordinates(
            block.global_kpoints, boltz_setup.kmesh)
        Enk, Vnk, EVec = _fill_boltz_eigensystem!(
            eigensystem_workspace, material, coordinates,
            Val(boltz_setup.decomp))

        if use_decomp_moments
            _boltz_accumulate_tdf_and_decomp_moments_block!(
                local_TDF, local_L0, local_L1, kernel0, kernel1,
                boltz_setup, TDF_Energy,
                Enk, EVec, Vnk, block.vertex_positions)
        else
            _boltz_accumulate_tdf_block!(
                local_TDF, boltz_setup, TDF_Energy,
                Enk, Vnk, block.vertex_positions)
        end
        if keep_decomp_tdf
            _boltz_accumulate_tdf_decomp_block!(
                local_TDF_decomp, boltz_setup, TDF_Energy,
                Enk, EVec, Vnk, block.vertex_positions)
        end

        nblocks += 1
        largest_block_vertices =
            max(largest_block_vertices, length(block.global_kpoints))

        # The large k-point arrays remain owned by the reusable workspace.
        # Only the small block metadata and coordinate list become dead here.
        Enk = nothing
        Vnk = nothing
        EVec = nothing
        coordinates = nothing
        block = nothing
    end

    eigensystem_workspace = nothing
    GC.gc()

    max_blocks = MPI.Reduce(nblocks, MPI.MAX, comm; root=0)
    max_stored_vertices = MPI.Reduce(
        largest_block_vertices, MPI.MAX, comm; root=0)
    legacy_k_bytes = kpoints.MPI_Nkpt * spinsize *
        (16 * Nwann * Nwann + 32 * Nwann +
         (boltz_setup.decomp ? 4 * Nwann * Nwann : 0))
    streaming_k_bytes = largest_block_vertices * spinsize *
        (32 * Nwann +
         (boltz_setup.decomp ? 4 * Nwann * Nwann : 0))
    max_legacy_k_bytes = MPI.Reduce(legacy_k_bytes, MPI.MAX, comm; root=0)
    max_streaming_k_bytes = MPI.Reduce(
        streaming_k_bytes, MPI.MAX, comm; root=0)
    fourier_workspace_bytes = Threads.maxthreadid() * 96 * Nwann * Nwann
    max_fourier_workspace_bytes = MPI.Reduce(
        fourier_workspace_bytes, MPI.MAX, comm; root=0)
    if myrank == 0
        println("<Boltz streaming memory>")
        println("\tmaximum blocks per rank: $max_blocks")
        println("\tmaximum resident k-point vertices per rank: " *
                "$max_stored_vertices")
        println("\tconfigured k-block payload: " *
                "$(boltz_setup.kblock_memory_mb) MiB")
        println("\tlegacy k-resolved payload estimate per rank: " *
                "$(round(max_legacy_k_bytes / 2.0^20; digits=3)) MiB")
        println("\tstreamed k-resolved payload maximum per rank: " *
                "$(round(max_streaming_k_bytes / 2.0^20; digits=3)) MiB")
        println("\tFourier thread workspace per rank: " *
                "$(round(max_fourier_workspace_bytes / 2.0^20; digits=3)) MiB")
    end

    MPI.Reduce!(local_TDF, MPI.SUM, comm; root=0)
    TDF = myrank == 0 ? local_TDF : nothing
    local_TDF = nothing

    TDF_decomp = keep_decomp_tdf && myrank == 0 ?
        Vector{Array{Float64,3}}(undef, Nwann) : nothing
    if keep_decomp_tdf
        for orbital = 1:Nwann
            MPI.Reduce!(local_TDF_decomp[orbital], MPI.SUM, comm; root=0)
            if myrank == 0
                TDF_decomp[orbital] = local_TDF_decomp[orbital]
            else
                local_TDF_decomp[orbital] =
                    Array{Float64}(undef, 0, 0, 0)
            end
        end
        local_TDF_decomp = nothing
    end

    if use_decomp_moments
        MPI.Reduce!(local_L0, MPI.SUM, comm; root=0)
        MPI.Reduce!(local_L1, MPI.SUM, comm; root=0)
        reduced_L0 = myrank == 0 ? local_L0 : nothing
        reduced_L1 = myrank == 0 ? local_L1 : nothing
        local_L0 = nothing
        local_L1 = nothing
    end
    GC.gc()

    if myrank == 0
        cell_volume_Ang = abs(det(material.Latvecs)) / Ang_to_bohr^3
        factor = boltz_setup.tau /
                 (prod(boltz_setup.kmesh) * 6 * cell_volume_Ang)
        TDF .*= factor
        if keep_decomp_tdf
            for orbital = 1:Nwann
                TDF_decomp[orbital] .*= factor
            end
        elseif use_decomp_moments
            reduced_L0 .*= factor
            reduced_L1 .*= factor
            # The established decomposed-Seebeck path constructs its
            # conductivity denominator by summing the orbital TDFs.  Use
            # the already reduced orbital moments here as well; integrating
            # the ordinary TDF instead can differ after Float32 orbital
            # weights and becomes visible when the tensor is ill-conditioned.
            total_L0 = zeros(
                Float64, spinsize, ncomponents, Ntemp, Nmu)
            @inbounds for orbital = 1:Nwann, spin = 1:spinsize,
                          component = 1:ncomponents, itemp = 1:Ntemp,
                          imu = 1:Nmu
                total_L0[spin, component, itemp, imu] +=
                    reduced_L0[orbital, spin, component, itemp, imu]
            end
            TDF_decomp = BoltzDecompMoments(
                reduced_L0, reduced_L1, total_L0)
        end

        if boltz_setup.Write_TDF
            if boltz_setup.plane_type
                Write_TDF_3elements(
                    boltz_setup.filename, boltz_setup.mat_type, SpinPol,
                    TDF_Energy, boltz_setup.tau, TDF)
            else
                Write_TDF_6elements(
                    boltz_setup.filename, boltz_setup.mat_type, SpinPol,
                    TDF_Energy, boltz_setup.tau, TDF)
            end

            if keep_decomp_tdf
                println("Write $(boltz_setup.filename).TDFdecomp.jld2")
                jldopen("$(boltz_setup.filename).TDFdecomp.jld2", "w") do file
                    file["Dates"] = now()
                    file["filepath"] = boltz_setup.filepath
                    file["SpinPol"] = SpinPol
                    file["spinsize"] = spinsize
                    file["Nwann"] = Nwann
                    file["kmesh"] = boltz_setup.kmesh
                    file["decomp"] = boltz_setup.decomp
                    file["plane_type"] = boltz_setup.plane_type
                    file["tau"] = boltz_setup.tau
                    file["TDF_Erange"] = boltz_setup.TDF_Erange
                    file["TDF_dE"] = boltz_setup.TDF_dE
                    file["TDF_Energy"] = TDF_Energy
                    file["TDFdecomp"] = TDF_decomp
                end
            end
        end
    end

    MPI.Barrier(comm)
    return TDF_Energy, TDF, TDF_decomp
end
