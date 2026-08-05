function Calc_PDosMain(filename::String, material::LCPAO_model, Enk, EVec,
                       halo_plan::DosHaloPlan, neg, kmesh, Dos_Erange;
                       de_Dos=0.01)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        println("<PDosMain>  Generate Projected Density of States using the tetrahedron method")
    end

    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Atoms_pao = material.Atoms_pao
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol == "off", 1, 2)

    Nkpt = prod(kmesh)
    DosEmin, DosEmax = Dos_Erange
    de = de_Dos / eV2Hartree
    Dos_N = floor(Int, (DosEmax - DosEmin) / de)
    Dos_N >= 2 || error("the PDOS energy grid must contain at least two points")

    DosE = zeros(Float64, Dos_N)
    @inbounds for ie = 1:Dos_N
        DosE[ie] = DosEmin +
                   (DosEmax - DosEmin) * (ie - 1) / (Dos_N - 1)
    end

    Spe_orb = Get_Atoms_data(Atoms_pao)[3]
    Spe_MaxL_Basis = zeros(Int, Nspecies)
    Spe_Num_Basis = [zeros(Int, 4) for _ = 1:Nspecies]
    @inbounds for spe = 1:Nspecies
        Spe_MaxL_Basis[spe], Num_Basis = get_ialpha_index(Spe_orb[spe])
        for angular = 1:Spe_MaxL_Basis[spe] + 1
            Spe_Num_Basis[spe][angular] = Num_Basis[angular]
        end
    end

    Spe_Num_Relation = [zeros(Int, maximum(Total_NumOrbs))
                        for _ = 1:Nspecies]
    @inbounds for spe = 1:Nspecies
        id = 0
        for L = 0:3
            if Spe_Num_Basis[spe][L + 1] > 0
                for _ = 0:Spe_Num_Basis[spe][L + 1] - 1, M = 0:2 * L
                    id += 1
                    Spe_Num_Relation[spe][id] = 100 * L + M + 1
                end
            end
        end
    end

    # Only rank zero holds the complete output. Other ranks reuse one energy
    # vector while reducing each projected orbital to the root.
    Dos_dense = myrank == 0 ? zeros(Float64, Dos_N, fsize, spinsize) : nothing
    local_Dos = zeros(Float64, Dos_N)
    cell_e = zeros(Float64, 8)
    cell_a = zeros(Float64, 8)
    tetra_e = zeros(Float64, 4)
    tetra_a = zeros(Float64, 4)
    tetra_id = (1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
               (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8)

    for spin = 1:spinsize, orbital = 1:fsize
        fill!(local_Dos, 0.0)
        @inbounds for ieg = 1:neg, cell_position = 1:halo_plan.nlocal
            for vertex = 1:8
                k_position = halo_plan.vertex_positions[vertex, cell_position]
                cell_e[vertex] = _dos_halo_get(Enk, ieg, spin, k_position)
                cell_a[vertex] = _dos_halo_get(
                    EVec, ieg, orbital, spin, k_position)
            end

            for vertices in tetra_id
                for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                    tetra_a[ic] = cell_a[vertices[ic]]
                end
                OrderE!(tetra_e, tetra_a, 4)

                first_energy = trunc(Int,
                    (tetra_e[1] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) - 1)
                last_energy = trunc(Int,
                    (tetra_e[4] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) + 1)
                first_energy = max(first_energy, 0)
                last_energy = min(last_energy, Dos_N - 1)

                if first_energy <= last_energy
                    for ie = first_energy:last_energy
                        local_Dos[ie + 1] +=
                            ATM_Spectrum(tetra_e, tetra_a, DosE[ie + 1])
                    end
                end
            end
        end

        if myrank == 0
            MPI.Reduce!(local_Dos, @view(Dos_dense[:, orbital, spin]),
                        MPI.SUM, comm; root=0)
        else
            MPI.Reduce!(local_Dos, nothing, MPI.SUM, comm; root=0)
        end
    end

    if myrank == 0
        Dos_dense .*= 1 / Nkpt / 6 / eV2Hartree

        # Keep the established nested output interface local to rank zero.
        Dos = Vector{Vector{Vector{Vector{Float64}}}}(undef, spinsize)
        for spin = 1:spinsize
            Dos[spin] = Vector{Vector{Vector{Float64}}}(undef, Natom)
            for atom = 1:Natom
                Dos[spin][atom] = Vector{Vector{Float64}}(
                    undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    orbital = MP[atom] + ist
                    Dos[spin][atom][ist] =
                        collect(@view Dos_dense[:, orbital, spin])
                end
            end
        end

        Calc_PDos_Atom_proj(filename, material, Dos_Erange, DosE, Dos)
        Calc_PDos_Orbital_proj(filename, material, Spe_Num_Relation,
                              Spe_Num_Basis, Dos_Erange, DosE, Dos)
        Write_PDos_gnuplot(filename, Natom, Dos_Erange)
    end

    MPI.Barrier(comm)
    return nothing
end


function Calc_PDosMain(filename::String, material::CWF_model, Enk, EVec,
                       halo_plan::DosHaloPlan, neg, kmesh, Dos_Erange;
                       de_Dos=0.01)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        println("<PDosMain>  Generate Projected Density of States using the tetrahedron method")
    end

    gsize = Int(material.gsize)
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Nkpt = prod(kmesh)

    DosEmin, DosEmax = Dos_Erange
    de = de_Dos / eV2Hartree
    Dos_N = floor(Int, (DosEmax - DosEmin) / de)
    Dos_N >= 2 || error("the PDOS energy grid must contain at least two points")
    DosE = zeros(Float64, Dos_N)
    @inbounds for ie = 1:Dos_N
        DosE[ie] = DosEmin +
                   (DosEmax - DosEmin) * (ie - 1) / (Dos_N - 1)
    end

    Dos_dense = myrank == 0 ? zeros(Float64, Dos_N, gsize, spinsize) : nothing
    local_Dos = zeros(Float64, Dos_N)
    cell_e = zeros(Float64, 8)
    cell_a = zeros(Float64, 8)
    tetra_e = zeros(Float64, 4)
    tetra_a = zeros(Float64, 4)
    tetra_id = (1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
               (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8)

    for spin = 1:spinsize, orbital = 1:gsize
        fill!(local_Dos, 0.0)
        @inbounds for ieg = 1:neg, cell_position = 1:halo_plan.nlocal
            for vertex = 1:8
                k_position = halo_plan.vertex_positions[vertex, cell_position]
                cell_e[vertex] = _dos_halo_get(Enk, ieg, spin, k_position)
                cell_a[vertex] = _dos_halo_get(
                    EVec, ieg, orbital, spin, k_position)
            end

            for vertices in tetra_id
                for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                    tetra_a[ic] = cell_a[vertices[ic]]
                end
                OrderE!(tetra_e, tetra_a, 4)

                first_energy = trunc(Int,
                    (tetra_e[1] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) - 1)
                last_energy = trunc(Int,
                    (tetra_e[4] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) + 1)
                first_energy = max(first_energy, 0)
                last_energy = min(last_energy, Dos_N - 1)

                if first_energy <= last_energy
                    for ie = first_energy:last_energy
                        local_Dos[ie + 1] +=
                            ATM_Spectrum(tetra_e, tetra_a, DosE[ie + 1])
                    end
                end
            end
        end

        if myrank == 0
            MPI.Reduce!(local_Dos, @view(Dos_dense[:, orbital, spin]),
                        MPI.SUM, comm; root=0)
        else
            MPI.Reduce!(local_Dos, nothing, MPI.SUM, comm; root=0)
        end
    end

    if myrank == 0
        Dos_dense .*= 1 / Nkpt / 6 / eV2Hartree
        Dos = Vector{Vector{Vector{Float64}}}(undef, spinsize)
        for spin = 1:spinsize
            Dos[spin] = Vector{Vector{Float64}}(undef, gsize)
            for orbital = 1:gsize
                Dos[spin][orbital] =
                    collect(@view Dos_dense[:, orbital, spin])
            end
        end

        Calc_PDos_Orbital_proj(filename, material, Dos_Erange, DosE, Dos)
    end

    MPI.Barrier(comm)
    return nothing
end
