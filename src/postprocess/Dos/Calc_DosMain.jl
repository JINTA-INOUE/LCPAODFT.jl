function Calc_DosMain(filename::String,
                      material::Union{LCPAO_model,CWF_model}, Enk,
                      halo_plan::DosHaloPlan, neg, kmesh, Dos_Erange;
                      de_Dos=0.01)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    if myrank == 0
        println("<DosMain>  Generate Density of States using the tetrahedron method")
    end

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    Nkpt = prod(kmesh)

    DosEmin, DosEmax = Dos_Erange
    de = de_Dos / eV2Hartree
    Dos_N = floor(Int, (DosEmax - DosEmin) / de)
    Dos_N >= 2 || error("the DOS energy grid must contain at least two points")

    DosE = zeros(Float64, Dos_N)
    @inbounds for ie = 1:Dos_N
        DosE[ie] = DosEmin +
                   (DosEmax - DosEmin) * (ie - 1) / (Dos_N - 1)
    end

    Dos = myrank == 0 ? zeros(Float64, Dos_N, spinsize) : nothing
    local_Dos = zeros(Float64, Dos_N)
    cell_e = zeros(Float64, 8)
    tetra_e = zeros(Float64, 4)
    tetra_id = (1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
               (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8)

    for spin = 1:spinsize
        fill!(local_Dos, 0.0)
        for ieg = 1:neg, cell_position = 1:halo_plan.nlocal
            @inbounds for vertex = 1:8
                k_position = halo_plan.vertex_positions[vertex, cell_position]
                cell_e[vertex] = _dos_halo_get(Enk, ieg, spin, k_position)
            end

            @inbounds for vertices in tetra_id
                for ic = 1:4
                    tetra_e[ic] = cell_e[vertices[ic]]
                end
                OrderE0!(tetra_e, 4)

                first_energy = trunc(Int,
                    (tetra_e[1] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) - 1)
                last_energy = trunc(Int,
                    (tetra_e[4] - DosEmin) / (DosEmax - DosEmin) *
                    (Dos_N - 1) + 1)

                # Preserve the serial DOS grid convention (one-based indices).
                first_energy = max(first_energy, 1)
                last_energy = min(last_energy, Dos_N - 1)
                if first_energy <= last_energy
                    for ie = first_energy:last_energy
                        local_Dos[ie] += ATM_Dos(tetra_e, DosE[ie])
                    end
                end
            end
        end

        if myrank == 0
            MPI.Reduce!(local_Dos, @view(Dos[:, spin]), MPI.SUM, comm;
                        root=0)
        else
            MPI.Reduce!(local_Dos, nothing, MPI.SUM, comm; root=0)
        end
    end

    if myrank == 0
        h = (DosEmax - DosEmin) / (Dos_N - 1) * eV2Hartree
        ssum = zeros(Float64, Dos_N, spinsize)
        factor = 1 / Nkpt / 6 / eV2Hartree

        @inbounds for spin = 1:spinsize
            for ie = 1:Dos_N
                Dos[ie, spin] *= factor
            end

            # Retain the original Simpson summation order for reproducibility.
            for q = 1:Dos_N
                s1 = 0.0
                s2 = 0.0
                for ie = 2:2:q - 1
                    s1 += Dos[ie, spin]
                end
                for ie = 3:2:q - 1
                    s2 += Dos[ie, spin]
                end
                ssum[q, spin] =
                    (Dos[1, spin] + 4 * s1 + 2 * s2 + Dos[q, spin]) * h / 3
            end
        end

        Write_Dos_Tetrahedron(filename, SpinPol, Dos_N, DosE, Dos, ssum)
        Write_Dos_gnuplot(filename, SpinPol, Dos_Erange)
    end

    MPI.Barrier(comm)
    return nothing
end
