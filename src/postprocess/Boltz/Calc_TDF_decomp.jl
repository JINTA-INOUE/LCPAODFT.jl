@timeit timer "Calc_TDF_decomp" function Calc_TDF_decomp(
    boltz_setup::Boltz_Setup, Enk::BoltzHaloData,
    EVec::BoltzHaloData, Vnk::BoltzHaloData,
    halo_plan::BoltzHaloPlan)

    material = boltz_setup.material
    SpinPol = material.SpinPol
    spinsize = SpinPol == "on" ? 2 : 1
    Nwann = Int(material.Ngsize)
    plane_type = boltz_setup.plane_type
    ncomponents = plane_type ? 3 : 6
    velocity_pairs = plane_type ?
        ((1, 1), (1, 2), (2, 2)) :
        ((1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3))

    TDF_Emin, TDF_Emax = boltz_setup.TDF_Erange
    TDF_EneNum = floor(Int, (TDF_Emax - TDF_Emin) / boltz_setup.TDF_dE)
    TDF_EneNum >= 2 || error("the TDF energy grid must contain at least two points")
    TDF_Energy = zeros(Float64, TDF_EneNum)
    for ie = 1:TDF_EneNum
        TDF_Energy[ie] = TDF_Emin +
            (TDF_Emax - TDF_Emin) * (ie - 1) / (TDF_EneNum - 1)
    end

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    local_TDF = zeros(Float64, TDF_EneNum, spinsize, ncomponents)
    TDF_decomp = myrank == 0 ?
        [similar(local_TDF) for _ = 1:Nwann] : nothing

    cell_e = zeros(Float64, 8)
    cell_a = zeros(Float64, 8, ncomponents)
    tetra_e = zeros(Float64, 4)
    tetra_a = zeros(Float64, 4)
    tetra_id = ((1, 2, 3, 6), (2, 3, 4, 6), (3, 4, 6, 8),
                (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8))

    # Reduce one orbital at a time so non-root ranks never hold the complete
    # decomposed output array.
    for orbital = 1:Nwann
        fill!(local_TDF, 0.0)
        for spin = 1:spinsize, band = 1:Nwann
            for cell_position = 1:halo_plan.nlocal
                @inbounds for vertex = 1:8
                    k_position =
                        halo_plan.vertex_positions[vertex, cell_position]
                    cell_e[vertex] =
                        _boltz_halo_get(Enk, band, spin, k_position)
                    weight = _boltz_halo_get(
                        EVec, orbital, band, spin, k_position)
                    for component = 1:ncomponents
                        xyz1, xyz2 = velocity_pairs[component]
                        velocity1 = _boltz_halo_get(
                            Vnk, band, xyz1, spin, k_position)
                        velocity2 = _boltz_halo_get(
                            Vnk, band, xyz2, spin, k_position)
                        cell_a[vertex, component] =
                            weight * velocity1 * velocity2
                    end
                end

                for vertices in tetra_id, component = 1:ncomponents
                    @inbounds for ic = 1:4
                        tetra_e[ic] = cell_e[vertices[ic]]
                        tetra_a[ic] = cell_a[vertices[ic], component]
                    end
                    OrderE!(tetra_e, tetra_a, 4)

                    # Keep the established decomposition indexing convention.
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
                        for ie = first_energy:last_energy
                            local_TDF[ie, spin, component] +=
                                ATM_Spectrum(tetra_e, tetra_a,
                                             TDF_Energy[ie + 1])
                        end
                    end
                end
            end
        end

        if myrank == 0
            MPI.Reduce!(local_TDF, TDF_decomp[orbital],
                        MPI.SUM, comm; root=0)
        else
            MPI.Reduce!(local_TDF, nothing, MPI.SUM, comm; root=0)
        end
    end

    if myrank == 0
        cell_volume_Ang = abs(det(material.Latvecs)) / Ang_to_bohr^3
        factor = boltz_setup.tau /
                 (prod(boltz_setup.kmesh) * 6 * cell_volume_Ang)
        for orbital = 1:Nwann
            TDF_decomp[orbital] .*= factor
        end

        if boltz_setup.Write_TDF
            println("Write $(boltz_setup.filename).TDFdecomp.jld2")
            jldopen("$(boltz_setup.filename).TDFdecomp.jld2", "w") do file
                file["Dates"] = now()
                file["filepath"] = boltz_setup.filepath
                file["SpinPol"] = SpinPol
                file["spinsize"] = spinsize
                file["Nwann"] = Nwann
                file["kmesh"] = boltz_setup.kmesh
                file["decomp"] = boltz_setup.decomp
                file["plane_type"] = plane_type
                file["tau"] = boltz_setup.tau
                file["TDF_Erange"] = boltz_setup.TDF_Erange
                file["TDF_dE"] = boltz_setup.TDF_dE
                file["TDF_Energy"] = TDF_Energy
                file["TDFdecomp"] = TDF_decomp
            end
        end
    end

    MPI.Barrier(comm)


    return TDF_Energy, TDF_decomp
end
