@timeit timer "Set_OLPexp" function Set_OLPexp(material::LCPAO_model, tot_bvector, bvector)

    Nspin = material.Nspin
    TCpyCell = material.TCpyCell
    Latvecs = material.Latvecs
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    FNAN = material.FNAN
    natn = material.natn
    Gxyz = material.Gxyz
    Atoms_pao = material.Atoms_pao
    Total_NumOrbs = material.Total_NumOrbs    
    Ngrid = material.Ngrid
    Atoms_Cut1 =  material.Atoms_Cut1
    Grid_Origin = material.Grid_Origin


    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    pao = Vector{PAO}(undef, Nspecies)
    for spe = 1:Nspecies
        pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end
        
    ucell = UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs)
    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)

    
    OLPexp = Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef, tot_bvector)
	for ib = 1:tot_bvector
		OLPexp[ib] = Vector{Vector{Matrix{ComplexF64}}}(undef, Natom)
		for atom = 1:Natom
			OLPexp[ib][atom] = Vector{Matrix{ComplexF64}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				OLPexp[ib][atom][Rn] = zeros(ComplexF64, Total_NumOrbs[atom], Total_NumOrbs[natn[atom][Rn]])
			end
		end
	end
    Set_OLPexp_opt!(tot_bvector, bvector, OLPexp, Orbs_Grid, ucell)


    return OLPexp
end


function Set_OLPexp_opt!(tot_bvector, bvector, OLPexp, Orbs_Grid, ucell::UCell)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    MPI_atom = system_grid.MPI_atom
	Total_NumOrbs = system_grid.Total_NumOrbs
    Ngrid1, Ngrid2, Ngrid3 = ucell.Ngrid
    Grid_Origin = system_grid.Grid_Origin
    atv = system_grid.atv
    Gxyz = system_grid.Gxyz
    GridVol = system_grid.GridVol
    
    MPI_size = system_grid.MPI_size
    Total_Hsize = system_grid.Total_Hsize
    MPI_natn = system_grid.MPI_natn
    MPHks = system_grid.MPHks
    Hks_Num = MPHks[myrank+1]
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3
    
    
    hst = 0
    block_starts = zeros(Int, MPI_size)
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        block_starts[loop] = hst
        hst += NO0*NO1
    end

    relx = Vector{Vector{Float64}}(undef, MPI_size)
    rely = Vector{Vector{Float64}}(undef, MPI_size)
    relz = Vector{Vector{Float64}}(undef, MPI_size)
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        count = MPI_NumOLG[loop]
        relx[loop] = zeros(Float64, count)
        rely[loop] = zeros(Float64, count)
        relz[loop] = zeros(Float64, count)
        glist = GridListAtom[atom]
        clist = CellListAtom[atom]
        list1 = MPI_GListTAtoms1[loop]
        @inbounds for Nog = 1:count
            Nc = list1[Nog] + 1
            GN = glist[Nc]
            cell = clist[Nc] + 1
            n1 = div(GN, Ngrid2*Ngrid3)
            n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
            n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3
            Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
            Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
            Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]
            relx[loop][Nog] = Cx + atv[cell][1] - Gxyz[atom][1]
            rely[loop][Nog] = Cy + atv[cell][2] - Gxyz[atom][2]
            relz[loop][Nog] = Cz + atv[cell][3] - Gxyz[atom][3]
        end
    end

    tmp = zeros(ComplexF64, tot_bvector, Total_Hsize)
    bx = [bvector[ib][1] for ib = 1:tot_bvector]
    by = [bvector[ib][2] for ib = 1:tot_bvector]
    bz = [bvector[ib][3] for ib = 1:tot_bvector]
    max_orbs = maximum(Total_NumOrbs)
    max_olg = isempty(MPI_NumOLG) ? 0 : maximum(MPI_NumOLG)
    grid_block_size = min(max_olg, 4096)
    phi1 = zeros(ComplexF64, max_orbs, grid_block_size)
    phi2 = zeros(Float64, max_orbs, grid_block_size)
    weighted_phi2 = zeros(ComplexF64, max_orbs, grid_block_size)
    phase = zeros(ComplexF64, tot_bvector, grid_block_size)
    block_work = zeros(ComplexF64, max_orbs, max_orbs, tot_bvector)

    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        list1 = MPI_GListTAtoms1[loop]
        list2 = MPI_GListTAtoms2[loop]
        block_start = block_starts[loop]
        count = MPI_NumOLG[loop]
        fill!(block_work, 0.0)
        for first_grid = 1:4096:count
            block_count = min(4096, count-first_grid+1)
            @inbounds for local_grid = 1:block_count
                Nog = first_grid + local_grid - 1
                Nc = list1[Nog] + 1
                Nh = list2[Nog] + 1
                x = relx[loop][Nog]
                y = rely[loop][Nog]
                z = relz[loop][Nog]
                for ib = 1:tot_bvector
                    phase[ib,local_grid] = cis(-(bx[ib]*x + by[ib]*y + bz[ib]*z))*GridVol
                end
                orb1 = Orbs_Grid[atom][Nc]
                orb2 = Orbs_Grid[jatom][Nh]
                for ist = 1:NO0
                    phi1[ist,local_grid] = orb1[ist]
                end
                for jst = 1:NO1
                    phi2[jst,local_grid] = orb2[jst]
                end
            end

            phi1_block = view(phi1, 1:NO0, 1:block_count)
            phi2_block = view(phi2, 1:NO1, 1:block_count)
            weighted_block = view(weighted_phi2, 1:NO1, 1:block_count)
            @inbounds for ib = 1:tot_bvector
                for local_grid = 1:block_count
                    p = phase[ib,local_grid]
                    for jst = 1:NO1
                        weighted_block[jst,local_grid] = p*phi2_block[jst,local_grid]
                    end
                end
                result_block = view(block_work, 1:NO0, 1:NO1, ib)
                mul!(result_block, phi1_block, transpose(weighted_block), 1.0, 1.0)
            end
        end

        @inbounds for ib = 1:tot_bvector
            result_block = view(block_work, 1:NO0, 1:NO1, ib)
            for ist = 1:NO0, jst = 1:NO1
                idx = Hks_Num + block_start + (ist-1)*NO1 + jst
                tmp[ib,idx] = result_block[ist,jst]
            end
        end
    end

    MPI.Allreduce!(tmp, MPI.SUM, comm)
    for ib = 1:tot_bvector
        _Set_OLPexp!(view(tmp,ib,:), OLPexp[ib], Natom, FNAN, natn, Total_NumOrbs)
    end
end


function _Set_OLPexp!(OLPexp_tmp, OLPexp, Natom, FNAN, natn, Total_NumOrbs)
    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        hst += 1
        OLPexp[atom][Rn][ist,jst] = OLPexp_tmp[hst]
    end
end
