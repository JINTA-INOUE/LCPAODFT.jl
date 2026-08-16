@timeit timer "Set_OLPexp" function Set_OLPexp(material::LCPAO_model, tot_bvector, bvector)

    Nspin = material.Nspin
    Latvecs = material.Latvecs
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    TCpyCell = material.TCpyCell
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

    
    
    OLPexp = Vector{Vector{Vector{Vector{Vector{ComplexF64}}}}}(undef, tot_bvector)
	for ib = 1:tot_bvector
		OLPexp[ib] = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Natom)
		for atom = 1:Natom
			OLPexp[ib][atom] = Vector{Vector{Vector{ComplexF64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				OLPexp[ib][atom][Rn] = Vector{Vector{ComplexF64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					OLPexp[ib][atom][Rn][ist] = zeros(ComplexF64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end
    Set_OLPexp!(tot_bvector, bvector, OLPexp, Orbs_Grid, ucell)

    
    return OLPexp
end


function Set_OLPexp!(tot_bvector, bvector, OLPexp, Orbs_Grid, ucell::UCell)

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
    
    system_grid = ucell.system_grid
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


    MPI_hst = Vector{Vector{Vector{Int32}}}(undef, MPI_size)
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        MPI_hst[loop] = Vector{Vector{Int32}}(undef, NO0)
        for ist = 1:NO0
            MPI_hst[loop][ist] = zeros(Int32, NO1)
        end
    end

    hst = 0
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        for ist = 1:NO0, jst = 1:NO1
            hst += 1
            MPI_hst[loop][ist][jst] = hst
        end
    end


    OLPexp_tmp = zeros(ComplexF64, Total_Hsize)

    
    for ib = 1:tot_bvector
        fill!(OLPexp_tmp, 0.0)
        for loop = 1:MPI_size

            atom = MPI_atom[loop]
            jatom = MPI_natn[loop]
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]

            _GridListAtom = GridListAtom[atom]
            _CellListAtom = CellListAtom[atom]
            _MPI_GListTAtoms1 = MPI_GListTAtoms1[loop]
            _MPI_GListTAtoms2 = MPI_GListTAtoms2[loop]

            bx, by, bz = bvector[ib]

            for Nog = 1:MPI_NumOLG[loop]

                Nc = _MPI_GListTAtoms1[Nog]+1
                Nh = _MPI_GListTAtoms2[Nog]+1
                GN = _GridListAtom[Nc]
                cell = _CellListAtom[Nc]+1
                    
                n1 = div(GN, Ngrid2*Ngrid3)
                n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
                n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

                Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
                Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
                Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

                x = Cx + atv[cell][1] - Gxyz[atom][1]
                y = Cy + atv[cell][2] - Gxyz[atom][2]
                z = Cz + atv[cell][3] - Gxyz[atom][3]

                ex = exp(-im*(bx*x + by*y + bz*z))
                fac = ex*GridVol

                Orbs_Grid1 = Orbs_Grid[atom][Nc]
                Orbs_Grid2 = Orbs_Grid[jatom][Nh]
                    
                for ist = 1:NO0
                    phi1 = Orbs_Grid1[ist]
                    @inbounds for jst = 1:NO1
                        hst = MPI_hst[loop][ist][jst]
                        OLPexp_tmp[Hks_Num+hst] += fac*Orbs_Grid2[jst]*phi1
                    end
                end
            end
        end


        MPI.Allreduce!(OLPexp_tmp, MPI.SUM, comm)
        _Set_OLPexp!(OLPexp_tmp, OLPexp[ib], Natom, FNAN, natn, Total_NumOrbs)
    end
end



function _Set_OLPexp!(OLPexp_tmp, OLPexp, Natom, FNAN, natn, Total_NumOrbs)
    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        hst += 1
        OLPexp[atom][Rn][ist][jst] = OLPexp_tmp[hst]
    end
end
