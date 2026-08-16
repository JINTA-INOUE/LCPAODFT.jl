@timeit timer "Set_OLPpos" function Set_OLPpos(Orbs_Grid, ucell::UCell)

    Natom = ucell.system_grid.Natom
	FNAN = ucell.system_grid.FNAN
	natn = ucell.system_grid.natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs

    OLPpos = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
    for xyz = 1:3
        OLPpos[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
        for atom = 1:Natom
            OLPpos[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                OLPpos[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    OLPpos[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end
    Set_OLPpos!(OLPpos, Orbs_Grid, ucell)


    return OLPpos
end


function Set_OLPpos!(OLPpos, Orbs_Grid, ucell::UCell)

    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
	natn = system_grid.natn
	Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_size = system_grid.MPI_size
    Grid_Origin = system_grid.Grid_Origin
    atv = system_grid.atv
    Gxyz = system_grid.Gxyz
    GridVol = system_grid.GridVol
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Ngrid1, Ngrid2, Ngrid3 = ucell.Ngrid
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
    gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
    gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        _GridListAtom = GridListAtom[atom]
        _CellListAtom = CellListAtom[atom]
        _MPI_GListTAtoms1 = MPI_GListTAtoms1[loop]
        _MPI_GListTAtoms2 = MPI_GListTAtoms2[loop]

        OLPpos1 = OLPpos[1][atom][Rn]
        OLPpos2 = OLPpos[2][atom][Rn]
        OLPpos3 = OLPpos[3][atom][Rn]

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

            Orbs_Grid1 = Orbs_Grid[atom][Nc]
            Orbs_Grid2 = Orbs_Grid[jatom][Nh]

            tmpx = x*GridVol
            tmpy = y*GridVol
            tmpz = z*GridVol

            for ist = 1:NO0
                phi1 = Orbs_Grid1[ist]
                @inbounds for jst = 1:NO1
                    OLPpos1[ist][jst] += tmpx*Orbs_Grid2[jst]*phi1
                    OLPpos2[ist][jst] += tmpy*Orbs_Grid2[jst]*phi1
                    OLPpos3[ist][jst] += tmpz*Orbs_Grid2[jst]*phi1
                end
            end
        end
    end
end
