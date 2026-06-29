@timeit timer "Set_CWF_Grid" function Set_CWF_Grid(CWF_ExpnCoef, Orbs_Grid, ucell::UCell, cwf_setup::Union{CWF_Setup,CWF_Setup_MO})

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    SpinPol = material.SpinPol


    if SpinPol ∈ ("off", "on")
        Set_CWF_Grid_Col(cwf_setup, ucell, CWF_ExpnCoef, Orbs_Grid)
    elseif SpinPol == "nc"
        Set_CWF_Grid_NonCol(cwf_setup, ucell, CWF_ExpnCoef, Orbs_Grid) 
    end
    MPI.Barrier(comm)
end


function Set_CWF_Grid_Col(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, ucell::UCell, CWF_ExpnCoef, Orbs_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    Natom = material.Natom
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Atoms_symbol = material.Atoms_symbol
    Gxyz_scf = material.Gxyz
    
    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    Ngrid = system_grid.Ngrid
    Grid_Origin = system_grid.Grid_Origin
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    spinsize = cwf_setup.spinsize
    CWF_Plot_Cube = cwf_setup.CWF_Plot_Cube
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_TNumGrid = prod(Ngrid)*Plot_NCell
    filename = cwf_setup.filename


    CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid = CWF_UCell(cwf_setup, ucell)
    MPI.Barrier(comm)

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    N_Polt_Cube = length(CWF_Plot_Cube)
    Nloop = spinsize*N_Polt_Cube
    OneD2spin = zeros(Int32, Nloop)
    OneD2proj = zeros(Int32, Nloop)

    counts = 1
    for spin = 1:spinsize, proj in CWF_Plot_Cube
        OneD2spin[counts] = spin
        OneD2proj[counts] = proj
        counts += 1
    end
    
    myrange = split_evenly(1:Nloop, nprocs)
    MPI_CWF_size = length(myrange[myrank+1])
    MPI_spin = OneD2spin[myrange[myrank+1]]
    MPI_proj = OneD2proj[myrange[myrank+1]]


    ExpnCoef = zeros(Float64, fsize)
    Wannier_Orbs_Grid = zeros(Float64, CWF_TNumGrid)

    for loop = 1:MPI_CWF_size
        spin = MPI_spin[loop]
        proj = MPI_proj[loop]
        println("  Write myrank = $myrank   $(filename)_CWF$(spin)_$(proj).cube")
        MPI.Barrier(comm)
        
        fill!(Wannier_Orbs_Grid, 0.0)
        for cell = 1:Plot_NCell
            @. ExpnCoef = CWF_ExpnCoef[spin][proj][cell]
            for atom = 1:Natom
                cwf_proj = MP[atom]
                NO0 = Total_NumOrbs[atom]
                _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid[atom], ExpnCoef, Wannier_Orbs_Grid)
            end
        end

        data = open("$(filename)_CWF$(spin)_$(proj).cube", "w")
        Write_Wannier_CubeInfo(data, Atoms_symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Gxyz_scf, Latvecs, gLatvecs, Ngrid)
        Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, Wannier_Orbs_Grid)
        close(data)
    end
end


function Set_CWF_Grid_NonCol(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, ucell::UCell, CWF_ExpnCoef, Orbs_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    Natom = material.Natom
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    Nfsize = 2*sum(Total_NumOrbs)
    Atoms_symbol = material.Atoms_symbol
    Gxyz_scf = material.Gxyz
    
    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    Ngrid = system_grid.Ngrid
    Grid_Origin = system_grid.Grid_Origin
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    CWF_Plot_Cube = cwf_setup.CWF_Plot_Cube
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_TNumGrid = prod(Ngrid)*Plot_NCell
    filename = cwf_setup.filename


    CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid = CWF_UCell(cwf_setup, ucell)

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    N_Polt_Cube = length(CWF_Plot_Cube)
    Nloop = N_Polt_Cube
    OneD2proj = zeros(Int32, Nloop)

    counts = 1
    for proj in CWF_Plot_Cube
        OneD2proj[counts] = proj
        counts += 1
    end
    
    myrange = split_evenly(1:Nloop, nprocs)
    MPI_CWF_size = length(myrange[myrank+1])
    MPI_proj = OneD2proj[myrange[myrank+1]]



    ExpnCoef = zeros(ComplexF64, Nfsize)
    Wannier_Orbs_Grid = zeros(ComplexF64, CWF_TNumGrid)

    for loop = 1:MPI_CWF_size
        proj = MPI_proj[loop]
        println("  Write myrank = $myrank   $(filename)_CWF$(proj).nccube")
        data = open("$(filename)_CWF$(proj).nccube", "w")
        Write_Wannier_CubeInfo(data, Atoms_symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Gxyz_scf, Latvecs, gLatvecs, Ngrid)
        for spin = 1:2
            fill!(Wannier_Orbs_Grid, 0.0)
            spin_site = ifelse(spin==1, 0, fsize)
            for cell = 1:Plot_NCell
                @. ExpnCoef = CWF_ExpnCoef[proj][cell]
                for atom = 1:Natom
                    cwf_proj = MP[atom]
                    NO0 = Total_NumOrbs[atom]
                    _Calc_CWF_Grid8!(spin_site+cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid[atom], ExpnCoef, Wannier_Orbs_Grid)
                end
            end

            Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, real(Wannier_Orbs_Grid))
            Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, imag(Wannier_Orbs_Grid))
        end
        close(data)
    end
end


function CWF_UCell(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, ucell::UCell)

    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    CWF_Plot_SuperCells1, CWF_Plot_SuperCells2, CWF_Plot_SuperCells3 = CWF_Plot_SuperCells

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    atv_ijk = system_grid.atv_ijk
    Ngrid = system_grid.Ngrid
    Grid_Origin = system_grid.Grid_Origin
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        Plot_cell_ijk[cell] = zeros(Int32, 3)
    end
    cell = 0
    for l1 = -CWF_Plot_SuperCells1:CWF_Plot_SuperCells1, l2 = -CWF_Plot_SuperCells2:CWF_Plot_SuperCells2, l3 = -CWF_Plot_SuperCells3:CWF_Plot_SuperCells3
        cell += 1
        Plot_cell_ijk[cell] = [l1, l2, l3]
    end


    CWF_GridN_Atom = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        CWF_GridN_Atom[cell] = zeros(Int32, Natom)
    end

    CWF_GridOrbs_Grid = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    CWF_GridListAtom = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        CWF_GridOrbs_Grid[cell] = Vector{Vector{Int32}}(undef, Natom)
        CWF_GridListAtom[cell] = Vector{Vector{Int32}}(undef, Natom)
    end


    for cell = 1:Plot_NCell

        l1, l2, l3 = Plot_cell_ijk[cell]
        for atom = 1:Natom

            gridn = 0
            for Nc = 1:GridN_Atom[atom]

                GRc = CellListAtom[atom][Nc]+1

                m1 = l1 + atv_ijk[GRc][1]
                m2 = l2 + atv_ijk[GRc][2]
                m3 = l3 + atv_ijk[GRc][3]
                                                
                if abs(m1)<=CWF_Plot_SuperCells1 && abs(m2)<=CWF_Plot_SuperCells2 && abs(m3)<=CWF_Plot_SuperCells3
                    gridn += 1
                end
            end
            
            CWF_GridN_Atom[cell][atom] = gridn
            CWF_GridOrbs_Grid[cell][atom] = zeros(Int32, gridn)
            CWF_GridListAtom[cell][atom] = zeros(Int32, gridn)

            gridn = 0
            for Nc = 1:GridN_Atom[atom]

                GNc = GridListAtom[atom][Nc]
                GRc = CellListAtom[atom][Nc]+1

                n1 = div(GNc, Ngrid2*Ngrid3)
                n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
                n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3

                m1 = l1 + atv_ijk[GRc][1]
                m2 = l2 + atv_ijk[GRc][2]
                m3 = l3 + atv_ijk[GRc][3]
                        
                if abs(m1)<=CWF_Plot_SuperCells1 && abs(m2)<=CWF_Plot_SuperCells2 && abs(m3)<=CWF_Plot_SuperCells3

                    p1 = n1 + (m1 + CWF_Plot_SuperCells1)*Ngrid1
                    p2 = n2 + (m2 + CWF_Plot_SuperCells2)*Ngrid2
                    p3 = n3 + (m3 + CWF_Plot_SuperCells3)*Ngrid3

                    GN1 = p1*Ngrid2*(2*CWF_Plot_SuperCells2+1)*Ngrid3*(2*CWF_Plot_SuperCells3+1) + p2*Ngrid3*(2*CWF_Plot_SuperCells3+1) + p3 + 1

                    gridn += 1
        
                    CWF_GridOrbs_Grid[cell][atom][gridn] = Nc
                    CWF_GridListAtom[cell][atom][gridn] = GN1
                end
            end
        end
    end


    return CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid
end


function _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom, CWF_GridOrbs_Grid, CWF_GridListAtom, Orbs_Grid, ExpnCoef, Wannier_Orbs_Grid)

    for Noc = 1:8:CWF_GridN_Atom-7

        Nc0 = CWF_GridOrbs_Grid[Noc]
        Nc1 = CWF_GridOrbs_Grid[Noc+1]
        Nc2 = CWF_GridOrbs_Grid[Noc+2]
        Nc3 = CWF_GridOrbs_Grid[Noc+3]
        Nc4 = CWF_GridOrbs_Grid[Noc+4]
        Nc5 = CWF_GridOrbs_Grid[Noc+5]
        Nc6 = CWF_GridOrbs_Grid[Noc+6]
        Nc7 = CWF_GridOrbs_Grid[Noc+7]

        GN0 = CWF_GridListAtom[Noc]
        GN1 = CWF_GridListAtom[Noc+1]
        GN2 = CWF_GridListAtom[Noc+2]
        GN3 = CWF_GridListAtom[Noc+3]
        GN4 = CWF_GridListAtom[Noc+4]
        GN5 = CWF_GridListAtom[Noc+5]
        GN6 = CWF_GridListAtom[Noc+6]
        GN7 = CWF_GridListAtom[Noc+7]

        temp0 = 0.0
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        temp4 = 0.0
        temp5 = 0.0
        temp6 = 0.0
        temp7 = 0.0
        for ist = 1:NO0
            Coef = ExpnCoef[cwf_proj+ist]
            temp0 += Coef*Orbs_Grid[Nc0][ist]
            temp1 += Coef*Orbs_Grid[Nc1][ist]
            temp2 += Coef*Orbs_Grid[Nc2][ist]
            temp3 += Coef*Orbs_Grid[Nc3][ist]
            temp4 += Coef*Orbs_Grid[Nc4][ist]
            temp5 += Coef*Orbs_Grid[Nc5][ist]
            temp6 += Coef*Orbs_Grid[Nc6][ist]
            temp7 += Coef*Orbs_Grid[Nc7][ist]
        end
        Wannier_Orbs_Grid[GN0] += temp0
        Wannier_Orbs_Grid[GN1] += temp1
        Wannier_Orbs_Grid[GN2] += temp2
        Wannier_Orbs_Grid[GN3] += temp3
        Wannier_Orbs_Grid[GN4] += temp4
        Wannier_Orbs_Grid[GN5] += temp5
        Wannier_Orbs_Grid[GN6] += temp6
        Wannier_Orbs_Grid[GN7] += temp7
    end



    Nog1 = 8*div(CWF_GridN_Atom, 8)
	rem_CWF_GridN_Atom = rem(CWF_GridN_Atom, 8)
    for Nog = 1:rem_CWF_GridN_Atom
        Nc = CWF_GridOrbs_Grid[Nog+Nog1]
		GN = CWF_GridListAtom[Nog+Nog1]

        Sum = 0.0
        for ist = 1:NO0
            Sum += ExpnCoef[cwf_proj+ist]*Orbs_Grid[Nc][ist]
        end
        Wannier_Orbs_Grid[GN] += Sum
    end
end