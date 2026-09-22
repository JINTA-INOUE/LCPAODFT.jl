@timeit timer "Set_MLWF_Grid" function Set_MLWF_Grid(mlwf_setup::MLWF_Setup, MLWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    TCpyCell = material.TCpyCell
    Nspin = material.Nspin
    SpinPol = material.SpinPol
    Latvecs = material.Latvecs
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Gxyz = material.Gxyz
    Atoms_pao = material.Atoms_pao
    Total_NumOrbs = material.Total_NumOrbs
    Ngrid = material.Ngrid
    Atoms_Cut1 = material.Atoms_Cut1
    Grid_Origin = material.Grid_Origin
    MLWF_Plot_Cube = mlwf_setup.MLWF_Plot_Cube
    MLWF_Plot_SuperCells = mlwf_setup.MLWF_Plot_SuperCells
    filename = mlwf_setup.filename

    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    pao = Vector{PAO}(undef, Nspecies)
    for spe = 1:Nspecies
        pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end
        
    ucell = UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs)
    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)


    if SpinPol ∈ ("off", "on")
        Set_CWF_Grid_Col(filename, material, ucell, MLWF_Plot_Cube, MLWF_Plot_SuperCells, MLWF_ExpnCoef, Orbs_Grid)
    elseif SpinPol == "nc"
        Set_CWF_Grid_NonCol(filename, material, ucell, MLWF_Plot_Cube, MLWF_Plot_SuperCells, MLWF_ExpnCoef, Orbs_Grid)
    end
    MPI.Barrier(comm)
end