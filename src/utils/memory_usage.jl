function memory_usage(ucell::UCell)

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn


end


function memory_usage_UCell(ucell::UCell)
    system_grid =  ucell.system_grid
    Natom = system_grid.Natom
    GridN_Atom = ucell.GridN_Atom
    MPI_NumOLG = ucell.MPI_NumOLG
    size_GridCellListAtom = sum(GridN_Atom)
    size_GListTAtoms = sum(MPI_NumOLG)

    
end


function memory_usage_System_Grid()

end


function memory_usage_Mixing()

end


function memory_usage_electron()

end


function memory_usage_Matrix()

end


function memory_usage_Potentials()
    
end