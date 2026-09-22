@timeit timer "Set_Density_Grid" function Set_Density_Grid!(SpinPol::String, ucell::UCell, Orbs_Grid, DM, Density_Grid)
    if SpinPol == "off"
        Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM, Density_Grid)
    elseif SpinPol == "on"
        Set_Density_Grid_pol!(ucell, Orbs_Grid, DM, Density_Grid)
    elseif SpinPol == "nc"
        Set_Density_Grid_nc!(ucell, Orbs_Grid, DM, Density_Grid)
    end
end


"""A zero-copy view of one atom-pair block in the flattened density matrices."""
struct _DensityBlock{N,T,V<:AbstractVector{T}}
    data::NTuple{N,V}
    offset::Int
    nrows::Int
end


@inline _DensityBlock(data::NTuple{N,V}, offset::Integer, nrows::Integer) where {N,T,V<:AbstractVector{T}} =
    _DensityBlock{N,T,V}(data, Int(offset), Int(nrows))

@inline function Base.getindex(block::_DensityBlock, j::Int, i::Int)
    @inbounds return block.data[1][block.offset + (i - 1)*block.nrows + j]
end

@inline function Base.getindex(block::_DensityBlock, j::Int, i::Int, spin::Int)
    @inbounds return block.data[spin][block.offset + (i - 1)*block.nrows + j]
end


"""Map atom-local grid indices directly onto a rank-local global density grid."""
struct _MappedDensityGrid{N,T,V<:AbstractVector{T},I<:AbstractVector{<:Integer}}
    data::NTuple{N,V}
    grid_indices::I
end

@inline function Base.getindex(grid::_MappedDensityGrid, local_index::Int)
    @inbounds return grid.data[1][grid.grid_indices[local_index] + 1]
end


@inline function Base.setindex!(grid::_MappedDensityGrid, value, local_index::Int)
    @inbounds grid.data[1][grid.grid_indices[local_index] + 1] = value
    return value
end


@inline function Base.getindex(grid::_MappedDensityGrid, local_index::Int, spin::Int)
    @inbounds return grid.data[spin][grid.grid_indices[local_index] + 1]
end


@inline function Base.setindex!(grid::_MappedDensityGrid, value, local_index::Int, spin::Int)
    @inbounds grid.data[spin][grid.grid_indices[local_index] + 1] = value
    return value
end


function Set_Density_Grid_nonpol!(ucell::UCell, Orbs_Grid, DM, Density_Grid)
   
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPI_Hoffset = ucell.system_grid.MPI_Hoffset

    DM1 = DM[1]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch

    fill!(Density_Grid[1], 0.0)

    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        block_size = NO0*NO1
        dm_index = MPI_Hoffset[loop]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            dm_index += 1
            density_matrix_scratch[jst, ist, 1] = DM1[dm_index]
        end
        @inbounds for xyz = 1:GridN_Atom[atom]
            density_scratch[xyz, 1] = 0.0
        end
        orbitals1 = atom_matrix(Orbs_Grid, atom)
        orbitals2 = atom_matrix(Orbs_Grid, jatom)
        _Calc_Density_Blocked!(density_scratch, 1, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbitals1, orbitals2, density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            Density_Grid[1][GridListAtom[atom][xyz] + 1] += density_scratch[xyz, 1]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
end


function Set_Density_Grid_pol!(ucell::UCell, Orbs_Grid, DM, Density_Grid)
   
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPI_Hoffset = ucell.system_grid.MPI_Hoffset

    DM1 = DM[1]
    DM2 = DM[2]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch

    fill!(Density_Grid[1], 0.0)
    fill!(Density_Grid[2], 0.0)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        block_size = NO0*NO1
        dm_index = MPI_Hoffset[loop]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            dm_index += 1
            density_matrix_scratch[jst, ist, 1] = DM1[dm_index]
            density_matrix_scratch[jst, ist, 2] = DM2[dm_index]
        end
        @inbounds for spin = 1:2, xyz = 1:GridN_Atom[atom]
            density_scratch[xyz, spin] = 0.0
        end
        orbitals1 = atom_matrix(Orbs_Grid, atom)
        orbitals2 = atom_matrix(Orbs_Grid, jatom)
        _Calc_Density_Blocked!(density_scratch, 2, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbitals1, orbitals2, density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            global_grid = GridListAtom[atom][xyz] + 1
            Density_Grid[1][global_grid] += density_scratch[xyz, 1]
            Density_Grid[2][global_grid] += density_scratch[xyz, 2]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[2], MPI.SUM, comm)
end


function Set_Density_Grid_nc!(ucell::UCell, Orbs_Grid, DM, Density_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    MPI_atom = ucell.system_grid.MPI_atom
    MPI_natn = ucell.system_grid.MPI_natn
    MPI_size = ucell.system_grid.MPI_size
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
    Total_NumOrbs = ucell.system_grid.Total_NumOrbs
    MPI_Hoffset = ucell.system_grid.MPI_Hoffset

    DM1 = DM[1]
    DM2 = DM[2]
    DM3 = DM[3]
    DM4 = DM[4]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch

    fill!(Density_Grid[1], 0.0)
    fill!(Density_Grid[2], 0.0)
    fill!(Density_Grid[3], 0.0)
    fill!(Density_Grid[4], 0.0)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        
        block_size = NO0*NO1
        dm_index = MPI_Hoffset[loop]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            dm_index += 1
            density_matrix_scratch[jst, ist, 1] = DM1[dm_index]
            density_matrix_scratch[jst, ist, 2] = DM2[dm_index]
            density_matrix_scratch[jst, ist, 3] = DM3[dm_index]
            density_matrix_scratch[jst, ist, 4] = DM4[dm_index]
        end
        @inbounds for spin = 1:4, xyz = 1:GridN_Atom[atom]
            density_scratch[xyz, spin] = 0.0
        end
        orbitals1 = atom_matrix(Orbs_Grid, atom)
        orbitals2 = atom_matrix(Orbs_Grid, jatom)
        _Calc_Density_Blocked!(density_scratch, 4, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbitals1, orbitals2, density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            global_grid = GridListAtom[atom][xyz] + 1
            Density_Grid[1][global_grid] += density_scratch[xyz, 1]
            Density_Grid[2][global_grid] += density_scratch[xyz, 2]
            Density_Grid[3][global_grid] += density_scratch[xyz, 3]
            Density_Grid[4][global_grid] += density_scratch[xyz, 4]
        end
    end

    MPI.Allreduce!(Density_Grid[1], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[2], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[3], MPI.SUM, comm)
    MPI.Allreduce!(Density_Grid[4], MPI.SUM, comm)
end


function diagonalize_nc_density!(Density_Grid)

    Ngrid = length(Density_Grid[1])

    @inbounds for i = 1:Ngrid

        Re11 = Density_Grid[1][i]
        Re22 = Density_Grid[2][i]
        Re12 = Density_Grid[3][i]
        Im12 = Density_Grid[4][i]

        Nup, Ndown, theta, phi = EulerAngle_Spin(Re11, Re22, Re12, Im12)

        Density_Grid[1][i] = Nup
        Density_Grid[2][i] = Ndown
        Density_Grid[3][i] = theta
        Density_Grid[4][i] = phi 
    end
end


function _Calc_Density_Blocked!(density, nspin, NO0, NO1, NumOLG,
                                GListTAtoms1, GListTAtoms2,
                                orbitals1::Matrix{Float64}, orbitals2::Matrix{Float64},
                                density_matrix::Array{Float64,3}, orbital_scratch)
    orbital_block1, orbital_block2, product_block = orbital_scratch
    block_size = size(orbital_block1, 2)

    for block_start = 1:block_size:NumOLG
        points = min(block_size, NumOLG - block_start + 1)
        @inbounds for point = 1:points
            overlap_index = block_start + point - 1
            grid1 = GListTAtoms1[overlap_index] + 1
            grid2 = GListTAtoms2[overlap_index] + 1
            @inbounds for orbital = 1:NO0
                orbital_block1[orbital, point] = orbitals1[orbital, grid1]
            end
            @inbounds for orbital = 1:NO1
                orbital_block2[orbital, point] = orbitals2[orbital, grid2]
            end
        end

        for spin = 1:nspin
            @views mul!(product_block[1:NO0, 1:points],
                        transpose(density_matrix[1:NO1, 1:NO0, spin]),
                        orbital_block2[1:NO1, 1:points])
            @inbounds for point = 1:points
                value = 0.0
                @inbounds for orbital = 1:NO0
                    value += orbital_block1[orbital, point]*product_block[orbital, point]
                end
                local_grid = GListTAtoms1[block_start + point - 1] + 1
                density[local_grid, spin] += value
            end
        end
    end
    return nothing
end
