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
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]
    DM1 = DM[1]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch
    orbital_data = Orbs_Grid.data

    DMsum = 0
    fill!(Density_Grid[1], 0.0)

    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        block_size = NO0*NO1
        dm_index = DMnum + DMsum
        @inbounds for ist = 1:NO0, jst = 1:NO1
            dm_index += 1
            density_matrix_scratch[jst, ist, 1] = DM1[dm_index]
        end
        @inbounds for xyz = 1:GridN_Atom[atom]
            density_scratch[xyz, 1] = 0.0
        end
        _Calc_Density_Blocked!(density_scratch, 1, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbital_data[atom], orbital_data[jatom], density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            Density_Grid[1][GridListAtom[atom][xyz] + 1] += density_scratch[xyz, 1]
        end
        DMsum += block_size
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
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]
    DM1 = DM[1]
    DM2 = DM[2]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch
    orbital_data = Orbs_Grid.data

    DMsum = 0
    fill!(Density_Grid[1], 0.0)
    fill!(Density_Grid[2], 0.0)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        block_size = NO0*NO1
        dm_index = DMnum + DMsum
        @inbounds for ist = 1:NO0, jst = 1:NO1
            dm_index += 1
            density_matrix_scratch[jst, ist, 1] = DM1[dm_index]
            density_matrix_scratch[jst, ist, 2] = DM2[dm_index]
        end
        @inbounds for spin = 1:2, xyz = 1:GridN_Atom[atom]
            density_scratch[xyz, spin] = 0.0
        end
        _Calc_Density_Blocked!(density_scratch, 2, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbital_data[atom], orbital_data[jatom], density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            global_grid = GridListAtom[atom][xyz] + 1
            Density_Grid[1][global_grid] += density_scratch[xyz, 1]
            Density_Grid[2][global_grid] += density_scratch[xyz, 2]
        end
        DMsum += block_size
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
    MPHks = ucell.system_grid.MPHks
    
    DMnum = MPHks[myrank+1]
    DM1 = DM[1]
    DM2 = DM[2]
    DM3 = DM[3]
    DM4 = DM[4]
    density_scratch = ucell.density_scratch
    density_matrix_scratch = ucell.density_matrix_scratch
    density_orbital_scratch = ucell.density_orbital_scratch
    orbital_data = Orbs_Grid.data

    DMsum = 0
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
        dm_index = DMnum + DMsum
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
        _Calc_Density_Blocked!(density_scratch, 4, NO0, NO1, MPI_NumOLG[loop], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], orbital_data[atom], orbital_data[jatom], density_matrix_scratch, density_orbital_scratch)
        @inbounds for xyz = 1:GridN_Atom[atom]
            global_grid = GridListAtom[atom][xyz] + 1
            Density_Grid[1][global_grid] += density_scratch[xyz, 1]
            Density_Grid[2][global_grid] += density_scratch[xyz, 2]
            Density_Grid[3][global_grid] += density_scratch[xyz, 3]
            Density_Grid[4][global_grid] += density_scratch[xyz, 4]
        end
        DMsum += block_size
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


@inline _density_matrix_value(matrix, j::Int, i::Int) = @inbounds matrix[j, i]
@inline _density_matrix_value(matrix::AbstractArray{<:Any,3}, j::Int, i::Int) = @inbounds matrix[j, i, 1]


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
            for orbital = 1:NO0
                orbital_block1[orbital, point] = orbitals1[orbital, grid1]
            end
            for orbital = 1:NO1
                orbital_block2[orbital, point] = orbitals2[orbital, grid2]
            end
        end

        for spin = 1:nspin
            @views mul!(product_block[1:NO0, 1:points],
                        transpose(density_matrix[1:NO1, 1:NO0, spin]),
                        orbital_block2[1:NO1, 1:points])
            @inbounds for point = 1:points
                value = 0.0
                @simd for orbital = 1:NO0
                    value += orbital_block1[orbital, point]*product_block[orbital, point]
                end
                local_grid = GListTAtoms1[block_start + point - 1] + 1
                density[local_grid, spin] += value
            end
        end
    end
    return nothing
end


function _Calc_Den16_nonpol!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:16:NumOLG-15

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1
        Nc2 = GListTAtoms1[Nog+2]+1
        Nc3 = GListTAtoms1[Nog+3]+1
        Nc4 = GListTAtoms1[Nog+4]+1
        Nc5 = GListTAtoms1[Nog+5]+1
        Nc6 = GListTAtoms1[Nog+6]+1
        Nc7 = GListTAtoms1[Nog+7]+1
        Nc8 = GListTAtoms1[Nog+8]+1
        Nc9 = GListTAtoms1[Nog+9]+1
        Nc10 = GListTAtoms1[Nog+10]+1
        Nc11 = GListTAtoms1[Nog+11]+1
        Nc12 = GListTAtoms1[Nog+12]+1
        Nc13 = GListTAtoms1[Nog+13]+1
        Nc14 = GListTAtoms1[Nog+14]+1
        Nc15 = GListTAtoms1[Nog+15]+1

        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1
        Nh2 = GListTAtoms2[Nog+2]+1
        Nh3 = GListTAtoms2[Nog+3]+1
        Nh4 = GListTAtoms2[Nog+4]+1
        Nh5 = GListTAtoms2[Nog+5]+1
        Nh6 = GListTAtoms2[Nog+6]+1
        Nh7 = GListTAtoms2[Nog+7]+1
        Nh8 = GListTAtoms2[Nog+8]+1
        Nh9 = GListTAtoms2[Nog+9]+1
        Nh10 = GListTAtoms2[Nog+10]+1
        Nh11 = GListTAtoms2[Nog+11]+1
        Nh12 = GListTAtoms2[Nog+12]+1
        Nh13 = GListTAtoms2[Nog+13]+1
        Nh14 = GListTAtoms2[Nog+14]+1
        Nh15 = GListTAtoms2[Nog+15]+1
        
        Sum0 = 0.0
        Sum1 = 0.0
        Sum2 = 0.0
        Sum3 = 0.0
        Sum4 = 0.0
        Sum5 = 0.0
        Sum6 = 0.0
        Sum7 = 0.0
        Sum8 = 0.0
        Sum9 = 0.0
        Sum10 = 0.0
        Sum11 = 0.0
        Sum12 = 0.0
        Sum13 = 0.0
        Sum14 = 0.0
        Sum15 = 0.0
        for ist = 1:NO0
            temp0 = 0.0
            temp1 = 0.0
            temp2 = 0.0
            temp3 = 0.0
            temp4 = 0.0
            temp5 = 0.0
            temp6 = 0.0
            temp7 = 0.0
            temp8 = 0.0
            temp9 = 0.0
            temp10 = 0.0
            temp11 = 0.0
            temp12 = 0.0
            temp13 = 0.0
            temp14 = 0.0
            temp15 = 0.0
            @inbounds for jst = 1:NO1
                tmp = _density_matrix_value(DM, jst, ist)
                temp0  += Orbs_Grid2[Nh0][jst]*tmp
                temp1  += Orbs_Grid2[Nh1][jst]*tmp
                temp2  += Orbs_Grid2[Nh2][jst]*tmp
                temp3  += Orbs_Grid2[Nh3][jst]*tmp
                temp4  += Orbs_Grid2[Nh4][jst]*tmp
                temp5  += Orbs_Grid2[Nh5][jst]*tmp
                temp6  += Orbs_Grid2[Nh6][jst]*tmp
                temp7  += Orbs_Grid2[Nh7][jst]*tmp
                temp8  += Orbs_Grid2[Nh8][jst]*tmp
                temp9  += Orbs_Grid2[Nh9][jst]*tmp
                temp10 += Orbs_Grid2[Nh10][jst]*tmp
                temp11 += Orbs_Grid2[Nh11][jst]*tmp
                temp12 += Orbs_Grid2[Nh12][jst]*tmp
                temp13 += Orbs_Grid2[Nh13][jst]*tmp
                temp14 += Orbs_Grid2[Nh14][jst]*tmp
                temp15 += Orbs_Grid2[Nh15][jst]*tmp
            end
            Sum0  += Orbs_Grid1[Nc0][ist]*temp0
            Sum1  += Orbs_Grid1[Nc1][ist]*temp1
            Sum2  += Orbs_Grid1[Nc2][ist]*temp2
            Sum3  += Orbs_Grid1[Nc3][ist]*temp3
            Sum4  += Orbs_Grid1[Nc4][ist]*temp4
            Sum5  += Orbs_Grid1[Nc5][ist]*temp5
            Sum6  += Orbs_Grid1[Nc6][ist]*temp6
            Sum7  += Orbs_Grid1[Nc7][ist]*temp7
            Sum8  += Orbs_Grid1[Nc8][ist]*temp8
            Sum9  += Orbs_Grid1[Nc9][ist]*temp9
            Sum10 += Orbs_Grid1[Nc10][ist]*temp10
            Sum11 += Orbs_Grid1[Nc11][ist]*temp11
            Sum12 += Orbs_Grid1[Nc12][ist]*temp12
            Sum13 += Orbs_Grid1[Nc13][ist]*temp13
            Sum14 += Orbs_Grid1[Nc14][ist]*temp14
            Sum15 += Orbs_Grid1[Nc15][ist]*temp15
        end
        ai_tempDGs[Nc0] += Sum0
        ai_tempDGs[Nc1] += Sum1
        ai_tempDGs[Nc2] += Sum2
        ai_tempDGs[Nc3] += Sum3
        ai_tempDGs[Nc4] += Sum4
        ai_tempDGs[Nc5] += Sum5
        ai_tempDGs[Nc6] += Sum6
        ai_tempDGs[Nc7] += Sum7
        ai_tempDGs[Nc8] += Sum8
        ai_tempDGs[Nc9] += Sum9
        ai_tempDGs[Nc10] += Sum10
        ai_tempDGs[Nc11] += Sum11
        ai_tempDGs[Nc12] += Sum12
        ai_tempDGs[Nc13] += Sum13
        ai_tempDGs[Nc14] += Sum14
        ai_tempDGs[Nc15] += Sum15
    end


    Nog1 = 16*div(NumOLG, 16)
    rem_NumOLG = rem(NumOLG, 16)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum = 0.0
        for ist = 1:NO0
            temp = 0.0
            @inbounds for jst = 1:NO1
                temp += Orbs_Grid2[Nh][jst]*_density_matrix_value(DM, jst, ist)
            end
            Sum += Orbs_Grid1[Nc][ist]*temp
        end
        ai_tempDGs[Nc] += Sum
    end
end


function _Calc_Den16_pol!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:16:NumOLG-15

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1
        Nc2 = GListTAtoms1[Nog+2]+1
        Nc3 = GListTAtoms1[Nog+3]+1
        Nc4 = GListTAtoms1[Nog+4]+1
        Nc5 = GListTAtoms1[Nog+5]+1
        Nc6 = GListTAtoms1[Nog+6]+1
        Nc7 = GListTAtoms1[Nog+7]+1
        Nc8 = GListTAtoms1[Nog+8]+1
        Nc9 = GListTAtoms1[Nog+9]+1
        Nc10 = GListTAtoms1[Nog+10]+1
        Nc11 = GListTAtoms1[Nog+11]+1
        Nc12 = GListTAtoms1[Nog+12]+1
        Nc13 = GListTAtoms1[Nog+13]+1
        Nc14 = GListTAtoms1[Nog+14]+1
        Nc15 = GListTAtoms1[Nog+15]+1

        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1
        Nh2 = GListTAtoms2[Nog+2]+1
        Nh3 = GListTAtoms2[Nog+3]+1
        Nh4 = GListTAtoms2[Nog+4]+1
        Nh5 = GListTAtoms2[Nog+5]+1
        Nh6 = GListTAtoms2[Nog+6]+1
        Nh7 = GListTAtoms2[Nog+7]+1
        Nh8 = GListTAtoms2[Nog+8]+1
        Nh9 = GListTAtoms2[Nog+9]+1
        Nh10 = GListTAtoms2[Nog+10]+1
        Nh11 = GListTAtoms2[Nog+11]+1
        Nh12 = GListTAtoms2[Nog+12]+1
        Nh13 = GListTAtoms2[Nog+13]+1
        Nh14 = GListTAtoms2[Nog+14]+1
        Nh15 = GListTAtoms2[Nog+15]+1

        Sum0_up = 0.0
        Sum1_up = 0.0
        Sum2_up = 0.0
        Sum3_up = 0.0
        Sum4_up = 0.0
        Sum5_up = 0.0
        Sum6_up = 0.0
        Sum7_up = 0.0
        Sum8_up = 0.0
        Sum9_up = 0.0
        Sum10_up = 0.0
        Sum11_up = 0.0
        Sum12_up = 0.0
        Sum13_up = 0.0
        Sum14_up = 0.0
        Sum15_up = 0.0

        Sum0_dn = 0.0
        Sum1_dn = 0.0
        Sum2_dn = 0.0
        Sum3_dn = 0.0
        Sum4_dn = 0.0
        Sum5_dn = 0.0
        Sum6_dn = 0.0
        Sum7_dn = 0.0
        Sum8_dn = 0.0
        Sum9_dn = 0.0
        Sum10_dn = 0.0
        Sum11_dn = 0.0
        Sum12_dn = 0.0
        Sum13_dn = 0.0
        Sum14_dn = 0.0
        Sum15_dn = 0.0

        for ist = 1:NO0
            temp0_up = 0.0
            temp1_up = 0.0
            temp2_up = 0.0
            temp3_up = 0.0
            temp4_up = 0.0
            temp5_up = 0.0
            temp6_up = 0.0
            temp7_up = 0.0
            temp8_up = 0.0
            temp9_up = 0.0
            temp10_up = 0.0
            temp11_up = 0.0
            temp12_up = 0.0
            temp13_up = 0.0
            temp14_up = 0.0
            temp15_up = 0.0

            temp0_dn = 0.0
            temp1_dn = 0.0
            temp2_dn = 0.0
            temp3_dn = 0.0
            temp4_dn = 0.0
            temp5_dn = 0.0
            temp6_dn = 0.0
            temp7_dn = 0.0
            temp8_dn = 0.0
            temp9_dn = 0.0
            temp10_dn = 0.0
            temp11_dn = 0.0
            temp12_dn = 0.0
            temp13_dn = 0.0
            temp14_dn = 0.0
            temp15_dn = 0.0

            @inbounds for jst = 1:NO1
                orbs2_0  = Orbs_Grid2[Nh0][jst]
                orbs2_1  = Orbs_Grid2[Nh1][jst]
                orbs2_2  = Orbs_Grid2[Nh2][jst]
                orbs2_3  = Orbs_Grid2[Nh3][jst]
                orbs2_4  = Orbs_Grid2[Nh4][jst]
                orbs2_5  = Orbs_Grid2[Nh5][jst]
                orbs2_6  = Orbs_Grid2[Nh6][jst]
                orbs2_7  = Orbs_Grid2[Nh7][jst]
                orbs2_8  = Orbs_Grid2[Nh8][jst]
                orbs2_9  = Orbs_Grid2[Nh9][jst]
                orbs2_10 = Orbs_Grid2[Nh10][jst]
                orbs2_11 = Orbs_Grid2[Nh11][jst]
                orbs2_12 = Orbs_Grid2[Nh12][jst]
                orbs2_13 = Orbs_Grid2[Nh13][jst]
                orbs2_14 = Orbs_Grid2[Nh14][jst]
                orbs2_15 = Orbs_Grid2[Nh15][jst]

                DM_up = DM[jst,ist,1]
                DM_dn = DM[jst,ist,2]
                
                temp0_up += orbs2_0*DM_up
                temp1_up += orbs2_1*DM_up
                temp2_up += orbs2_2*DM_up
                temp3_up += orbs2_3*DM_up
                temp4_up += orbs2_4*DM_up
                temp5_up += orbs2_5*DM_up
                temp6_up += orbs2_6*DM_up
                temp7_up += orbs2_7*DM_up
                temp8_up += orbs2_8*DM_up
                temp9_up += orbs2_9*DM_up
                temp10_up += orbs2_10*DM_up
                temp11_up += orbs2_11*DM_up
                temp12_up += orbs2_12*DM_up
                temp13_up += orbs2_13*DM_up
                temp14_up += orbs2_14*DM_up
                temp15_up += orbs2_15*DM_up

                temp0_dn += orbs2_0*DM_dn
                temp1_dn += orbs2_1*DM_dn
                temp2_dn += orbs2_2*DM_dn
                temp3_dn += orbs2_3*DM_dn
                temp4_dn += orbs2_4*DM_dn
                temp5_dn += orbs2_5*DM_dn
                temp6_dn += orbs2_6*DM_dn
                temp7_dn += orbs2_7*DM_dn
                temp8_dn += orbs2_8*DM_dn
                temp9_dn += orbs2_9*DM_dn
                temp10_dn += orbs2_10*DM_dn
                temp11_dn += orbs2_11*DM_dn
                temp12_dn += orbs2_12*DM_dn
                temp13_dn += orbs2_13*DM_dn
                temp14_dn += orbs2_14*DM_dn
                temp15_dn += orbs2_15*DM_dn
            end

            orbs1_0  = Orbs_Grid1[Nc0][ist]
            orbs1_1  = Orbs_Grid1[Nc1][ist]
            orbs1_2  = Orbs_Grid1[Nc2][ist]
            orbs1_3  = Orbs_Grid1[Nc3][ist]
            orbs1_4  = Orbs_Grid1[Nc4][ist]
            orbs1_5  = Orbs_Grid1[Nc5][ist]
            orbs1_6  = Orbs_Grid1[Nc6][ist]
            orbs1_7  = Orbs_Grid1[Nc7][ist]
            orbs1_8  = Orbs_Grid1[Nc8][ist]
            orbs1_9  = Orbs_Grid1[Nc9][ist]
            orbs1_10 = Orbs_Grid1[Nc10][ist]
            orbs1_11 = Orbs_Grid1[Nc11][ist]
            orbs1_12 = Orbs_Grid1[Nc12][ist]
            orbs1_13 = Orbs_Grid1[Nc13][ist]
            orbs1_14 = Orbs_Grid1[Nc14][ist]
            orbs1_15 = Orbs_Grid1[Nc15][ist]

            Sum0_up += orbs1_0*temp0_up
            Sum1_up += orbs1_1*temp1_up
            Sum2_up += orbs1_2*temp2_up
            Sum3_up += orbs1_3*temp3_up
            Sum4_up += orbs1_4*temp4_up
            Sum5_up += orbs1_5*temp5_up
            Sum6_up += orbs1_6*temp6_up
            Sum7_up += orbs1_7*temp7_up
            Sum8_up += orbs1_8*temp8_up
            Sum9_up += orbs1_9*temp9_up
            Sum10_up += orbs1_10*temp10_up
            Sum11_up += orbs1_11*temp11_up
            Sum12_up += orbs1_12*temp12_up
            Sum13_up += orbs1_13*temp13_up
            Sum14_up += orbs1_14*temp14_up
            Sum15_up += orbs1_15*temp15_up

            Sum0_dn += orbs1_0*temp0_dn
            Sum1_dn += orbs1_1*temp1_dn
            Sum2_dn += orbs1_2*temp2_dn
            Sum3_dn += orbs1_3*temp3_dn
            Sum4_dn += orbs1_4*temp4_dn
            Sum5_dn += orbs1_5*temp5_dn
            Sum6_dn += orbs1_6*temp6_dn
            Sum7_dn += orbs1_7*temp7_dn
            Sum8_dn += orbs1_8*temp8_dn
            Sum9_dn += orbs1_9*temp9_dn
            Sum10_dn += orbs1_10*temp10_dn
            Sum11_dn += orbs1_11*temp11_dn
            Sum12_dn += orbs1_12*temp12_dn
            Sum13_dn += orbs1_13*temp13_dn
            Sum14_dn += orbs1_14*temp14_dn
            Sum15_dn += orbs1_15*temp15_dn
        end
        ai_tempDGs[Nc0,1] += Sum0_up
        ai_tempDGs[Nc1,1] += Sum1_up
        ai_tempDGs[Nc2,1] += Sum2_up
        ai_tempDGs[Nc3,1] += Sum3_up
        ai_tempDGs[Nc4,1] += Sum4_up
        ai_tempDGs[Nc5,1] += Sum5_up
        ai_tempDGs[Nc6,1] += Sum6_up
        ai_tempDGs[Nc7,1] += Sum7_up
        ai_tempDGs[Nc8,1] += Sum8_up
        ai_tempDGs[Nc9,1] += Sum9_up
        ai_tempDGs[Nc10,1] += Sum10_up
        ai_tempDGs[Nc11,1] += Sum11_up
        ai_tempDGs[Nc12,1] += Sum12_up
        ai_tempDGs[Nc13,1] += Sum13_up
        ai_tempDGs[Nc14,1] += Sum14_up
        ai_tempDGs[Nc15,1] += Sum15_up

        ai_tempDGs[Nc0,2] += Sum0_dn
        ai_tempDGs[Nc1,2] += Sum1_dn
        ai_tempDGs[Nc2,2] += Sum2_dn
        ai_tempDGs[Nc3,2] += Sum3_dn
        ai_tempDGs[Nc4,2] += Sum4_dn
        ai_tempDGs[Nc5,2] += Sum5_dn
        ai_tempDGs[Nc6,2] += Sum6_dn
        ai_tempDGs[Nc7,2] += Sum7_dn
        ai_tempDGs[Nc8,2] += Sum8_dn
        ai_tempDGs[Nc9,2] += Sum9_dn
        ai_tempDGs[Nc10,2] += Sum10_dn
        ai_tempDGs[Nc11,2] += Sum11_dn
        ai_tempDGs[Nc12,2] += Sum12_dn
        ai_tempDGs[Nc13,2] += Sum13_dn
        ai_tempDGs[Nc14,2] += Sum14_dn
        ai_tempDGs[Nc15,2] += Sum15_dn
    end


    Nog1 = 16*div(NumOLG, 16)
    rem_NumOLG = rem(NumOLG, 16)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum_up = 0.0
        Sum_dn = 0.0
        for ist = 1:NO0
            temp_up = 0.0
            temp_dn = 0.0
            @inbounds for jst = 1:NO1
                temp_up += Orbs_Grid2[Nh][jst]*DM[jst,ist,1]
                temp_dn += Orbs_Grid2[Nh][jst]*DM[jst,ist,2]
            end
            Sum_up += Orbs_Grid1[Nc][ist]*temp_up
            Sum_dn += Orbs_Grid1[Nc][ist]*temp_dn
        end
        ai_tempDGs[Nc,1] += Sum_up
        ai_tempDGs[Nc,2] += Sum_dn
    end
end


function _Calc_Den16_nc!(ai_tempDGs, NO0, NO1, NumOLG, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, DM)

    for Nog = 1:16:NumOLG-15

        Nc0 = GListTAtoms1[Nog]+1
        Nc1 = GListTAtoms1[Nog+1]+1
        Nc2 = GListTAtoms1[Nog+2]+1
        Nc3 = GListTAtoms1[Nog+3]+1
        Nc4 = GListTAtoms1[Nog+4]+1
        Nc5 = GListTAtoms1[Nog+5]+1
        Nc6 = GListTAtoms1[Nog+6]+1
        Nc7 = GListTAtoms1[Nog+7]+1
        Nc8 = GListTAtoms1[Nog+8]+1
        Nc9 = GListTAtoms1[Nog+9]+1
        Nc10 = GListTAtoms1[Nog+10]+1
        Nc11 = GListTAtoms1[Nog+11]+1
        Nc12 = GListTAtoms1[Nog+12]+1
        Nc13 = GListTAtoms1[Nog+13]+1
        Nc14 = GListTAtoms1[Nog+14]+1
        Nc15 = GListTAtoms1[Nog+15]+1

        Nh0 = GListTAtoms2[Nog]+1
        Nh1 = GListTAtoms2[Nog+1]+1
        Nh2 = GListTAtoms2[Nog+2]+1
        Nh3 = GListTAtoms2[Nog+3]+1
        Nh4 = GListTAtoms2[Nog+4]+1
        Nh5 = GListTAtoms2[Nog+5]+1
        Nh6 = GListTAtoms2[Nog+6]+1
        Nh7 = GListTAtoms2[Nog+7]+1
        Nh8 = GListTAtoms2[Nog+8]+1
        Nh9 = GListTAtoms2[Nog+9]+1
        Nh10 = GListTAtoms2[Nog+10]+1
        Nh11 = GListTAtoms2[Nog+11]+1
        Nh12 = GListTAtoms2[Nog+12]+1
        Nh13 = GListTAtoms2[Nog+13]+1
        Nh14 = GListTAtoms2[Nog+14]+1
        Nh15 = GListTAtoms2[Nog+15]+1
        
        Sum0_uu = 0.0
        Sum1_uu = 0.0
        Sum2_uu = 0.0
        Sum3_uu = 0.0
        Sum4_uu = 0.0
        Sum5_uu = 0.0
        Sum6_uu = 0.0
        Sum7_uu = 0.0
        Sum8_uu = 0.0
        Sum9_uu = 0.0
        Sum10_uu = 0.0
        Sum11_uu = 0.0
        Sum12_uu = 0.0
        Sum13_uu = 0.0
        Sum14_uu = 0.0
        Sum15_uu = 0.0

        Sum0_dd = 0.0
        Sum1_dd = 0.0
        Sum2_dd = 0.0
        Sum3_dd = 0.0
        Sum4_dd = 0.0
        Sum5_dd = 0.0
        Sum6_dd = 0.0
        Sum7_dd = 0.0
        Sum8_dd = 0.0
        Sum9_dd = 0.0
        Sum10_dd = 0.0
        Sum11_dd = 0.0
        Sum12_dd = 0.0
        Sum13_dd = 0.0
        Sum14_dd = 0.0
        Sum15_dd = 0.0

        Sum0_ud_r = 0.0
        Sum1_ud_r = 0.0
        Sum2_ud_r = 0.0
        Sum3_ud_r = 0.0
        Sum4_ud_r = 0.0
        Sum5_ud_r = 0.0
        Sum6_ud_r = 0.0
        Sum7_ud_r = 0.0
        Sum8_ud_r = 0.0
        Sum9_ud_r = 0.0
        Sum10_ud_r = 0.0
        Sum11_ud_r = 0.0
        Sum12_ud_r = 0.0
        Sum13_ud_r = 0.0
        Sum14_ud_r = 0.0
        Sum15_ud_r = 0.0

        Sum0_ud_i = 0.0
        Sum1_ud_i = 0.0
        Sum2_ud_i = 0.0
        Sum3_ud_i = 0.0
        Sum4_ud_i = 0.0
        Sum5_ud_i = 0.0
        Sum6_ud_i = 0.0
        Sum7_ud_i = 0.0
        Sum8_ud_i = 0.0
        Sum9_ud_i = 0.0
        Sum10_ud_i = 0.0
        Sum11_ud_i = 0.0
        Sum12_ud_i = 0.0
        Sum13_ud_i = 0.0
        Sum14_ud_i = 0.0
        Sum15_ud_i = 0.0

        for ist = 1:NO0
            temp0_uu = 0.0
            temp1_uu = 0.0
            temp2_uu = 0.0
            temp3_uu = 0.0
            temp4_uu = 0.0
            temp5_uu = 0.0
            temp6_uu = 0.0
            temp7_uu = 0.0
            temp8_uu = 0.0
            temp9_uu = 0.0
            temp10_uu = 0.0
            temp11_uu = 0.0
            temp12_uu = 0.0
            temp13_uu = 0.0
            temp14_uu = 0.0
            temp15_uu = 0.0

            temp0_dd = 0.0
            temp1_dd = 0.0
            temp2_dd = 0.0
            temp3_dd = 0.0
            temp4_dd = 0.0
            temp5_dd = 0.0
            temp6_dd = 0.0
            temp7_dd = 0.0
            temp8_dd = 0.0
            temp9_dd = 0.0
            temp10_dd = 0.0
            temp11_dd = 0.0
            temp12_dd = 0.0
            temp13_dd = 0.0
            temp14_dd = 0.0
            temp15_dd = 0.0

            temp0_ud_r = 0.0
            temp1_ud_r = 0.0
            temp2_ud_r = 0.0
            temp3_ud_r = 0.0
            temp4_ud_r = 0.0
            temp5_ud_r = 0.0
            temp6_ud_r = 0.0
            temp7_ud_r = 0.0
            temp8_ud_r = 0.0
            temp9_ud_r = 0.0
            temp10_ud_r = 0.0
            temp11_ud_r = 0.0
            temp12_ud_r = 0.0
            temp13_ud_r = 0.0
            temp14_ud_r = 0.0
            temp15_ud_r = 0.0

            temp0_ud_i = 0.0
            temp1_ud_i = 0.0
            temp2_ud_i = 0.0
            temp3_ud_i = 0.0
            temp4_ud_i = 0.0
            temp5_ud_i = 0.0
            temp6_ud_i = 0.0
            temp7_ud_i = 0.0
            temp8_ud_i = 0.0
            temp9_ud_i = 0.0
            temp10_ud_i = 0.0
            temp11_ud_i = 0.0
            temp12_ud_i = 0.0
            temp13_ud_i = 0.0
            temp14_ud_i = 0.0
            temp15_ud_i = 0.0

            @inbounds for jst = 1:NO1

                orbs2_0  = Orbs_Grid2[Nh0][jst]
                orbs2_1  = Orbs_Grid2[Nh1][jst]
                orbs2_2  = Orbs_Grid2[Nh2][jst]
                orbs2_3  = Orbs_Grid2[Nh3][jst]
                orbs2_4  = Orbs_Grid2[Nh4][jst]
                orbs2_5  = Orbs_Grid2[Nh5][jst]
                orbs2_6  = Orbs_Grid2[Nh6][jst]
                orbs2_7  = Orbs_Grid2[Nh7][jst]
                orbs2_8  = Orbs_Grid2[Nh8][jst]
                orbs2_9  = Orbs_Grid2[Nh9][jst]
                orbs2_10 = Orbs_Grid2[Nh10][jst]
                orbs2_11 = Orbs_Grid2[Nh11][jst]
                orbs2_12 = Orbs_Grid2[Nh12][jst]
                orbs2_13 = Orbs_Grid2[Nh13][jst]
                orbs2_14 = Orbs_Grid2[Nh14][jst]
                orbs2_15 = Orbs_Grid2[Nh15][jst]

                DM_uu = DM[jst,ist,1]
                DM_dd = DM[jst,ist,2]
                DM_ud_r = DM[jst,ist,3]
                DM_ud_i = DM[jst,ist,4]
                
                temp0_uu += orbs2_0*DM_uu
                temp1_uu += orbs2_1*DM_uu
                temp2_uu += orbs2_2*DM_uu
                temp3_uu += orbs2_3*DM_uu
                temp4_uu += orbs2_4*DM_uu
                temp5_uu += orbs2_5*DM_uu
                temp6_uu += orbs2_6*DM_uu
                temp7_uu += orbs2_7*DM_uu
                temp8_uu += orbs2_8*DM_uu
                temp9_uu += orbs2_9*DM_uu
                temp10_uu += orbs2_10*DM_uu
                temp11_uu += orbs2_11*DM_uu
                temp12_uu += orbs2_12*DM_uu
                temp13_uu += orbs2_13*DM_uu
                temp14_uu += orbs2_14*DM_uu
                temp15_uu += orbs2_15*DM_uu

                temp0_dd += orbs2_0*DM_dd
                temp1_dd += orbs2_1*DM_dd
                temp2_dd += orbs2_2*DM_dd
                temp3_dd += orbs2_3*DM_dd
                temp4_dd += orbs2_4*DM_dd
                temp5_dd += orbs2_5*DM_dd
                temp6_dd += orbs2_6*DM_dd
                temp7_dd += orbs2_7*DM_dd
                temp8_dd += orbs2_8*DM_dd
                temp9_dd += orbs2_9*DM_dd
                temp10_dd += orbs2_10*DM_dd
                temp11_dd += orbs2_11*DM_dd
                temp12_dd += orbs2_12*DM_dd
                temp13_dd += orbs2_13*DM_dd
                temp14_dd += orbs2_14*DM_dd
                temp15_dd += orbs2_15*DM_dd

                temp0_ud_r += orbs2_0*DM_ud_r
                temp1_ud_r += orbs2_1*DM_ud_r
                temp2_ud_r += orbs2_2*DM_ud_r
                temp3_ud_r += orbs2_3*DM_ud_r
                temp4_ud_r += orbs2_4*DM_ud_r
                temp5_ud_r += orbs2_5*DM_ud_r
                temp6_ud_r += orbs2_6*DM_ud_r
                temp7_ud_r += orbs2_7*DM_ud_r
                temp8_ud_r += orbs2_8*DM_ud_r
                temp9_ud_r += orbs2_9*DM_ud_r
                temp10_ud_r += orbs2_10*DM_ud_r
                temp11_ud_r += orbs2_11*DM_ud_r
                temp12_ud_r += orbs2_12*DM_ud_r
                temp13_ud_r += orbs2_13*DM_ud_r
                temp14_ud_r += orbs2_14*DM_ud_r
                temp15_ud_r += orbs2_15*DM_ud_r

                temp0_ud_i += orbs2_0*DM_ud_i
                temp1_ud_i += orbs2_1*DM_ud_i
                temp2_ud_i += orbs2_2*DM_ud_i
                temp3_ud_i += orbs2_3*DM_ud_i
                temp4_ud_i += orbs2_4*DM_ud_i
                temp5_ud_i += orbs2_5*DM_ud_i
                temp6_ud_i += orbs2_6*DM_ud_i
                temp7_ud_i += orbs2_7*DM_ud_i
                temp8_ud_i += orbs2_8*DM_ud_i
                temp9_ud_i += orbs2_9*DM_ud_i
                temp10_ud_i += orbs2_10*DM_ud_i
                temp11_ud_i += orbs2_11*DM_ud_i
                temp12_ud_i += orbs2_12*DM_ud_i
                temp13_ud_i += orbs2_13*DM_ud_i
                temp14_ud_i += orbs2_14*DM_ud_i
                temp15_ud_i += orbs2_15*DM_ud_i
            end

            orbs1_0  = Orbs_Grid1[Nc0][ist]
            orbs1_1  = Orbs_Grid1[Nc1][ist]
            orbs1_2  = Orbs_Grid1[Nc2][ist]
            orbs1_3  = Orbs_Grid1[Nc3][ist]
            orbs1_4  = Orbs_Grid1[Nc4][ist]
            orbs1_5  = Orbs_Grid1[Nc5][ist]
            orbs1_6  = Orbs_Grid1[Nc6][ist]
            orbs1_7  = Orbs_Grid1[Nc7][ist]
            orbs1_8  = Orbs_Grid1[Nc8][ist]
            orbs1_9  = Orbs_Grid1[Nc9][ist]
            orbs1_10 = Orbs_Grid1[Nc10][ist]
            orbs1_11 = Orbs_Grid1[Nc11][ist]
            orbs1_12 = Orbs_Grid1[Nc12][ist]
            orbs1_13 = Orbs_Grid1[Nc13][ist]
            orbs1_14 = Orbs_Grid1[Nc14][ist]
            orbs1_15 = Orbs_Grid1[Nc15][ist]

            Sum0_uu += orbs1_0*temp0_uu
            Sum1_uu += orbs1_1*temp1_uu
            Sum2_uu += orbs1_2*temp2_uu
            Sum3_uu += orbs1_3*temp3_uu
            Sum4_uu += orbs1_4*temp4_uu
            Sum5_uu += orbs1_5*temp5_uu
            Sum6_uu += orbs1_6*temp6_uu
            Sum7_uu += orbs1_7*temp7_uu
            Sum8_uu += orbs1_8*temp8_uu
            Sum9_uu += orbs1_9*temp9_uu
            Sum10_uu += orbs1_10*temp10_uu
            Sum11_uu += orbs1_11*temp11_uu
            Sum12_uu += orbs1_12*temp12_uu
            Sum13_uu += orbs1_13*temp13_uu
            Sum14_uu += orbs1_14*temp14_uu
            Sum15_uu += orbs1_15*temp15_uu

            Sum0_dd += orbs1_0*temp0_dd
            Sum1_dd += orbs1_1*temp1_dd
            Sum2_dd += orbs1_2*temp2_dd
            Sum3_dd += orbs1_3*temp3_dd
            Sum4_dd += orbs1_4*temp4_dd
            Sum5_dd += orbs1_5*temp5_dd
            Sum6_dd += orbs1_6*temp6_dd
            Sum7_dd += orbs1_7*temp7_dd
            Sum8_dd += orbs1_8*temp8_dd
            Sum9_dd += orbs1_9*temp9_dd
            Sum10_dd += orbs1_10*temp10_dd
            Sum11_dd += orbs1_11*temp11_dd
            Sum12_dd += orbs1_12*temp12_dd
            Sum13_dd += orbs1_13*temp13_dd
            Sum14_dd += orbs1_14*temp14_dd
            Sum15_dd += orbs1_15*temp15_dd

            Sum0_ud_r += orbs1_0*temp0_ud_r
            Sum1_ud_r += orbs1_1*temp1_ud_r
            Sum2_ud_r += orbs1_2*temp2_ud_r
            Sum3_ud_r += orbs1_3*temp3_ud_r
            Sum4_ud_r += orbs1_4*temp4_ud_r
            Sum5_ud_r += orbs1_5*temp5_ud_r
            Sum6_ud_r += orbs1_6*temp6_ud_r
            Sum7_ud_r += orbs1_7*temp7_ud_r
            Sum8_ud_r += orbs1_8*temp8_ud_r
            Sum9_ud_r += orbs1_9*temp9_ud_r
            Sum10_ud_r += orbs1_10*temp10_ud_r
            Sum11_ud_r += orbs1_11*temp11_ud_r
            Sum12_ud_r += orbs1_12*temp12_ud_r
            Sum13_ud_r += orbs1_13*temp13_ud_r
            Sum14_ud_r += orbs1_14*temp14_ud_r
            Sum15_ud_r += orbs1_15*temp15_ud_r

            Sum0_ud_i += orbs1_0*temp0_ud_i
            Sum1_ud_i += orbs1_1*temp1_ud_i
            Sum2_ud_i += orbs1_2*temp2_ud_i
            Sum3_ud_i += orbs1_3*temp3_ud_i
            Sum4_ud_i += orbs1_4*temp4_ud_i
            Sum5_ud_i += orbs1_5*temp5_ud_i
            Sum6_ud_i += orbs1_6*temp6_ud_i
            Sum7_ud_i += orbs1_7*temp7_ud_i
            Sum8_ud_i += orbs1_8*temp8_ud_i
            Sum9_ud_i += orbs1_9*temp9_ud_i
            Sum10_ud_i += orbs1_10*temp10_ud_i
            Sum11_ud_i += orbs1_11*temp11_ud_i
            Sum12_ud_i += orbs1_12*temp12_ud_i
            Sum13_ud_i += orbs1_13*temp13_ud_i
            Sum14_ud_i += orbs1_14*temp14_ud_i
            Sum15_ud_i += orbs1_15*temp15_ud_i
        end
        ai_tempDGs[Nc0,1] += Sum0_uu
        ai_tempDGs[Nc1,1] += Sum1_uu
        ai_tempDGs[Nc2,1] += Sum2_uu
        ai_tempDGs[Nc3,1] += Sum3_uu
        ai_tempDGs[Nc4,1] += Sum4_uu
        ai_tempDGs[Nc5,1] += Sum5_uu
        ai_tempDGs[Nc6,1] += Sum6_uu
        ai_tempDGs[Nc7,1] += Sum7_uu
        ai_tempDGs[Nc8,1] += Sum8_uu
        ai_tempDGs[Nc9,1] += Sum9_uu
        ai_tempDGs[Nc10,1] += Sum10_uu
        ai_tempDGs[Nc11,1] += Sum11_uu
        ai_tempDGs[Nc12,1] += Sum12_uu
        ai_tempDGs[Nc13,1] += Sum13_uu
        ai_tempDGs[Nc14,1] += Sum14_uu
        ai_tempDGs[Nc15,1] += Sum15_uu

        ai_tempDGs[Nc0,2] += Sum0_dd
        ai_tempDGs[Nc1,2] += Sum1_dd
        ai_tempDGs[Nc2,2] += Sum2_dd
        ai_tempDGs[Nc3,2] += Sum3_dd
        ai_tempDGs[Nc4,2] += Sum4_dd
        ai_tempDGs[Nc5,2] += Sum5_dd
        ai_tempDGs[Nc6,2] += Sum6_dd
        ai_tempDGs[Nc7,2] += Sum7_dd
        ai_tempDGs[Nc8,2] += Sum8_dd
        ai_tempDGs[Nc9,2] += Sum9_dd
        ai_tempDGs[Nc10,2] += Sum10_dd
        ai_tempDGs[Nc11,2] += Sum11_dd
        ai_tempDGs[Nc12,2] += Sum12_dd
        ai_tempDGs[Nc13,2] += Sum13_dd
        ai_tempDGs[Nc14,2] += Sum14_dd
        ai_tempDGs[Nc15,2] += Sum15_dd

        ai_tempDGs[Nc0,3] += Sum0_ud_r
        ai_tempDGs[Nc1,3] += Sum1_ud_r
        ai_tempDGs[Nc2,3] += Sum2_ud_r
        ai_tempDGs[Nc3,3] += Sum3_ud_r
        ai_tempDGs[Nc4,3] += Sum4_ud_r
        ai_tempDGs[Nc5,3] += Sum5_ud_r
        ai_tempDGs[Nc6,3] += Sum6_ud_r
        ai_tempDGs[Nc7,3] += Sum7_ud_r
        ai_tempDGs[Nc8,3] += Sum8_ud_r
        ai_tempDGs[Nc9,3] += Sum9_ud_r
        ai_tempDGs[Nc10,3] += Sum10_ud_r
        ai_tempDGs[Nc11,3] += Sum11_ud_r
        ai_tempDGs[Nc12,3] += Sum12_ud_r
        ai_tempDGs[Nc13,3] += Sum13_ud_r
        ai_tempDGs[Nc14,3] += Sum14_ud_r
        ai_tempDGs[Nc15,3] += Sum15_ud_r

        ai_tempDGs[Nc0,4] += Sum0_ud_i
        ai_tempDGs[Nc1,4] += Sum1_ud_i
        ai_tempDGs[Nc2,4] += Sum2_ud_i
        ai_tempDGs[Nc3,4] += Sum3_ud_i
        ai_tempDGs[Nc4,4] += Sum4_ud_i
        ai_tempDGs[Nc5,4] += Sum5_ud_i
        ai_tempDGs[Nc6,4] += Sum6_ud_i
        ai_tempDGs[Nc7,4] += Sum7_ud_i
        ai_tempDGs[Nc8,4] += Sum8_ud_i
        ai_tempDGs[Nc9,4] += Sum9_ud_i
        ai_tempDGs[Nc10,4] += Sum10_ud_i
        ai_tempDGs[Nc11,4] += Sum11_ud_i
        ai_tempDGs[Nc12,4] += Sum12_ud_i
        ai_tempDGs[Nc13,4] += Sum13_ud_i
        ai_tempDGs[Nc14,4] += Sum14_ud_i
        ai_tempDGs[Nc15,4] += Sum15_ud_i
    end


    Nog1 = 16*div(NumOLG, 16)
    rem_NumOLG = rem(NumOLG, 16)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
        Nh = GListTAtoms2[Nog+Nog1]+1

        Sum_uu = 0.0
        Sum_dd = 0.0
        Sum_ud_r = 0.0
        Sum_ud_i = 0.0
        for ist = 1:NO0
            temp_uu = 0.0
            temp_dd = 0.0
            temp_ud_r = 0.0
            temp_ud_i = 0.0
            @inbounds for jst = 1:NO1
                orbs2 = Orbs_Grid2[Nh][jst]
                temp_uu += orbs2*DM[jst,ist,1]
                temp_dd += orbs2*DM[jst,ist,2]
                temp_ud_r += orbs2*DM[jst,ist,3]
                temp_ud_i += orbs2*DM[jst,ist,4]
            end
            orbs1 = Orbs_Grid1[Nc][ist]
            Sum_uu += orbs1*temp_uu
            Sum_dd += orbs1*temp_dd
            Sum_ud_r += orbs1*temp_ud_r
            Sum_ud_i += orbs1*temp_ud_i
        end
        ai_tempDGs[Nc,1] += Sum_uu
        ai_tempDGs[Nc,2] += Sum_dd
        ai_tempDGs[Nc,3] += Sum_ud_r
        ai_tempDGs[Nc,4] += Sum_ud_i
    end
end
