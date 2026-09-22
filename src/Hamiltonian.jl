struct Hamiltonian
	SpinPol::String
	OLP::Vector{Float64}
	MPI_Hkin::Vector{Float64}
	MPI_HNL::Vector{Vector{Float64}}
	MPI_iHNL::Vector{Vector{Float64}}
	MPI_HVNA::Vector{Float64}
	MPI_NLPforce::Vector{Vector{Vector{Matrix{Float64}}}}
	MPI_DS_VNAforce::Vector{Vector{Matrix{Float64}}}
	MPI_HVNA2force::Vector{Vector{Float64}}
	MPI_HVNA3force::Vector{Vector{Float64}}
end


"""Zero-copy atom-pair block backed by the flattened MPI Hamiltonian arrays."""
struct _HamiltonianBlock{N,T,V<:AbstractVector{T}}
    data::NTuple{N,V}
    offset::Int
    nrows::Int
end

@inline _HamiltonianBlock(data::NTuple{N,V}, offset::Integer, nrows::Integer) where {N,T,V<:AbstractVector{T}} =
    _HamiltonianBlock{N,T,V}(data, Int(offset), Int(nrows))

@inline function Base.getindex(block::_HamiltonianBlock, j::Int, i::Int)
    @inbounds return block.data[1][block.offset + (i - 1)*block.nrows + j]
end
@inline function Base.setindex!(block::_HamiltonianBlock, value, j::Int, i::Int)
    @inbounds block.data[1][block.offset + (i - 1)*block.nrows + j] = value
    return value
end
@inline function Base.getindex(block::_HamiltonianBlock, j::Int, i::Int, spin::Int)
    @inbounds return block.data[spin][block.offset + (i - 1)*block.nrows + j]
end
@inline function Base.setindex!(block::_HamiltonianBlock, value, j::Int, i::Int, spin::Int)
    @inbounds block.data[spin][block.offset + (i - 1)*block.nrows + j] = value
    return value
end


@timeit timer "Hamiltonian" function Hamiltonian(cal_force::Bool, SpinPol::AbstractString, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)	
	
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	if SpinPol ∉ ("off", "on", "nc")
		error("please check SpinPol.")
	end

	Natom = system_grid.Natom
	atom2spe = system_grid.atom2spe
	Total_NumOrbs = system_grid.Total_NumOrbs
	MPI_atom = system_grid.MPI_atom
	MPI_natn = system_grid.MPI_natn
	MPI_size = system_grid.MPI_size
	MPI_Hsize = system_grid.MPI_Hsize
	myHsize = MPI_Hsize[myrank+1]
	Total_Hsize = system_grid.Total_Hsize

	
    # Overlap/Kinetic Matrix
	OLP = zeros(Float64, Total_Hsize)
    MPI_Hkin = zeros(Float64, myHsize)

	myrank == 0 && println("<Set_OLP_Kin>  Calculation of the overlap matrix")
	Set_OLP_Kin!(OLP, MPI_Hkin, pao, system_grid)



	Nspecies = length(pao)
	maxL = maximum([pao[spe].Spe_MaxL_Basis for spe = 1:Nspecies]) + BufferL_ProVNA
    VNATotal_Num = (maxL+1)^2 * maxM

	MPI_DS_VNAforce = Vector{Vector{Matrix{Float64}}}(undef, 4)
	for xyz = 1:4
		MPI_DS_VNAforce[xyz] = Vector{Matrix{Float64}}(undef, MPI_size)
		for loop = 1:MPI_size
			NO0 = Total_NumOrbs[MPI_atom[loop]]
			MPI_DS_VNAforce[xyz][loop] = zeros(Float64, NO0, VNATotal_Num)
		end
	end

	MPI_HVNA = zeros(Float64, myHsize)

	myHVNA2force = 0
	myHVNA3force = 0
    for loop = 1:MPI_size
		for _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:Total_NumOrbs[MPI_atom[loop]]
        	myHVNA2force += 1
		end
		for _ = 1:Total_NumOrbs[MPI_natn[loop]], _ = 1:Total_NumOrbs[MPI_natn[loop]]
        	myHVNA3force += 1
		end
    end

	MPI_HVNA2force = Vector{Vector{Float64}}(undef, 3)
    MPI_HVNA3force = Vector{Vector{Float64}}(undef, 3)
	for xyz = 1:3
		MPI_HVNA2force[xyz] = zeros(Float64, myHVNA2force)
		MPI_HVNA3force[xyz] = zeros(Float64, myHVNA3force)
	end

	myrank == 0 && println("<Set_ProExpn_VNA>  Calculation of the VNA projector matrix")
	Set_ProExpn_VNA!(MPI_DS_VNAforce, MPI_HVNA, MPI_HVNA2force, MPI_HVNA3force, pao, pspot, system_grid)
	if !cal_force
		MPI_DS_VNAforce = [[zeros(Float64,1,1)]]
		MPI_HVNA2force = [[0.0]]
		MPI_HVNA3force = [[0.0]]
	end




	VPS_j_Num = zeros(Int64, Natom)
    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot[spe].Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end

        VPS_j_Num[atom] = pspot[spe].VPS_j_dependency
        NLTotal_Num[atom] = tot
    end

    MPI_NLPforce = Vector{Vector{Vector{Matrix{Float64}}}}(undef, 4)
    for xyz = 1:4
        MPI_NLPforce[xyz] = Vector{Vector{Matrix{Float64}}}(undef, MPI_size)
        for loop = 1:MPI_size
			NO0 = Total_NumOrbs[MPI_atom[loop]]
			VPS_j_dependency = VPS_j_Num[MPI_natn[loop]]
			NO1 = NLTotal_Num[MPI_natn[loop]]
            MPI_NLPforce[xyz][loop] = Vector{Matrix{Float64}}(undef, VPS_j_dependency+1)
            for so = 1:VPS_j_dependency+1
                MPI_NLPforce[xyz][loop][so] = zeros(Float64, NO0, NO1)
            end
        end
    end


	if SpinPol ∈ ("off", "on")
		MPI_HNL = Vector{Vector{Float64}}(undef, 1)
		for spin = 1:1
			MPI_HNL[spin] = zeros(Float64, myHsize)
		end
		MPI_iHNL = [[1.0]]
	elseif SpinPol == "nc"
		MPI_HNL = Vector{Vector{Float64}}(undef, 3)
		MPI_iHNL = Vector{Vector{Float64}}(undef, 3)
		for spin = 1:3
			MPI_HNL[spin] = zeros(Float64, myHsize)
			MPI_iHNL[spin] = zeros(Float64, myHsize)
		end
	end
	    
	myrank == 0 && println("<Set_Nonlocal>  Calculation of the nonlocal matrix")
	Set_Nonlocal!(SpinPol, MPI_NLPforce, MPI_HNL, MPI_iHNL, pao, pspot, system_grid)
	if !cal_force
		MPI_NLPforce = [[[zeros(Float64, 1, 1)]]]
	end
	

	return Hamiltonian(SpinPol, OLP, MPI_Hkin, MPI_HNL, MPI_iHNL, MPI_HVNA, MPI_NLPforce, MPI_DS_VNAforce, MPI_HVNA2force, MPI_HVNA3force)
end



"""
Set_Hamiltonian!  
update Hamiltonian matrix using potentials and Orbs_Grid

Mandatory arguments:

- `Ham`: an instance of `Hamiltonian`
- `Orbs_Grid`: PAO Real Grid `ϕiα(r)`
- `Vpot_Grid`: `V(r)`, (periodic case when `V(r) = V(r+Rn)` satisfy, `V(r)` only in cell data (`Ngrid1*Ngrid2*Ngrid3`))
- `Hks`: calculate Kohn-Sham Hamiltonian matrix, `<ϕiα|Hks|ϕjβ>`
"""
function Set_Hamiltonian!(Ham::Hamiltonian, ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, MPI_Hks)
	
    comm = MPI.COMM_WORLD
	myrank = MPI.Comm_rank(comm)

	system_grid = ucell.system_grid
	myHsize = system_grid.MPI_Hsize[myrank+1]
	SpinPol = Ham.SpinPol
	MPI_Hkin = Ham.MPI_Hkin
	MPI_HVNA = Ham.MPI_HVNA
	MPI_HNL = Ham.MPI_HNL

	if SpinPol == "off"
		Calc_MatrixElements_dVH_Vxc_off!(ucell, Orbs_Grid, Vpot_Grid, MPI_Hks)
		Add_Hkin_HVNA_HNL!(MPI_Hks[1], MPI_Hkin, MPI_HVNA, MPI_HNL[1], myHsize)
	elseif SpinPol == "on"
		Calc_MatrixElements_dVH_Vxc_on!(ucell, Orbs_Grid, Vpot_Grid, MPI_Hks)
		Add_Hkin_HVNA_HNL!(MPI_Hks[1], MPI_Hkin, MPI_HVNA, MPI_HNL[1], myHsize)
		Add_Hkin_HVNA_HNL!(MPI_Hks[2], MPI_Hkin, MPI_HVNA, MPI_HNL[1], myHsize)
	elseif SpinPol == "nc"
		Calc_MatrixElements_dVH_Vxc_nc!(ucell, Orbs_Grid, Vpot_Grid, MPI_Hks)
		Add_Hkin_HVNA_HNL!(MPI_Hks[1], MPI_Hkin, MPI_HVNA, MPI_HNL[1], myHsize)
		Add_Hkin_HVNA_HNL!(MPI_Hks[2], MPI_Hkin, MPI_HVNA, MPI_HNL[2], myHsize)
		Add_HNL3!(MPI_Hks[3], MPI_HNL[3], myHsize)
	end
end


@inline function Add_HNL3!(Hks, HNL, myHsize)
	@inbounds for hst = 1:myHsize
		Hks[hst] += HNL[hst]
	end
end


@inline function Add_Hkin_HVNA_HNL!(Hks, Hkin, HVNA, HNL, myHsize)
	@inbounds for hst = 1:myHsize
		Hks[hst] += Hkin[hst] + HVNA[hst] + HNL[hst]
	end
end


@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_off!(
    ucell::UCell, orbitals_grid, potential_grid::Vector{Vector{Float64}}, local_hks)
    return _Calc_MatrixElements_dVH_Vxc!(Val(1), ucell, orbitals_grid,
                                         potential_grid, local_hks)
end

@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_on!(
    ucell::UCell, orbitals_grid, potential_grid::Vector{Vector{Float64}}, local_hks)
    return _Calc_MatrixElements_dVH_Vxc!(Val(2), ucell, orbitals_grid,
                                         potential_grid, local_hks)
end

@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_nc!(
    ucell::UCell, orbitals_grid, potential_grid::Vector{Vector{Float64}}, local_hks)
    return _Calc_MatrixElements_dVH_Vxc!(Val(4), ucell, orbitals_grid,
                                         potential_grid, local_hks)
end

"""Integrate all rank-local atom pairs for `N` potential channels."""
function _Calc_MatrixElements_dVH_Vxc!(::Val{N}, ucell::UCell,
                                       orbitals_grid, potential_grid,
                                       local_hks) where N
    grid = ucell.system_grid
    hks_channels = ntuple(channel -> local_hks[channel], Val(N))
    foreach(channel -> fill!(channel, 0.0), hks_channels)

    offset = 0
    for pair in eachindex(grid.MPI_atom)
        atom = grid.MPI_atom[pair]
        neighbor = grid.MPI_natn[pair]
        atom_orbitals = grid.Total_NumOrbs[atom]
        neighbor_orbitals = grid.Total_NumOrbs[neighbor]
        hks_block = _HamiltonianBlock(hks_channels, offset, neighbor_orbitals)

        _Calc_Hamiltonian_Blocked!(
            hks_block, N, atom_orbitals, neighbor_orbitals,
            ucell.MPI_NumOLG[pair], ucell.GridListAtom[atom],
            ucell.MPI_GListTAtoms1[pair], ucell.MPI_GListTAtoms2[pair],
            atom_matrix(orbitals_grid, atom),
            atom_matrix(orbitals_grid, neighbor), potential_grid,
            ucell.hamiltonian_orbital_scratch,
            ucell.hamiltonian_product_scratch,
        )
        offset += atom_orbitals*neighbor_orbitals
    end

    foreach(channel -> channel .*= grid.GridVol, hks_channels)
    return nothing
end


"""Blocked BLAS contraction for `<phi_neighbor|V_channel|phi_atom>`."""
function _Calc_Hamiltonian_Blocked!(hks_block, channel_count::Integer,
                                    atom_orbital_count::Integer,
                                    neighbor_orbital_count::Integer,
                                    overlap_count::Integer, atom_grid_points,
                                    atom_overlap_points, neighbor_overlap_points,
                                    atom_orbitals::Matrix{Float64},
                                    neighbor_orbitals::Matrix{Float64},
                                    potential_grid,
                                    orbital_scratch,
                                    product_scratch::Matrix{Float64})
    atom_block, neighbor_block, weighted_atom_block = orbital_scratch
    block_size = size(atom_block, 2)

    for block_start = 1:block_size:overlap_count
        points = min(block_size, overlap_count - block_start + 1)
        @inbounds for point = 1:points
            overlap_index = block_start + point - 1
            atom_grid_index = atom_overlap_points[overlap_index] + 1
            neighbor_grid_index = neighbor_overlap_points[overlap_index] + 1
            for orbital = 1:atom_orbital_count
                atom_block[orbital, point] =
                    atom_orbitals[orbital, atom_grid_index]
            end
            for orbital = 1:neighbor_orbital_count
                neighbor_block[orbital, point] =
                    neighbor_orbitals[orbital, neighbor_grid_index]
            end
        end

        for channel = 1:channel_count
            potential = potential_grid[channel]
            @inbounds for point = 1:points
                overlap_index = block_start + point - 1
                atom_grid_index = atom_overlap_points[overlap_index] + 1
                value = potential[atom_grid_points[atom_grid_index] + 1]
                @simd for orbital = 1:atom_orbital_count
                    weighted_atom_block[orbital, point] =
                        atom_block[orbital, point]*value
                end
            end

            @views mul!(
                product_scratch[1:neighbor_orbital_count, 1:atom_orbital_count],
                neighbor_block[1:neighbor_orbital_count, 1:points],
                transpose(weighted_atom_block[1:atom_orbital_count, 1:points]),
            )
            @inbounds for atom_orbital = 1:atom_orbital_count,
                          neighbor_orbital = 1:neighbor_orbital_count
                hks_block[neighbor_orbital, atom_orbital, channel] +=
                    product_scratch[neighbor_orbital, atom_orbital]
            end
        end
    end
    return nothing
end
