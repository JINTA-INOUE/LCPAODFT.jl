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
	FNAN = system_grid.FNAN
	natn = system_grid.natn
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

	myDS_VNAforce = 0
    for loop = 1:MPI_size, _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:VNATotal_Num
        myDS_VNAforce += 1
    end

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

	myNLPforce = 0
    for loop = 1:MPI_size, _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:NLTotal_Num[MPI_natn[loop]]
        myNLPforce += 1
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


@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_off!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, MPI_Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPI_size = ucell.system_grid.MPI_size

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
	orbital_scratch = ucell.hamiltonian_orbital_scratch
	product_scratch = ucell.hamiltonian_product_scratch
	fill!(MPI_Hks[1], 0.0)

	hst = 0
	for loop = 1:MPI_size
		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		block_size = NO0*NO1
		Hks_block = _HamiltonianBlock((MPI_Hks[1],), hst, NO1)
		if Orbs_Grid isa PackedOrbitalsGrid
			_Calc_Hamiltonian_Blocked!(Hks_block, 1, NO0, NO1,
				MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop],
				MPI_GListTAtoms2[loop], Orbs_Grid.data[atom],
				Orbs_Grid.data[jatom], Vpot_Grid, orbital_scratch,
				product_scratch)
		else
			_Calc_Ham8_off!(Hks_block, NO0, NO1, MPI_NumOLG[loop],
				GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop],
				Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid[1])
		end
		hst += block_size
    end
	@. MPI_Hks[1] *= GridVol
end


@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_on!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, MPI_Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPI_size = ucell.system_grid.MPI_size
	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
	orbital_scratch = ucell.hamiltonian_orbital_scratch
	product_scratch = ucell.hamiltonian_product_scratch

	fill!(MPI_Hks[1], 0.0)
	fill!(MPI_Hks[2], 0.0)

	hst = 0
	for loop = 1:MPI_size
		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		block_size = NO0*NO1
		Hks_block = _HamiltonianBlock((MPI_Hks[1], MPI_Hks[2]), hst, NO1)
		if Orbs_Grid isa PackedOrbitalsGrid
			_Calc_Hamiltonian_Blocked!(Hks_block, 2, NO0, NO1,
				MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop],
				MPI_GListTAtoms2[loop], Orbs_Grid.data[atom],
				Orbs_Grid.data[jatom], Vpot_Grid, orbital_scratch,
				product_scratch)
		else
			_Calc_Ham8_on!(Hks_block, NO0, NO1, MPI_NumOLG[loop],
				GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop],
				Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)
		end
		hst += block_size
    end
	@. MPI_Hks[1] *= GridVol
	@. MPI_Hks[2] *= GridVol
end


@timeit timer "Set_Hamiltonian" function Calc_MatrixElements_dVH_Vxc_nc!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, MPI_Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPI_size = ucell.system_grid.MPI_size
	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2
	orbital_scratch = ucell.hamiltonian_orbital_scratch
	product_scratch = ucell.hamiltonian_product_scratch

	fill!(MPI_Hks[1], 0.0)
	fill!(MPI_Hks[2], 0.0)	
	fill!(MPI_Hks[3], 0.0)	
	fill!(MPI_Hks[4], 0.0)	


	MPI_Hks1 = MPI_Hks[1]
	MPI_Hks2 = MPI_Hks[2]
	MPI_Hks3 = MPI_Hks[3]
	MPI_Hks4 = MPI_Hks[4]
	hst = 0
	for loop = 1:MPI_size
		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		block_size = NO0*NO1
		Hks_block = _HamiltonianBlock((MPI_Hks1, MPI_Hks2, MPI_Hks3, MPI_Hks4), hst, NO1)
		if Orbs_Grid isa PackedOrbitalsGrid
			_Calc_Hamiltonian_Blocked!(Hks_block, 4, NO0, NO1,
				MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop],
				MPI_GListTAtoms2[loop], Orbs_Grid.data[atom],
				Orbs_Grid.data[jatom], Vpot_Grid, orbital_scratch,
				product_scratch)
		else
			_Calc_Ham8_nc!(Hks_block, NO0, NO1, MPI_NumOLG[loop],
				GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop],
				Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)
		end
		hst += block_size
    end
	@. MPI_Hks1 *= GridVol
	@. MPI_Hks2 *= GridVol
	@. MPI_Hks3 *= GridVol
	@. MPI_Hks4 *= GridVol
end


"""Blocked BLAS contraction for `<phi_j|V_channel|phi_i>` on one atom pair."""
function _Calc_Hamiltonian_Blocked!(Hks_block, nchannels::Integer, NO0::Integer,
                                    NO1::Integer, NumOLG::Integer, GridListAtom,
                                    GListTAtoms1, GListTAtoms2,
                                    orbitals1::Matrix{Float64},
                                    orbitals2::Matrix{Float64}, Vpot_Grid,
                                    orbital_scratch,
                                    product_scratch::Matrix{Float64})
    orbital_block1, orbital_block2, weighted_block = orbital_scratch
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

        for channel = 1:nchannels
            potential = Vpot_Grid[channel]
            @inbounds for point = 1:points
                overlap_index = block_start + point - 1
                grid1 = GListTAtoms1[overlap_index] + 1
                value = potential[GridListAtom[grid1] + 1]
                @simd for orbital = 1:NO0
                    weighted_block[orbital, point] =
                        orbital_block1[orbital, point]*value
                end
            end

            @views mul!(product_scratch[1:NO1, 1:NO0],
                        orbital_block2[1:NO1, 1:points],
                        transpose(weighted_block[1:NO0, 1:points]))
            @inbounds for ist = 1:NO0, jst = 1:NO1
                Hks_block[jst, ist, channel] += product_scratch[jst, ist]
            end
        end
    end
    return nothing
end


function _Calc_Ham8_off!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	for Nog = 1:8:NumOLG-7

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		Nc2 = GListTAtoms1[Nog+2]+1
		Nc3 = GListTAtoms1[Nog+3]+1
		Nc4 = GListTAtoms1[Nog+4]+1
		Nc5 = GListTAtoms1[Nog+5]+1
		Nc6 = GListTAtoms1[Nog+6]+1
		Nc7 = GListTAtoms1[Nog+7]+1

		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		MN2 = GridListAtom[Nc2]+1
		MN3 = GridListAtom[Nc3]+1
		MN4 = GridListAtom[Nc4]+1
		MN5 = GridListAtom[Nc5]+1
		MN6 = GridListAtom[Nc6]+1
		MN7 = GridListAtom[Nc7]+1

		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		Nh2 = GListTAtoms2[Nog+2]+1
		Nh3 = GListTAtoms2[Nog+3]+1
		Nh4 = GListTAtoms2[Nog+4]+1
		Nh5 = GListTAtoms2[Nog+5]+1
		Nh6 = GListTAtoms2[Nog+6]+1
		Nh7 = GListTAtoms2[Nog+7]+1
		
		temp0 = Vpot_Grid[MN0]
		temp1 = Vpot_Grid[MN1]
		temp2 = Vpot_Grid[MN2]
		temp3 = Vpot_Grid[MN3]
		temp4 = Vpot_Grid[MN4]
		temp5 = Vpot_Grid[MN5]
		temp6 = Vpot_Grid[MN6]
		temp7 = Vpot_Grid[MN7]

		for ist = 1:NO0
			Sum0 = temp0 * Orbs_Grid1[Nc0][ist]
			Sum1 = temp1 * Orbs_Grid1[Nc1][ist]
			Sum2 = temp2 * Orbs_Grid1[Nc2][ist]
			Sum3 = temp3 * Orbs_Grid1[Nc3][ist]
			Sum4 = temp4 * Orbs_Grid1[Nc4][ist]
			Sum5 = temp5 * Orbs_Grid1[Nc5][ist]
			Sum6 = temp6 * Orbs_Grid1[Nc6][ist]
			Sum7 = temp7 * Orbs_Grid1[Nc7][ist]
			@inbounds for jst = 1:NO1
				Hks_temp[jst,ist] += Sum0 * Orbs_Grid2[Nh0][jst]
				Hks_temp[jst,ist] += Sum1 * Orbs_Grid2[Nh1][jst]
				Hks_temp[jst,ist] += Sum2 * Orbs_Grid2[Nh2][jst]
				Hks_temp[jst,ist] += Sum3 * Orbs_Grid2[Nh3][jst]
				Hks_temp[jst,ist] += Sum4 * Orbs_Grid2[Nh4][jst]
				Hks_temp[jst,ist] += Sum5 * Orbs_Grid2[Nh5][jst]
				Hks_temp[jst,ist] += Sum6 * Orbs_Grid2[Nh6][jst]
				Hks_temp[jst,ist] += Sum7 * Orbs_Grid2[Nh7][jst]
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp = Vpot_Grid[MN]
        for ist = 1:NO0
            Sum = temp * Orbs_Grid1[Nc][ist]
            @inbounds for jst = 1:NO1
                Hks_temp[jst,ist] += Sum * Orbs_Grid2[Nh][jst]
            end
        end
    end
end


function _Calc_Ham8_on!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	Vpot_Grid1 = Vpot_Grid[1]
	Vpot_Grid2 = Vpot_Grid[2]

	for Nog = 1:8:NumOLG-7

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		Nc2 = GListTAtoms1[Nog+2]+1
		Nc3 = GListTAtoms1[Nog+3]+1
		Nc4 = GListTAtoms1[Nog+4]+1
		Nc5 = GListTAtoms1[Nog+5]+1
		Nc6 = GListTAtoms1[Nog+6]+1
		Nc7 = GListTAtoms1[Nog+7]+1

		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		MN2 = GridListAtom[Nc2]+1
		MN3 = GridListAtom[Nc3]+1
		MN4 = GridListAtom[Nc4]+1
		MN5 = GridListAtom[Nc5]+1
		MN6 = GridListAtom[Nc6]+1
		MN7 = GridListAtom[Nc7]+1

		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		Nh2 = GListTAtoms2[Nog+2]+1
		Nh3 = GListTAtoms2[Nog+3]+1
		Nh4 = GListTAtoms2[Nog+4]+1
		Nh5 = GListTAtoms2[Nog+5]+1
		Nh6 = GListTAtoms2[Nog+6]+1
		Nh7 = GListTAtoms2[Nog+7]+1
		
		temp0_up = Vpot_Grid1[MN0]
		temp1_up = Vpot_Grid1[MN1]
		temp2_up = Vpot_Grid1[MN2]
		temp3_up = Vpot_Grid1[MN3]
		temp4_up = Vpot_Grid1[MN4]
		temp5_up = Vpot_Grid1[MN5]
		temp6_up = Vpot_Grid1[MN6]
		temp7_up = Vpot_Grid1[MN7]

		temp0_dn = Vpot_Grid2[MN0]
		temp1_dn = Vpot_Grid2[MN1]
		temp2_dn = Vpot_Grid2[MN2]
		temp3_dn = Vpot_Grid2[MN3]
		temp4_dn = Vpot_Grid2[MN4]
		temp5_dn = Vpot_Grid2[MN5]
		temp6_dn = Vpot_Grid2[MN6]
		temp7_dn = Vpot_Grid2[MN7]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[Nc0][ist]
			orbs1_1 = Orbs_Grid1[Nc1][ist]
			orbs1_2 = Orbs_Grid1[Nc2][ist]
			orbs1_3 = Orbs_Grid1[Nc3][ist]
			orbs1_4 = Orbs_Grid1[Nc4][ist]
			orbs1_5 = Orbs_Grid1[Nc5][ist]
			orbs1_6 = Orbs_Grid1[Nc6][ist]
			orbs1_7 = Orbs_Grid1[Nc7][ist]

			Sum0_up = temp0_up * orbs1_0
			Sum1_up = temp1_up * orbs1_1
			Sum2_up = temp2_up * orbs1_2
			Sum3_up = temp3_up * orbs1_3
			Sum4_up = temp4_up * orbs1_4
			Sum5_up = temp5_up * orbs1_5
			Sum6_up = temp6_up * orbs1_6
			Sum7_up = temp7_up * orbs1_7

			Sum0_dn = temp0_dn * orbs1_0
			Sum1_dn = temp1_dn * orbs1_1
			Sum2_dn = temp2_dn * orbs1_2
			Sum3_dn = temp3_dn * orbs1_3
			Sum4_dn = temp4_dn * orbs1_4
			Sum5_dn = temp5_dn * orbs1_5
			Sum6_dn = temp6_dn * orbs1_6
			Sum7_dn = temp7_dn * orbs1_7

			@inbounds for jst = 1:NO1
				orbs2_0 = Orbs_Grid2[Nh0][jst]
				orbs2_1 = Orbs_Grid2[Nh1][jst]
				orbs2_2 = Orbs_Grid2[Nh2][jst]
				orbs2_3 = Orbs_Grid2[Nh3][jst]
				orbs2_4 = Orbs_Grid2[Nh4][jst]
				orbs2_5 = Orbs_Grid2[Nh5][jst]
				orbs2_6 = Orbs_Grid2[Nh6][jst]
				orbs2_7 = Orbs_Grid2[Nh7][jst]

				tmp = 0.0
				tmp += Sum0_up * orbs2_0
				tmp += Sum1_up * orbs2_1
				tmp += Sum2_up * orbs2_2
				tmp += Sum3_up * orbs2_3
				tmp += Sum4_up * orbs2_4
				tmp += Sum5_up * orbs2_5
				tmp += Sum6_up * orbs2_6
				tmp += Sum7_up * orbs2_7
				Hks_temp[jst,ist,1] += tmp

				tmp = 0.0
				tmp += Sum0_dn * orbs2_0
				tmp += Sum1_dn * orbs2_1
				tmp += Sum2_dn * orbs2_2
				tmp += Sum3_dn * orbs2_3
				tmp += Sum4_dn * orbs2_4
				tmp += Sum5_dn * orbs2_5
				tmp += Sum6_dn * orbs2_6
				tmp += Sum7_dn * orbs2_7
				Hks_temp[jst,ist,2] += tmp
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_up = Vpot_Grid1[MN]
        temp_dn = Vpot_Grid2[MN]
        for ist = 1:NO0
			orbs1 = Orbs_Grid1[Nc][ist]
            Sum_up = temp_up * orbs1
            Sum_dn = temp_dn * orbs1
            @inbounds for jst = 1:NO1
				orbs2 = Orbs_Grid2[Nh][jst]
                Hks_temp[jst,ist,1] += Sum_up * orbs2
                Hks_temp[jst,ist,2] += Sum_dn * orbs2
            end
        end
    end
end


function _Calc_Ham8_nc!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	Vpot_Grid1 = Vpot_Grid[1]
	Vpot_Grid2 = Vpot_Grid[2]
	Vpot_Grid3 = Vpot_Grid[3]
	Vpot_Grid4 = Vpot_Grid[4]

	for Nog = 1:8:NumOLG-7

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		Nc2 = GListTAtoms1[Nog+2]+1
		Nc3 = GListTAtoms1[Nog+3]+1
		Nc4 = GListTAtoms1[Nog+4]+1
		Nc5 = GListTAtoms1[Nog+5]+1
		Nc6 = GListTAtoms1[Nog+6]+1
		Nc7 = GListTAtoms1[Nog+7]+1

		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		MN2 = GridListAtom[Nc2]+1
		MN3 = GridListAtom[Nc3]+1
		MN4 = GridListAtom[Nc4]+1
		MN5 = GridListAtom[Nc5]+1
		MN6 = GridListAtom[Nc6]+1
		MN7 = GridListAtom[Nc7]+1

		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		Nh2 = GListTAtoms2[Nog+2]+1
		Nh3 = GListTAtoms2[Nog+3]+1
		Nh4 = GListTAtoms2[Nog+4]+1
		Nh5 = GListTAtoms2[Nog+5]+1
		Nh6 = GListTAtoms2[Nog+6]+1
		Nh7 = GListTAtoms2[Nog+7]+1
		
		temp0_uu = Vpot_Grid1[MN0]
		temp1_uu = Vpot_Grid1[MN1]
		temp2_uu = Vpot_Grid1[MN2]
		temp3_uu = Vpot_Grid1[MN3]
		temp4_uu = Vpot_Grid1[MN4]
		temp5_uu = Vpot_Grid1[MN5]
		temp6_uu = Vpot_Grid1[MN6]
		temp7_uu = Vpot_Grid1[MN7]

		temp0_dd = Vpot_Grid2[MN0]
		temp1_dd = Vpot_Grid2[MN1]
		temp2_dd = Vpot_Grid2[MN2]
		temp3_dd = Vpot_Grid2[MN3]
		temp4_dd = Vpot_Grid2[MN4]
		temp5_dd = Vpot_Grid2[MN5]
		temp6_dd = Vpot_Grid2[MN6]
		temp7_dd = Vpot_Grid2[MN7]

		temp0_ud_r = Vpot_Grid3[MN0]
		temp1_ud_r = Vpot_Grid3[MN1]
		temp2_ud_r = Vpot_Grid3[MN2]
		temp3_ud_r = Vpot_Grid3[MN3]
		temp4_ud_r = Vpot_Grid3[MN4]
		temp5_ud_r = Vpot_Grid3[MN5]
		temp6_ud_r = Vpot_Grid3[MN6]
		temp7_ud_r = Vpot_Grid3[MN7]

		temp0_ud_i = Vpot_Grid4[MN0]
		temp1_ud_i = Vpot_Grid4[MN1]
		temp2_ud_i = Vpot_Grid4[MN2]
		temp3_ud_i = Vpot_Grid4[MN3]
		temp4_ud_i = Vpot_Grid4[MN4]
		temp5_ud_i = Vpot_Grid4[MN5]
		temp6_ud_i = Vpot_Grid4[MN6]
		temp7_ud_i = Vpot_Grid4[MN7]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[Nc0][ist]
			orbs1_1 = Orbs_Grid1[Nc1][ist]
			orbs1_2 = Orbs_Grid1[Nc2][ist]
			orbs1_3 = Orbs_Grid1[Nc3][ist]
			orbs1_4 = Orbs_Grid1[Nc4][ist]
			orbs1_5 = Orbs_Grid1[Nc5][ist]
			orbs1_6 = Orbs_Grid1[Nc6][ist]
			orbs1_7 = Orbs_Grid1[Nc7][ist]

			Sum0_uu = temp0_uu * orbs1_0
			Sum1_uu = temp1_uu * orbs1_1
			Sum2_uu = temp2_uu * orbs1_2
			Sum3_uu = temp3_uu * orbs1_3
			Sum4_uu = temp4_uu * orbs1_4
			Sum5_uu = temp5_uu * orbs1_5
			Sum6_uu = temp6_uu * orbs1_6
			Sum7_uu = temp7_uu * orbs1_7

			Sum0_dd = temp0_dd * orbs1_0
			Sum1_dd = temp1_dd * orbs1_1
			Sum2_dd = temp2_dd * orbs1_2
			Sum3_dd = temp3_dd * orbs1_3
			Sum4_dd = temp4_dd * orbs1_4
			Sum5_dd = temp5_dd * orbs1_5
			Sum6_dd = temp6_dd * orbs1_6
			Sum7_dd = temp7_dd * orbs1_7

			Sum0_ud_r = temp0_ud_r * orbs1_0
			Sum1_ud_r = temp1_ud_r * orbs1_1
			Sum2_ud_r = temp2_ud_r * orbs1_2
			Sum3_ud_r = temp3_ud_r * orbs1_3
			Sum4_ud_r = temp4_ud_r * orbs1_4
			Sum5_ud_r = temp5_ud_r * orbs1_5
			Sum6_ud_r = temp6_ud_r * orbs1_6
			Sum7_ud_r = temp7_ud_r * orbs1_7

			Sum0_ud_i = temp0_ud_i * orbs1_0
			Sum1_ud_i = temp1_ud_i * orbs1_1
			Sum2_ud_i = temp2_ud_i * orbs1_2
			Sum3_ud_i = temp3_ud_i * orbs1_3
			Sum4_ud_i = temp4_ud_i * orbs1_4
			Sum5_ud_i = temp5_ud_i * orbs1_5
			Sum6_ud_i = temp6_ud_i * orbs1_6
			Sum7_ud_i = temp7_ud_i * orbs1_7

			@inbounds for jst = 1:NO1

				orbs2_0 = Orbs_Grid2[Nh0][jst]
				orbs2_1 = Orbs_Grid2[Nh1][jst]
				orbs2_2 = Orbs_Grid2[Nh2][jst]
				orbs2_3 = Orbs_Grid2[Nh3][jst]
				orbs2_4 = Orbs_Grid2[Nh4][jst]
				orbs2_5 = Orbs_Grid2[Nh5][jst]
				orbs2_6 = Orbs_Grid2[Nh6][jst]
				orbs2_7 = Orbs_Grid2[Nh7][jst]

				tmp = 0.0
				tmp += Sum0_uu * orbs2_0
				tmp += Sum1_uu * orbs2_1
				tmp += Sum2_uu * orbs2_2
				tmp += Sum3_uu * orbs2_3
				tmp += Sum4_uu * orbs2_4
				tmp += Sum5_uu * orbs2_5
				tmp += Sum6_uu * orbs2_6
				tmp += Sum7_uu * orbs2_7
				Hks_temp[jst,ist,1] += tmp

				tmp = 0.0
				tmp += Sum0_dd * orbs2_0
				tmp += Sum1_dd * orbs2_1
				tmp += Sum2_dd * orbs2_2
				tmp += Sum3_dd * orbs2_3
				tmp += Sum4_dd * orbs2_4
				tmp += Sum5_dd * orbs2_5
				tmp += Sum6_dd * orbs2_6
				tmp += Sum7_dd * orbs2_7
				Hks_temp[jst,ist,2] += tmp

				tmp = 0.0
				tmp += Sum0_ud_r * orbs2_0
				tmp += Sum1_ud_r * orbs2_1
				tmp += Sum2_ud_r * orbs2_2
				tmp += Sum3_ud_r * orbs2_3
				tmp += Sum4_ud_r * orbs2_4
				tmp += Sum5_ud_r * orbs2_5
				tmp += Sum6_ud_r * orbs2_6
				tmp += Sum7_ud_r * orbs2_7
				Hks_temp[jst,ist,3] += tmp

				tmp = 0.0
				tmp += Sum0_ud_i * orbs2_0
				tmp += Sum1_ud_i * orbs2_1
				tmp += Sum2_ud_i * orbs2_2
				tmp += Sum3_ud_i * orbs2_3
				tmp += Sum4_ud_i * orbs2_4
				tmp += Sum5_ud_i * orbs2_5
				tmp += Sum6_ud_i * orbs2_6
				tmp += Sum7_ud_i * orbs2_7
				Hks_temp[jst,ist,4] += tmp
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_uu = Vpot_Grid1[MN]
        temp_dd = Vpot_Grid2[MN]
        temp_ud_r = Vpot_Grid3[MN]
        temp_ud_i = Vpot_Grid4[MN]
        for ist = 1:NO0
			orbs1 = Orbs_Grid1[Nc][ist]
            Sum_uu = temp_uu * orbs1
            Sum_dd = temp_dd * orbs1
            Sum_ud_r = temp_ud_r * orbs1
            Sum_ud_i = temp_ud_i * orbs1
            @inbounds for jst = 1:NO1
				orbs2 = Orbs_Grid2[Nh][jst]
                Hks_temp[jst,ist,1] += Sum_uu * orbs2
                Hks_temp[jst,ist,2] += Sum_dd * orbs2
                Hks_temp[jst,ist,3] += Sum_ud_r * orbs2
                Hks_temp[jst,ist,4] += Sum_ud_i * orbs2
            end
        end
    end
end
