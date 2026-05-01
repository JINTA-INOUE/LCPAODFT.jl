struct Hamiltonian
	SpinPol::String
	Hkin::Vector{Float64}
	HNL::Vector{Vector{Float64}}
	iHNL::Union{Vector{Vector{Float64}}, Nothing}
	HVNA::Vector{Float64}
    OLP::Vector{Float64}
end




"""
    Ham = Hamiltonian(...)

Create an instance of `Hamiltonian`.

Mandatory arguments:

- `atoms`: an instance of `Atoms`
- `LatVecs`: Lattice Vectors with unit
- `system`: System name (`Atom`, `Cluster`, `Crystal`)

The following is the most commonly used optional arguments:
- `Ecut`: energy cutoff for real space grids with unit (default `Ecut = 100.0(Ry)`)
"""
@timeit timer "Hamiltonian" function Hamiltonian(SpinPol::AbstractString, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)	
	
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	Total_Hsize = system_grid.Total_Hsize

	
    # Overlap/Kinetic Matrix
    OLP = zeros(Float64, Total_Hsize)
    Hkin = zeros(Float64, Total_Hsize)

	myrank == 0 && println("<Set_OLP_Kin>  Calculation of the overlap matrix")
	Set_OLP_Kin!(OLP, Hkin, pao, system_grid)


	# Neutral Potentials Matrix
    HVNA = zeros(Float64, Total_Hsize)

	
	myrank == 0 && println("<Set_ProExpn_VNA>  Calculation of the VNA projector matrix")
	Set_ProExpn_VNA!(HVNA, pao, pspot, system_grid)



	# Nonlocal Potentials Matrix
	if SpinPol ∈ ("off", "on")
		spinmax = 1
	elseif SpinPol == "nc"
		spinmax = 3
	else
		println("Now SpinPol is $SpinPol")
		error("please check SpinPol")
	end

	if SpinPol ∈ ("off", "on")
		HNL = Vector{Vector{Float64}}(undef, 1)
		HNL[1] = zeros(Float64, Total_Hsize)
		iHNL = nothing
	elseif SpinPol == "nc"
		HNL = Vector{Vector{Float64}}(undef, 3)
		iHNL = Vector{Vector{Float64}}(undef, 3)
		for spin = 1:3
			HNL[spin] = zeros(Float64, Total_Hsize)
			iHNL[spin] = zeros(Float64, Total_Hsize)
		end
	end
    
	myrank == 0 && println("<Set_Nonlocal>  Calculation of the nonlocal matrix")
	Set_Nonlocal!(SpinPol, HNL, iHNL, pao, pspot, system_grid)


	return Hamiltonian(SpinPol, Hkin, HNL, iHNL, HVNA, OLP)
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
@timeit timer "Set_Hamiltonian" function Set_Hamiltonian!(Ham::Hamiltonian, ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, Hks)
	
    comm = MPI.COMM_WORLD
	myrank = MPI.Comm_rank(comm)

	SpinPol = Ham.SpinPol
	Hkin = Ham.Hkin
	HVNA = Ham.HVNA
	HNL = Ham.HNL
	Total_Hsize = ucell.system_grid.Total_Hsize


	if SpinPol == "off"
		Calc_MatrixElements_dVH_Vxc_off!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1],MPI.SUM,comm)
	elseif SpinPol == "on"
		Calc_MatrixElements_dVH_Vxc_on!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1],MPI.SUM,comm)
		MPI.Allreduce!(Hks[2],MPI.SUM,comm)
	elseif SpinPol == "nc"
		Calc_MatrixElements_dVH_Vxc_nc!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1],MPI.SUM,comm)
		MPI.Allreduce!(Hks[2],MPI.SUM,comm)
		MPI.Allreduce!(Hks[3],MPI.SUM,comm)
		MPI.Allreduce!(Hks[4],MPI.SUM,comm)
	end
	
	

	if SpinPol == "off"
		for hst = 1:Total_Hsize
			Hks[1][hst] += Hkin[hst] + HVNA[hst] + HNL[1][hst]
		end
	elseif SpinPol == "on"
		for hst = 1:Total_Hsize
			tmp = Hkin[hst] + HVNA[hst] + HNL[1][hst]
			Hks[1][hst] += tmp
			Hks[2][hst] += tmp
		end
	elseif SpinPol == "nc"
		for hst = 1:Total_Hsize
			tmp = Hkin[hst] + HVNA[hst]
			Hks[1][hst] += tmp + HNL[1][hst]
			Hks[2][hst] += tmp + HNL[2][hst]
			Hks[3][hst] += HNL[3][hst]
		end
	end
end


function Calc_MatrixElements_dVH_Vxc_off!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPHks = ucell.system_grid.MPHks
	MPI_size = ucell.system_grid.MPI_size
	HksNum = MPHks[myrank+1]

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

	Hks_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs))

	

	fill!(Hks[1], 0.0)

	hst = 0
	for loop = 1:MPI_size

		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		fill!(Hks_temp, 0.0)
		_Calc_Ham8_off!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid[1])

		for ist = 1:NO0, jst = 1:NO1
			hst += 1
			Hks[1][HksNum+hst] = GridVol*Hks_temp[jst,ist]
		end
    end
end


function Calc_MatrixElements_dVH_Vxc_on!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPHks = ucell.system_grid.MPHks
	MPI_size = ucell.system_grid.MPI_size

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

	Hks_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 2)


	HksNum = MPHks[myrank+1]


	hst = 0
	for loop = 1:MPI_size

		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		fill!(Hks_temp, 0.0)
		_Calc_Ham8_on!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)

		for ist = 1:NO0, jst = 1:NO1
			hst += 1
			Hks[1][HksNum+hst] = GridVol*Hks_temp[jst,ist,1]
			Hks[2][HksNum+hst] = GridVol*Hks_temp[jst,ist,2]
		end
    end
end


function Calc_MatrixElements_dVH_Vxc_nc!(ucell::UCell, Orbs_Grid, Vpot_Grid::Vector{Vector{Float64}}, Hks)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

	MPI_atom = ucell.system_grid.MPI_atom
	MPI_natn = ucell.system_grid.MPI_natn
	Total_NumOrbs = ucell.system_grid.Total_NumOrbs
	MPHks = ucell.system_grid.MPHks
	MPI_size = ucell.system_grid.MPI_size

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

	Hks_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 4)


	HksNum = MPHks[myrank+1]


	hst = 0
	for loop = 1:MPI_size

		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		fill!(Hks_temp, 0.0)
		_Calc_Ham8_nc!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)

		for ist = 1:NO0, jst = 1:NO1
			hst += 1
			Hks[1][HksNum+hst] = GridVol*Hks_temp[jst,ist,1]
			Hks[2][HksNum+hst] = GridVol*Hks_temp[jst,ist,2]
			Hks[3][HksNum+hst] = GridVol*Hks_temp[jst,ist,3]
			Hks[4][HksNum+hst] = GridVol*Hks_temp[jst,ist,4]
		end
    end
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
			Sum0 = temp0 * Orbs_Grid1[ist][Nc0]
			Sum1 = temp1 * Orbs_Grid1[ist][Nc1]
			Sum2 = temp2 * Orbs_Grid1[ist][Nc2]
			Sum3 = temp3 * Orbs_Grid1[ist][Nc3]
			Sum4 = temp4 * Orbs_Grid1[ist][Nc4]
			Sum5 = temp5 * Orbs_Grid1[ist][Nc5]
			Sum6 = temp6 * Orbs_Grid1[ist][Nc6]
			Sum7 = temp7 * Orbs_Grid1[ist][Nc7]
			for jst = 1:NO1
				Hks_temp[jst,ist] += Sum0 * Orbs_Grid2[jst][Nh0]
				Hks_temp[jst,ist] += Sum1 * Orbs_Grid2[jst][Nh1]
				Hks_temp[jst,ist] += Sum2 * Orbs_Grid2[jst][Nh2]
				Hks_temp[jst,ist] += Sum3 * Orbs_Grid2[jst][Nh3]
				Hks_temp[jst,ist] += Sum4 * Orbs_Grid2[jst][Nh4]
				Hks_temp[jst,ist] += Sum5 * Orbs_Grid2[jst][Nh5]
				Hks_temp[jst,ist] += Sum6 * Orbs_Grid2[jst][Nh6]
				Hks_temp[jst,ist] += Sum7 * Orbs_Grid2[jst][Nh7]
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
            Sum = temp * Orbs_Grid1[ist][Nc]
            for jst = 1:NO1
                Hks_temp[jst,ist] += Sum * Orbs_Grid2[jst][Nh]
            end
        end
    end
end




function _Calc_Ham8_on!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

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
		
		temp0_up = Vpot_Grid[1][MN0]
		temp1_up = Vpot_Grid[1][MN1]
		temp2_up = Vpot_Grid[1][MN2]
		temp3_up = Vpot_Grid[1][MN3]
		temp4_up = Vpot_Grid[1][MN4]
		temp5_up = Vpot_Grid[1][MN5]
		temp6_up = Vpot_Grid[1][MN6]
		temp7_up = Vpot_Grid[1][MN7]

		temp0_dn = Vpot_Grid[2][MN0]
		temp1_dn = Vpot_Grid[2][MN1]
		temp2_dn = Vpot_Grid[2][MN2]
		temp3_dn = Vpot_Grid[2][MN3]
		temp4_dn = Vpot_Grid[2][MN4]
		temp5_dn = Vpot_Grid[2][MN5]
		temp6_dn = Vpot_Grid[2][MN6]
		temp7_dn = Vpot_Grid[2][MN7]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[ist][Nc0]
			orbs1_1 = Orbs_Grid1[ist][Nc1]
			orbs1_2 = Orbs_Grid1[ist][Nc2]
			orbs1_3 = Orbs_Grid1[ist][Nc3]
			orbs1_4 = Orbs_Grid1[ist][Nc4]
			orbs1_5 = Orbs_Grid1[ist][Nc5]
			orbs1_6 = Orbs_Grid1[ist][Nc6]
			orbs1_7 = Orbs_Grid1[ist][Nc7]

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

			for jst = 1:NO1
				orbs2_0 = Orbs_Grid2[jst][Nh0]
				orbs2_1 = Orbs_Grid2[jst][Nh1]
				orbs2_2 = Orbs_Grid2[jst][Nh2]
				orbs2_3 = Orbs_Grid2[jst][Nh3]
				orbs2_4 = Orbs_Grid2[jst][Nh4]
				orbs2_5 = Orbs_Grid2[jst][Nh5]
				orbs2_6 = Orbs_Grid2[jst][Nh6]
				orbs2_7 = Orbs_Grid2[jst][Nh7]

				Hks_temp[jst,ist,1] += Sum0_up * orbs2_0
				Hks_temp[jst,ist,1] += Sum1_up * orbs2_1
				Hks_temp[jst,ist,1] += Sum2_up * orbs2_2
				Hks_temp[jst,ist,1] += Sum3_up * orbs2_3
				Hks_temp[jst,ist,1] += Sum4_up * orbs2_4
				Hks_temp[jst,ist,1] += Sum5_up * orbs2_5
				Hks_temp[jst,ist,1] += Sum6_up * orbs2_6
				Hks_temp[jst,ist,1] += Sum7_up * orbs2_7

				Hks_temp[jst,ist,2] += Sum0_dn * orbs2_0
				Hks_temp[jst,ist,2] += Sum1_dn * orbs2_1
				Hks_temp[jst,ist,2] += Sum2_dn * orbs2_2
				Hks_temp[jst,ist,2] += Sum3_dn * orbs2_3
				Hks_temp[jst,ist,2] += Sum4_dn * orbs2_4
				Hks_temp[jst,ist,2] += Sum5_dn * orbs2_5
				Hks_temp[jst,ist,2] += Sum6_dn * orbs2_6
				Hks_temp[jst,ist,2] += Sum7_dn * orbs2_7
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_up = Vpot_Grid[1][MN]
        temp_dn = Vpot_Grid[2][MN]
        for ist = 1:NO0
			orbs1 = Orbs_Grid1[ist][Nc]
            Sum_up = temp_up * orbs1
            Sum_dn = temp_dn * orbs1
            for jst = 1:NO1
				orbs2 = Orbs_Grid2[jst][Nh]
                Hks_temp[jst,ist,1] += Sum_up * orbs2
                Hks_temp[jst,ist,2] += Sum_dn * orbs2
            end
        end
    end
end




function _Calc_Ham8_nc!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

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
		
		temp0_uu = Vpot_Grid[1][MN0]
		temp1_uu = Vpot_Grid[1][MN1]
		temp2_uu = Vpot_Grid[1][MN2]
		temp3_uu = Vpot_Grid[1][MN3]
		temp4_uu = Vpot_Grid[1][MN4]
		temp5_uu = Vpot_Grid[1][MN5]
		temp6_uu = Vpot_Grid[1][MN6]
		temp7_uu = Vpot_Grid[1][MN7]

		temp0_dd = Vpot_Grid[2][MN0]
		temp1_dd = Vpot_Grid[2][MN1]
		temp2_dd = Vpot_Grid[2][MN2]
		temp3_dd = Vpot_Grid[2][MN3]
		temp4_dd = Vpot_Grid[2][MN4]
		temp5_dd = Vpot_Grid[2][MN5]
		temp6_dd = Vpot_Grid[2][MN6]
		temp7_dd = Vpot_Grid[2][MN7]

		temp0_ud_r = Vpot_Grid[3][MN0]
		temp1_ud_r = Vpot_Grid[3][MN1]
		temp2_ud_r = Vpot_Grid[3][MN2]
		temp3_ud_r = Vpot_Grid[3][MN3]
		temp4_ud_r = Vpot_Grid[3][MN4]
		temp5_ud_r = Vpot_Grid[3][MN5]
		temp6_ud_r = Vpot_Grid[3][MN6]
		temp7_ud_r = Vpot_Grid[3][MN7]

		temp0_ud_i = Vpot_Grid[4][MN0]
		temp1_ud_i = Vpot_Grid[4][MN1]
		temp2_ud_i = Vpot_Grid[4][MN2]
		temp3_ud_i = Vpot_Grid[4][MN3]
		temp4_ud_i = Vpot_Grid[4][MN4]
		temp5_ud_i = Vpot_Grid[4][MN5]
		temp6_ud_i = Vpot_Grid[4][MN6]
		temp7_ud_i = Vpot_Grid[4][MN7]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[ist][Nc0]
			orbs1_1 = Orbs_Grid1[ist][Nc1]
			orbs1_2 = Orbs_Grid1[ist][Nc2]
			orbs1_3 = Orbs_Grid1[ist][Nc3]
			orbs1_4 = Orbs_Grid1[ist][Nc4]
			orbs1_5 = Orbs_Grid1[ist][Nc5]
			orbs1_6 = Orbs_Grid1[ist][Nc6]
			orbs1_7 = Orbs_Grid1[ist][Nc7]

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

			for jst = 1:NO1

				orbs2_0 = Orbs_Grid2[jst][Nh0]
				orbs2_1 = Orbs_Grid2[jst][Nh1]
				orbs2_2 = Orbs_Grid2[jst][Nh2]
				orbs2_3 = Orbs_Grid2[jst][Nh3]
				orbs2_4 = Orbs_Grid2[jst][Nh4]
				orbs2_5 = Orbs_Grid2[jst][Nh5]
				orbs2_6 = Orbs_Grid2[jst][Nh6]
				orbs2_7 = Orbs_Grid2[jst][Nh7]

				Hks_temp[jst,ist,1] += Sum0_uu * orbs2_0
				Hks_temp[jst,ist,1] += Sum1_uu * orbs2_1
				Hks_temp[jst,ist,1] += Sum2_uu * orbs2_2
				Hks_temp[jst,ist,1] += Sum3_uu * orbs2_3
				Hks_temp[jst,ist,1] += Sum4_uu * orbs2_4
				Hks_temp[jst,ist,1] += Sum5_uu * orbs2_5
				Hks_temp[jst,ist,1] += Sum6_uu * orbs2_6
				Hks_temp[jst,ist,1] += Sum7_uu * orbs2_7

				Hks_temp[jst,ist,2] += Sum0_dd * orbs2_0
				Hks_temp[jst,ist,2] += Sum1_dd * orbs2_1
				Hks_temp[jst,ist,2] += Sum2_dd * orbs2_2
				Hks_temp[jst,ist,2] += Sum3_dd * orbs2_3
				Hks_temp[jst,ist,2] += Sum4_dd * orbs2_4
				Hks_temp[jst,ist,2] += Sum5_dd * orbs2_5
				Hks_temp[jst,ist,2] += Sum6_dd * orbs2_6
				Hks_temp[jst,ist,2] += Sum7_dd * orbs2_7

				Hks_temp[jst,ist,3] += Sum0_ud_r * orbs2_0
				Hks_temp[jst,ist,3] += Sum1_ud_r * orbs2_1
				Hks_temp[jst,ist,3] += Sum2_ud_r * orbs2_2
				Hks_temp[jst,ist,3] += Sum3_ud_r * orbs2_3
				Hks_temp[jst,ist,3] += Sum4_ud_r * orbs2_4
				Hks_temp[jst,ist,3] += Sum5_ud_r * orbs2_5
				Hks_temp[jst,ist,3] += Sum6_ud_r * orbs2_6
				Hks_temp[jst,ist,3] += Sum7_ud_r * orbs2_7

				Hks_temp[jst,ist,4] += Sum0_ud_i * orbs2_0
				Hks_temp[jst,ist,4] += Sum1_ud_i * orbs2_1
				Hks_temp[jst,ist,4] += Sum2_ud_i * orbs2_2
				Hks_temp[jst,ist,4] += Sum3_ud_i * orbs2_3
				Hks_temp[jst,ist,4] += Sum4_ud_i * orbs2_4
				Hks_temp[jst,ist,4] += Sum5_ud_i * orbs2_5
				Hks_temp[jst,ist,4] += Sum6_ud_i * orbs2_6
				Hks_temp[jst,ist,4] += Sum7_ud_i * orbs2_7
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_uu = Vpot_Grid[1][MN]
        temp_dd = Vpot_Grid[2][MN]
        temp_ud_r = Vpot_Grid[3][MN]
        temp_ud_i = Vpot_Grid[4][MN]
        for ist = 1:NO0
			orbs1 = Orbs_Grid1[ist][Nc]
            Sum_uu = temp_uu * orbs1
            Sum_dd = temp_dd * orbs1
            Sum_ud_r = temp_ud_r * orbs1
            Sum_ud_i = temp_ud_i * orbs1
            for jst = 1:NO1
				orbs2 = Orbs_Grid2[jst][Nh]
                Hks_temp[jst,ist,1] += Sum_uu * orbs2
                Hks_temp[jst,ist,2] += Sum_dd * orbs2
                Hks_temp[jst,ist,3] += Sum_ud_r * orbs2
                Hks_temp[jst,ist,4] += Sum_ud_i * orbs2
            end
        end
    end
end