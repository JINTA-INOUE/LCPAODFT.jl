struct Hamiltonian
	SpinPol::String
	OLP::Vector{Float64}
	Hkin::Vector{Vector{Vector{Vector{Float64}}}}
	HNL::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
	iHNL::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
	HVNA::Vector{Vector{Vector{Vector{Float64}}}}
	NLPforce::Vector{Vector{Vector{Vector{Matrix{Float64}}}}}
	DS_VNAforce::Vector{Vector{Vector{Matrix{Float64}}}}
	HVNA2force::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
	HVNA3force::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
end


@timeit timer "Hamiltonian" function Hamiltonian(SpinPol::AbstractString, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)	
	
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
	Total_Hsize = system_grid.Total_Hsize

	
    # Overlap/Kinetic Matrix
	OLP = zeros(Float64, Total_Hsize)
    Hkin = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
	for atom = 1:Natom
		Hkin[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
		for Rn = 1:FNAN[atom]+1
			Hkin[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
			for ist = 1:Total_NumOrbs[atom]
				Hkin[atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
			end
		end
	end

	myrank == 0 && println("<Set_OLP_Kin>  Calculation of the overlap matrix")
	Set_OLP_Kin!(OLP, Hkin, pao, system_grid)



	Nspecies = length(pao)
	maxFNAN = maximum(FNAN)
	maxTotal_NumOrbs = maximum(Total_NumOrbs)
	maxL = maximum([pao[spe].Spe_MaxL_Basis for spe = 1:Nspecies]) + BufferL_ProVNA
    VNATotal_Num = (maxL+1)^2 * maxM

	DS_VNAforce = Vector{Vector{Vector{Matrix{Float64}}}}(undef, 4)
	for xyz = 1:4
		DS_VNAforce[xyz] = Vector{Vector{Matrix{Float64}}}(undef, Natom+1)
		for atom = 1:Natom+1
			if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
			DS_VNAforce[xyz][atom] = Vector{Matrix{Float64}}(undef, fan)
			for Rn = 1:fan
				DS_VNAforce[xyz][atom][Rn] = zeros(Float64, NO0, VNATotal_Num)
			end
		end
	end

	HVNA = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
	for atom = 1:Natom
		HVNA[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
		for Rn = 1:FNAN[atom]+1
			HVNA[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
			for ist = 1:Total_NumOrbs[atom]
				HVNA[atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
			end
		end
	end

	HVNA2force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
    HVNA3force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
	for xyz = 1:3
		HVNA2force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		HVNA3force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			HVNA2force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			HVNA3force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
				HVNA2force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, NO0)
				for ist = 1:NO0
					HVNA2force[xyz][atom][Rn][ist] = zeros(Float64, NO0)
				end

				NO1 = Total_NumOrbs[natn[atom][Rn]]
                HVNA3force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, NO1)
                for ist = 1:NO1
					HVNA3force[xyz][atom][Rn][ist] = zeros(Float64, NO1)
				end
			end
		end
	end

	myrank == 0 && println("<Set_ProExpn_VNA>  Calculation of the VNA projector matrix")
	Set_ProExpn_VNA!(DS_VNAforce, HVNA, HVNA2force, HVNA3force, pao, pspot, system_grid)





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

	maxVPS_j_Num = maximum(VPS_j_Num)
	maxNLTotal_Num = maximum(NLTotal_Num)
    NLPforce = Vector{Vector{Vector{Vector{Matrix{Float64}}}}}(undef, 4)
    for xyz = 1:4
        NLPforce[xyz] = Vector{Vector{Vector{Matrix{Float64}}}}(undef, Natom+1)
        for atom = 1:Natom+1
            if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
            
            NLPforce[xyz][atom] = Vector{Vector{Matrix{Float64}}}(undef, fan)
            for Rn = 1:fan
                if atom == Natom+1
                    VPS_j_dependency = maxVPS_j_Num
                    NO1 = maxNLTotal_Num
                else
                    jatom = natn[atom][Rn]
                    VPS_j_dependency = VPS_j_Num[jatom]
                    NO1 = NLTotal_Num[jatom]
                end
                
                NLPforce[xyz][atom][Rn] = Vector{Matrix{Float64}}(undef, VPS_j_dependency+1)
                for so = 1:VPS_j_dependency+1
                    NLPforce[xyz][atom][Rn][so] = zeros(Float64, NO0, NO1)
                end
            end
        end
    end


	if SpinPol ∈ ("off", "on")
		HNL = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 1)
		for spin = 1:1
			HNL[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
			for atom = 1:Natom
				HNL[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
				for Rn = 1:FNAN[atom]+1
					HNL[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
					for ist = 1:Total_NumOrbs[atom]
						HNL[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
					end
				end
			end
		end
		iHNL = [[[[[1.0]]]]]
	elseif SpinPol == "nc"
		HNL = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
		iHNL = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
		for spin = 1:3
			HNL[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
			iHNL[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
			for atom = 1:Natom
				HNL[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
				iHNL[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
				for Rn = 1:FNAN[atom]+1
					HNL[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
					iHNL[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
					for ist = 1:Total_NumOrbs[atom]
						HNL[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
						iHNL[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
					end
				end
			end
		end
	end
	    
	myrank == 0 && println("<Set_Nonlocal>  Calculation of the nonlocal matrix")
	Set_Nonlocal!(SpinPol, NLPforce, HNL, iHNL, pao, pspot, system_grid)
	

	return Hamiltonian(SpinPol, OLP, Hkin, HNL, iHNL, HVNA, NLPforce, DS_VNAforce, HVNA2force, HVNA3force)
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

	system_grid = ucell.system_grid
	Natom = system_grid.Natom
	FNAN = system_grid.FNAN
	natn = system_grid.natn
	Total_NumOrbs = system_grid.Total_NumOrbs
	SpinPol = Ham.SpinPol
	Hkin = Ham.Hkin
	HVNA = Ham.HVNA
	HNL = Ham.HNL


	if SpinPol == "off"
		Calc_MatrixElements_dVH_Vxc_off!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1], MPI.SUM, comm)
	elseif SpinPol == "on"
		Calc_MatrixElements_dVH_Vxc_on!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1], MPI.SUM, comm)
		MPI.Allreduce!(Hks[2], MPI.SUM, comm)
	elseif SpinPol == "nc"
		Calc_MatrixElements_dVH_Vxc_nc!(ucell, Orbs_Grid, Vpot_Grid, Hks)
		MPI.Allreduce!(Hks[1], MPI.SUM, comm)
		MPI.Allreduce!(Hks[2], MPI.SUM, comm)
		MPI.Allreduce!(Hks[3], MPI.SUM, comm)
		MPI.Allreduce!(Hks[4], MPI.SUM, comm)
	end
	
	

	if SpinPol == "off"
		Add_Hkin_HVNA_HNL!(Hks[1], Hkin, HVNA, HNL[1], Natom, FNAN, natn, Total_NumOrbs)
	elseif SpinPol == "on"
		Add_Hkin_HVNA_HNL!(Hks[1], Hkin, HVNA, HNL[1], Natom, FNAN, natn, Total_NumOrbs)
		Add_Hkin_HVNA_HNL!(Hks[2], Hkin, HVNA, HNL[1], Natom, FNAN, natn, Total_NumOrbs)
	elseif SpinPol == "nc"
		Add_Hkin_HVNA_HNL!(Hks[1], Hkin, HVNA, HNL[1], Natom, FNAN, natn, Total_NumOrbs)
		Add_Hkin_HVNA_HNL!(Hks[2], Hkin, HVNA, HNL[2], Natom, FNAN, natn, Total_NumOrbs)
		Add_HNL3!(Hks[3], HNL[3], Natom, FNAN, natn, Total_NumOrbs)
	end
end


@inline function Add_HNL3!(Hks, HNL, Natom, FNAN, natn, Total_NumOrbs)
	hst = 0
	@inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
		hst += 1
		Hks[hst] += HNL[atom][Rn][ist][jst]
	end
end


@inline function Add_Hkin_HVNA_HNL!(Hks, Hkin, HVNA, HNL, Natom, FNAN, natn, Total_NumOrbs)
	hst = 0
	@inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
		hst += 1
		Hks[hst] += Hkin[atom][Rn][ist][jst] + HVNA[atom][Rn][ist][jst] + HNL[atom][Rn][ist][jst]
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
		_Calc_Ham2_off!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid[1])

		@inbounds for ist = 1:NO0, jst = 1:NO1
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
	HksNum = MPHks[myrank+1]

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

	Hks_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 2)


	fill!(Hks[1], 0.0)
	fill!(Hks[2], 0.0)


	hst = 0
	for loop = 1:MPI_size

		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		fill!(Hks_temp, 0.0)
		_Calc_Ham2_on!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)

		@inbounds for ist = 1:NO0, jst = 1:NO1
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
	HksNum = MPHks[myrank+1]

	GridVol = ucell.system_grid.GridVol
	GridListAtom = ucell.GridListAtom
	MPI_NumOLG = ucell.MPI_NumOLG
	MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
	MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

	Hks_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 4)


	fill!(Hks[1], 0.0)
	fill!(Hks[2], 0.0)	
	fill!(Hks[3], 0.0)	
	fill!(Hks[4], 0.0)	


	hst = 0
	for loop = 1:MPI_size

		atom = MPI_atom[loop]
		jatom = MPI_natn[loop]
		NO0 = Total_NumOrbs[atom]
		NO1 = Total_NumOrbs[jatom]

		fill!(Hks_temp, 0.0)
		_Calc_Ham2_nc!(Hks_temp, NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], Orbs_Grid[atom], Orbs_Grid[jatom], Vpot_Grid)

		@inbounds for ist = 1:NO0, jst = 1:NO1
			hst += 1
			Hks[1][HksNum+hst] = GridVol*Hks_temp[jst,ist,1]
			Hks[2][HksNum+hst] = GridVol*Hks_temp[jst,ist,2]
			Hks[3][HksNum+hst] = GridVol*Hks_temp[jst,ist,3]
			Hks[4][HksNum+hst] = GridVol*Hks_temp[jst,ist,4]
		end
    end
end


function _Calc_Ham2_off!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	for Nog = 1:2:NumOLG-1

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		
		temp0 = Vpot_Grid[MN0]
		temp1 = Vpot_Grid[MN1]

		for ist = 1:NO0
			Sum0 = temp0 * Orbs_Grid1[Nc0][ist]
			Sum1 = temp1 * Orbs_Grid1[Nc1][ist]
			@inbounds for jst = 1:NO1
				Hks_temp[jst,ist] += Sum0 * Orbs_Grid2[Nh0][jst]
				Hks_temp[jst,ist] += Sum1 * Orbs_Grid2[Nh1][jst]
			end
		end
	end


	Nog1 = 2*div(NumOLG, 2)
	rem_NumOLG = rem(NumOLG, 2)
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


function _Calc_Ham2_on!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	for Nog = 1:2:NumOLG-1

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		
		temp0_up = Vpot_Grid[1][MN0]
		temp1_up = Vpot_Grid[1][MN1]

		temp0_dn = Vpot_Grid[2][MN0]
		temp1_dn = Vpot_Grid[2][MN1]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[Nc0][ist]
			orbs1_1 = Orbs_Grid1[Nc1][ist]

			Sum0_up = temp0_up * orbs1_0
			Sum1_up = temp1_up * orbs1_1

			Sum0_dn = temp0_dn * orbs1_0
			Sum1_dn = temp1_dn * orbs1_1

			@inbounds for jst = 1:NO1
				orbs2_0 = Orbs_Grid2[Nh0][jst]
				orbs2_1 = Orbs_Grid2[Nh1][jst]

				tmp = 0.0
				tmp += Sum0_up * orbs2_0
				tmp += Sum1_up * orbs2_1
				Hks_temp[jst,ist,1] += tmp

				tmp = 0.0
				tmp += Sum0_dn * orbs2_0
				tmp += Sum1_dn * orbs2_1
				Hks_temp[jst,ist,2] += tmp
			end
		end
	end


	Nog1 = 2*div(NumOLG, 2)
	rem_NumOLG = rem(NumOLG, 2)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_up = Vpot_Grid[1][MN]
        temp_dn = Vpot_Grid[2][MN]
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


function _Calc_Ham2_nc!(Hks_temp, NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, Orbs_Grid1, Orbs_Grid2, Vpot_Grid)

	for Nog = 1:2:NumOLG-1

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1

		temp0_uu = Vpot_Grid[1][MN0]
		temp1_uu = Vpot_Grid[1][MN1]

		temp0_dd = Vpot_Grid[2][MN0]
		temp1_dd = Vpot_Grid[2][MN1]

		temp0_ud_r = Vpot_Grid[3][MN0]
		temp1_ud_r = Vpot_Grid[3][MN1]

		temp0_ud_i = Vpot_Grid[4][MN0]
		temp1_ud_i = Vpot_Grid[4][MN1]

		for ist = 1:NO0

			orbs1_0 = Orbs_Grid1[Nc0][ist]
			orbs1_1 = Orbs_Grid1[Nc1][ist]

			Sum0_uu = temp0_uu * orbs1_0
			Sum1_uu = temp1_uu * orbs1_1

			Sum0_dd = temp0_dd * orbs1_0
			Sum1_dd = temp1_dd * orbs1_1

			Sum0_ud_r = temp0_ud_r * orbs1_0
			Sum1_ud_r = temp1_ud_r * orbs1_1

			Sum0_ud_i = temp0_ud_i * orbs1_0
			Sum1_ud_i = temp1_ud_i * orbs1_1

			@inbounds for jst = 1:NO1

				orbs2_0 = Orbs_Grid2[Nh0][jst]
				orbs2_1 = Orbs_Grid2[Nh1][jst]

				tmp = 0.0
				tmp += Sum0_uu * orbs2_0
				tmp += Sum1_uu * orbs2_1
				Hks_temp[jst,ist,1] += tmp

				tmp = 0.0
				tmp += Sum0_dd * orbs2_0
				tmp += Sum1_dd * orbs2_1
				Hks_temp[jst,ist,2] += tmp

				tmp = 0.0
				tmp += Sum0_ud_r * orbs2_0
				tmp += Sum1_ud_r * orbs2_1
				Hks_temp[jst,ist,3] += tmp

				tmp = 0.0
				tmp += Sum0_ud_i * orbs2_0
				tmp += Sum1_ud_i * orbs2_1
				Hks_temp[jst,ist,4] += tmp
			end
		end
	end


	Nog1 = 2*div(NumOLG, 2)
	rem_NumOLG = rem(NumOLG, 2)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp_uu = Vpot_Grid[1][MN]
        temp_dd = Vpot_Grid[2][MN]
        temp_ud_r = Vpot_Grid[3][MN]
        temp_ud_i = Vpot_Grid[4][MN]
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
