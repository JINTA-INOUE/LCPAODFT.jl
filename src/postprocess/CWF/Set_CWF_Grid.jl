@timeit timer "Set_CWF_Grid" function Set_CWF_Grid(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, CWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
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
    Ecut = cwf_setup.Ecut
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    Atoms_Cut1 = material.Atoms_Cut1
    Grid_Origin = material.Grid_Origin
    CWF_Plot_Cube = cwf_setup.CWF_Plot_Cube
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    filename = cwf_setup.filename

    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)
    pao = Vector{PAO}(undef, Nspecies)
    for spe = 1:Nspecies
        pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end
        
    ucell = UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs)
    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)


    if SpinPol ∈ ("off", "on")
        Set_CWF_Grid_Col(filename, material, ucell, CWF_Plot_Cube, CWF_Plot_SuperCells, CWF_ExpnCoef, Orbs_Grid)
    elseif SpinPol == "nc"
        Set_CWF_Grid_NonCol(filename, material, ucell, CWF_Plot_Cube, CWF_Plot_SuperCells, CWF_ExpnCoef, Orbs_Grid) 
    end
    MPI.Barrier(comm)
end


function Set_CWF_Grid_Col(filename::AbstractString, material::LCPAO_model, ucell::UCell, CWF_Plot_Cube, CWF_Plot_SuperCells, CWF_ExpnCoef, Orbs_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Natom = material.Natom
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    Atoms_symbol = material.Atoms_symbol
    Gxyz_scf = material.Gxyz
    
    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    Ngrid = system_grid.Ngrid
    Grid_Origin = system_grid.Grid_Origin
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_TNumGrid = prod(Ngrid)*Plot_NCell


    CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid = Set_CWF_UCell(CWF_Plot_SuperCells, ucell)
    MPI.Barrier(comm)

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    N_Polt_Cube = length(CWF_Plot_Cube)
    Nloop = spinsize*N_Polt_Cube
    myrange = split_evenly(1:Nloop, nprocs)
    Wannier_Orbs_Grid = zeros(Float64, CWF_TNumGrid)

    for loop in myrange[myrank+1]
        spin = div(loop-1, N_Polt_Cube)+1
        proj = CWF_Plot_Cube[mod(loop-1, N_Polt_Cube)+1]
        println("  Write myrank = $myrank   $(filename)_CWF$(spin)_$(proj).cube")
        
        fill!(Wannier_Orbs_Grid, 0.0)
        for cell = 1:Plot_NCell
            ExpnCoef = CWF_ExpnCoef[spin][proj][cell]
            for atom = 1:Natom
                cwf_proj = MP[atom]
                NO0 = Total_NumOrbs[atom]
                _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid.data[atom], ExpnCoef, Wannier_Orbs_Grid)
            end
        end

        open("$(filename)_CWF$(spin)_$(proj).cube", "w") do data
            Write_Wannier_CubeInfo(data, Atoms_symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Gxyz_scf, Latvecs, gLatvecs, Ngrid)
            Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, Wannier_Orbs_Grid)
        end
    end
end


function Set_CWF_Grid_NonCol(filename::AbstractString, material::LCPAO_model, ucell::UCell, CWF_Plot_Cube, CWF_Plot_SuperCells, CWF_ExpnCoef, Orbs_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

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
    
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_TNumGrid = prod(Ngrid)*Plot_NCell


    CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid = Set_CWF_UCell(CWF_Plot_SuperCells, ucell)

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    N_Polt_Cube = length(CWF_Plot_Cube)
    Nloop = N_Polt_Cube
    myrange = split_evenly(1:Nloop, nprocs)
    Wannier_Orbs_Grid = zeros(ComplexF64, CWF_TNumGrid)

    for loop in myrange[myrank+1]
        proj = CWF_Plot_Cube[loop]
        println("  Write myrank = $myrank   $(filename)_CWF$(proj).nccube")
        open("$(filename)_CWF$(proj).nccube", "w") do data
            Write_Wannier_CubeInfo(data, Atoms_symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Gxyz_scf, Latvecs, gLatvecs, Ngrid)
            for spin = 1:2
                fill!(Wannier_Orbs_Grid, 0.0)
                spin_site = ifelse(spin==1, 0, fsize)
                for cell = 1:Plot_NCell
                    ExpnCoef = CWF_ExpnCoef[proj][cell]
                    for atom = 1:Natom
                        cwf_proj = MP[atom]
                        NO0 = Total_NumOrbs[atom]
                        _Calc_CWF_Grid8!(spin_site+cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid.data[atom], ExpnCoef, Wannier_Orbs_Grid)
                    end
                end

                Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, Wannier_Orbs_Grid, real)
                Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, Wannier_Orbs_Grid, imag)
            end
        end
    end
end


function Set_CWF_UCell(CWF_Plot_SuperCells, ucell::UCell)

    CWF_Plot_SuperCells1, CWF_Plot_SuperCells2, CWF_Plot_SuperCells3 = CWF_Plot_SuperCells

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    atv_ijk = system_grid.atv_ijk
    Ngrid = system_grid.Ngrid
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_GridN_Atom = Vector{Vector{Int32}}(undef, Plot_NCell)
    CWF_GridOrbs_Grid = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    CWF_GridListAtom = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    output_n2 = Ngrid2*(2*CWF_Plot_SuperCells2+1)
    output_n3 = Ngrid3*(2*CWF_Plot_SuperCells3+1)

    cell = 0
    for l1 = -CWF_Plot_SuperCells1:CWF_Plot_SuperCells1,
        l2 = -CWF_Plot_SuperCells2:CWF_Plot_SuperCells2,
        l3 = -CWF_Plot_SuperCells3:CWF_Plot_SuperCells3
        cell += 1
        grid_counts = zeros(Int32, Natom)
        grid_orbitals = Vector{Vector{Int32}}(undef, Natom)
        grid_indices = Vector{Vector{Int32}}(undef, Natom)

        for atom = 1:Natom
            atom_grid_count = Int(GridN_Atom[atom])
            orbital_indices = Vector{Int32}(undef, atom_grid_count)
            output_indices = Vector{Int32}(undef, atom_grid_count)
            gridn = 0
            @inbounds for Nc = 1:atom_grid_count
                GNc = GridListAtom[atom][Nc]
                GRc = CellListAtom[atom][Nc]+1
                m1 = l1 + atv_ijk[GRc][1]
                m2 = l2 + atv_ijk[GRc][2]
                m3 = l3 + atv_ijk[GRc][3]
                if abs(m1)<=CWF_Plot_SuperCells1 && abs(m2)<=CWF_Plot_SuperCells2 && abs(m3)<=CWF_Plot_SuperCells3
                    n1 = div(GNc, Ngrid2*Ngrid3)
                    n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
                    n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3
                    p1 = n1 + (m1 + CWF_Plot_SuperCells1)*Ngrid1
                    p2 = n2 + (m2 + CWF_Plot_SuperCells2)*Ngrid2
                    p3 = n3 + (m3 + CWF_Plot_SuperCells3)*Ngrid3

                    gridn += 1
                    orbital_indices[gridn] = Nc
                    output_indices[gridn] = p1*output_n2*output_n3 + p2*output_n3 + p3 + 1
                end
            end

            resize!(orbital_indices, gridn)
            resize!(output_indices, gridn)
            grid_counts[atom] = gridn
            grid_orbitals[atom] = orbital_indices
            grid_indices[atom] = output_indices
        end

        CWF_GridN_Atom[cell] = grid_counts
        CWF_GridOrbs_Grid[cell] = grid_orbitals
        CWF_GridListAtom[cell] = grid_indices
    end


    return CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid
end


function _Calc_CWF_Grid8!(
    cwf_proj::Integer,
    NO0::Integer,
    CWF_GridN_Atom::Integer,
    CWF_GridOrbs_Grid::AbstractVector{<:Integer},
    CWF_GridListAtom::AbstractVector{<:Integer},
    Orbs_Grid::AbstractMatrix{Float64},
    ExpnCoef::AbstractVector{T},
    Wannier_Orbs_Grid::AbstractVector{T},
) where {T<:Union{Float64,ComplexF64}}

    @inbounds for Noc = 1:8:CWF_GridN_Atom-7

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

        temp0 = zero(T)
        temp1 = zero(T)
        temp2 = zero(T)
        temp3 = zero(T)
        temp4 = zero(T)
        temp5 = zero(T)
        temp6 = zero(T)
        temp7 = zero(T)
        for ist = 1:NO0
            Coef = ExpnCoef[cwf_proj+ist]
            temp0 += Coef*Orbs_Grid[ist,Nc0]
            temp1 += Coef*Orbs_Grid[ist,Nc1]
            temp2 += Coef*Orbs_Grid[ist,Nc2]
            temp3 += Coef*Orbs_Grid[ist,Nc3]
            temp4 += Coef*Orbs_Grid[ist,Nc4]
            temp5 += Coef*Orbs_Grid[ist,Nc5]
            temp6 += Coef*Orbs_Grid[ist,Nc6]
            temp7 += Coef*Orbs_Grid[ist,Nc7]
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
    @inbounds for Nog = 1:rem_CWF_GridN_Atom
        Nc = CWF_GridOrbs_Grid[Nog+Nog1]
		GN = CWF_GridListAtom[Nog+Nog1]

        Sum = zero(T)
        for ist = 1:NO0
            Sum += ExpnCoef[cwf_proj+ist]*Orbs_Grid[ist,Nc]
        end
        Wannier_Orbs_Grid[GN] += Sum
    end
    return nothing
end
