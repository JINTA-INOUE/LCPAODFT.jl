function Print_psi(filename::String, kpts, material::LCPAO_model, ucell::UCell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)

    SpinPol = material.SpinPol
    Latvecs = material.Latvecs
    Natom = material.Natom
    Atoms_symbol = material.Atoms_symbol
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    atv_ijk = material.atv_ijk
    ChemP = material.ChemP
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Grid_Origin = ucell.system_grid.Grid_Origin
    Ngrid = ucell.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    Nkpt = length(kpts)

    if SpinPol ∈ ("on", "nc")
        Spindeg = 2
    elseif SpinPol == "off"
        Spindeg = 1
    end

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3



    RMO_Grid = zeros(ComplexF64, NN)    
    HOMOs_Coef_tmp = Vector{Vector{ComplexF64}}(undef, Natom)
    for atom = 1:Natom
        HOMOs_Coef_tmp[atom] = zeros(ComplexF64, Total_NumOrbs[atom])
    end

    for ik = 1:Nkpt, spin = 1:Spindeg, μ = 1:1
        μ1 = Bulk_HOMO[ik][spin]
        @. HOMOs_Coef_tmp = HOMOs_Coef[ik][spin][μ]
        fill!(RMO_Grid, 0.0)
        for atom = 1:Natom, Nc = 1:GridN_Atom[atom]

            GN = GridListAtom[atom][Nc]+1
            Rn = CellListAtom[atom][Nc]+1

            l1, l2, l3 = atv_ijk[Rn]

            kRn = kpts[ik][1]*l1 + kpts[ik][2]*l2 + kpts[ik][3]*l3
            ex = exp(-im*2*pi*kRn)

            for ist = 1:Total_NumOrbs[atom]
                RMO_Grid[GN] += ex*HOMOs_Coef_tmp[atom][ist]*Orbs_Grid[atom][Nc][ist]
            end
        end

        if SpinPol == "nc"
            EigenValues = Enk[1][ik][μ1]
        else
            EigenValues = Enk[spin][ik][μ1]
        end

        
        file1 = open(filename*".homo$(ik-1)_$(spin-1)_$(μ-1)_r.cube", "w")
        file2 = open(filename*".homo$(ik-1)_$(spin-1)_$(μ-1)_i.cube", "w")
        Print_CubeTitle_psi(file1, kpts[ik], Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_symbol, Ngrid, EigenValues, ChemP)
        Print_CubeTitle_psi(file2, kpts[ik], Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_symbol, Ngrid, EigenValues, ChemP)
        Print_CubeCData_MO(file1, file2, RMO_Grid, Ngrid)
        close(file1)
        close(file2)
    end
end


function Print_hwfs(filename::String, kpts, material::LCPAO_model, ucell::UCell, Orbs_Grid, HOMOs_Coef)

    SpinPol = material.SpinPol
    Latvecs = material.Latvecs
    Natom = material.Natom
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    atv_ijk = material.atv_ijk
    Atoms_symbol = material.Atoms_symbol
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Grid_Origin = ucell.system_grid.Grid_Origin
    Ngrid = ucell.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    Nkpt = length(kpts)

    if SpinPol ∈ ("on", "nc")
        Spindeg = 2
    elseif SpinPol == "off"
        Spindeg = 1
    end

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    elseif SpinPol == "nc"
        Nocc = Int64(Valence_Electrons)
    end
    NHOMO = Nocc


    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3



    RMO_Grid = zeros(ComplexF64, NN)

    for ik = 1:Nkpt, spin = 1:Spindeg, μ = 1:NHOMO
        fill!(RMO_Grid, 0.0)
        for atom = 1:Natom, Nc = 1:GridN_Atom[atom]
            GN = GridListAtom[atom][Nc]+1
            Rn = CellListAtom[atom][Nc]+1

            l1, l2, l3 = atv_ijk[Rn]

            kRn = kpts[1]*l1 + kpts[2]*l2 + kpts[3]*l3
            ex = exp(-im*2*pi*kRn)

            for ist = 1:Total_NumOrbs[atom]
                RMO_Grid[GN] += ex*HOMOs_Coef[spin][μ][atom][ist]*Orbs_Grid[atom][Nc][ist]
            end
        end

        
        file1 = open(filename*"_hwfs.homo$(ik-1)_$(spin-1)_$(μ-1)_r.cube", "w")
        file2 = open(filename*"_hwfs.homo$(ik-1)_$(spin-1)_$(μ-1)_i.cube", "w")
        Print_CubeTitle(file1, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_symbol, Ngrid)
        Print_CubeTitle(file2, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_symbol, Ngrid)
        Print_CubeCData_MO(file1, file2, RMO_Grid, Ngrid)
        close(file1)
        close(file2)
    end
end