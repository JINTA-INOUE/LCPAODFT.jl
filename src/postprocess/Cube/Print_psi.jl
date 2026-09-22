@timeit timer "Print_psi" function Print_psi(filename::String, kpts, material::LCPAO_model, ucell::UCell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)

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


@timeit timer "Print_hwfs" function Print_hwfs(filename::String, kpts, material::LCPAO_model, ucell::UCell, Orbs_Grid, HOMOs_Coef)

    SpinPol = material.SpinPol
    Spindeg = ifelse(SpinPol=="off", 1, 2)
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

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    elseif SpinPol == "nc"
        Nocc = Int64(Valence_Electrons)
    end
    NHOMO = Nocc
    @show Nkpt, Spindeg, NHOMO


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

            kRn = kpts[ik][1]*l1 + kpts[ik][2]*l2 + kpts[ik][3]*l3
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


function Print_psi_1D(filename::String, kpts, material::LCPAO_model, ucell::UCell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)

    SpinPol = material.SpinPol
    Latvecs = material.Latvecs
    Natom = material.Natom
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



    RMO_Grid_3D = zeros(ComplexF64, Ngrid1, Ngrid2, Ngrid3)
    HOMOs_Coef_tmp = Vector{Vector{ComplexF64}}(undef, Natom)
    for atom = 1:Natom
        HOMOs_Coef_tmp[atom] = zeros(ComplexF64, Total_NumOrbs[atom])
    end

    for ik = 1:Nkpt, spin = 1:Spindeg, μ = 1:1
        μ1 = Bulk_HOMO[ik][spin]
        @. HOMOs_Coef_tmp = HOMOs_Coef[ik][spin][μ]
        fill!(RMO_Grid_3D, 0.0)
        for atom = 1:Natom, Nc = 1:GridN_Atom[atom]
            GN = GridListAtom[atom][Nc]
            Rn = CellListAtom[atom][Nc]+1

            n1 = div(GN, Ngrid2*Ngrid3)
            n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
            n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

            l1 = -atv_ijk[Rn][1]
            l2 = -atv_ijk[Rn][2]
            l3 = -atv_ijk[Rn][3]

            kRn = kpts[ik][1]*l1 + kpts[ik][2]*l2 + kpts[ik][3]*l3
            ex = exp(im*2*pi*kRn)

            for ist = 1:Total_NumOrbs[atom]
                RMO_Grid_3D[n1+1,n2+1,n3+1] += ex*HOMOs_Coef_tmp[atom][ist]*Orbs_Grid[atom][Nc][ist]
            end
        end

        
        file = open(filename*".homo$(ik-1)_$(spin-1)_$(μ-1)_a.dat", "w")
        Print_CubeData_1DTitle(file, kpts[ik], Enk[spin][ik][μ1], ChemP)
        Print_CubeData_1D_a(file, RMO_Grid_3D, gLatvecs, Ngrid)
        close(file)
        file = open(filename*".homo$(ik-1)_$(spin-1)_$(μ-1)_b.dat", "w")
        Print_CubeData_1DTitle(file, kpts[ik], Enk[spin][ik][μ1], ChemP)
        Print_CubeData_1D_b(file, RMO_Grid_3D, gLatvecs, Ngrid)
        close(file)
        file = open(filename*".homo$(ik-1)_$(spin-1)_$(μ-1)_c.dat", "w")
        Print_CubeData_1DTitle(file, kpts[ik], Enk[spin][ik][μ1], ChemP)
        Print_CubeData_1D_c(file, RMO_Grid_3D, gLatvecs, Ngrid)
        close(file)
    end
end


function Print_CubeData_1D_a(file, RMO_Grid::Array{ComplexF64,3}, gLatvecs, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    for i = 1:Ngrid1
        Sum = 0.0
        for j = 1:Ngrid2, k = 1:Ngrid3
            Sum += conj(RMO_Grid[i,j,k])*RMO_Grid[i,j,k]
        end

        a = (i-1)*gLatvecs[1,1]
        @printf(file, "%5.15f %5.15f\n", a, Sum/Ngrid2/Ngrid3)
    end
end


function Print_CubeData_1D_b(file, RMO_Grid::Array{ComplexF64,3}, gLatvecs, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    for j = 1:Ngrid2
        Sum = 0.0
        for i = 1:Ngrid1, k = 1:Ngrid3
            Sum += conj(RMO_Grid[i,j,k])*RMO_Grid[i,j,k]
        end

        b = (j-1)*gLatvecs[2,2]
        @printf(file, "%5.15f %5.15f\n", b, Sum/Ngrid1/Ngrid3)
    end
end


function Print_CubeData_1D_c(file, RMO_Grid::Array{ComplexF64,3}, gLatvecs, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    for k = 1:Ngrid3
        Sum = 0.0
        for i = 1:Ngrid1, j = 1:Ngrid2
            Sum += conj(RMO_Grid[i,j,k])*RMO_Grid[i,j,k]
        end

        c = (k-1)*gLatvecs[3,3]
        @printf(file, "%5.15f %5.15f\n", c, Sum/Ngrid1/Ngrid2)
    end
end