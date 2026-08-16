function check_CWF_norm(CWF_Plot_SuperCells, GridVol, Ngrid, Wannier_Orbs_Grid)
    
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    CWF_norm = 0.0
    for l = 0:(2*CWF_Plot_SuperCells[1]+1)*Ngrid1-1
        for m = 0:(2*CWF_Plot_SuperCells[2]+1)*Ngrid2-1
            @inbounds for n = 0:(2*CWF_Plot_SuperCells[3]+1)*Ngrid3-1

                GN = l*Ngrid2*(2*CWF_Plot_SuperCells[2]+1)*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + m*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + n + 1

                wann = Wannier_Orbs_Grid[GN]
                CWF_norm += abs2(wann)
            end
        end
    end

    return CWF_norm*GridVol
end


function Calc_Pos_Omega(CWF_Plot_SuperCells, gLatvecs, GridVol, Ngrid, Wannier_Orbs_Grid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    Sumx = 0.0
    Sumy = 0.0
    Sumz = 0.0
    Sumr2 = 0.0

    for l = 0:(2*CWF_Plot_SuperCells[1]+1)*Ngrid1-1
        for m = 0:(2*CWF_Plot_SuperCells[2]+1)*Ngrid2-1
            @inbounds for n = 0:(2*CWF_Plot_SuperCells[3]+1)*Ngrid3-1
                GN = l*Ngrid2*(2*CWF_Plot_SuperCells[2]+1)*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + m*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + n + 1
            
                GNc = GN - 1
                temp2 = Ngrid2*(2*CWF_Plot_SuperCells[2]+1)
                temp3 = Ngrid3*(2*CWF_Plot_SuperCells[3]+1)
                n1 = div(GNc, temp2*temp3)
                n2 = div(GNc - n1*temp2*temp3, temp3)
                n3 = GNc - n1*temp2*temp3 - n2*temp3
                x = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1]
                y = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2]
                z = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3]
                
                wann = Wannier_Orbs_Grid[GN]
                wann2 = abs2(wann)
                Sumx += wann2 * x
                Sumy += wann2 * y
                Sumz += wann2 * z
                Sumr2 += wann2 * (x^2 + y^2 + z^2)
            end
        end
    end

    Pos_x = Sumx*GridVol
    Pos_y = Sumy*GridVol
    Pos_z = Sumz*GridVol
    Pos_r2 = Sumr2*GridVol
 
 
    return Pos_r2, Pos_x, Pos_y, Pos_z
end


function Calc_WannierCenter(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, filepath::String, CWF_Plot_Cube)

    data = jldopen(filepath, "r")
    SpinPol = data["SpinPol"]
    spinsize = data["spinsize"]
    Nfsize = data["Nfsize"]
    CWF_Plot_SuperCells = data["CWF_Plot_SuperCells"]
    CWF_ExpnCoef = data["CWF_ExpnCoef"]
    close(data)


    material = cwf_setup.material
    Nspin = material.Nspin
    Latvecs = material.Latvecs
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Gxyz = material.Gxyz
    TCpyCell = material.TCpyCell
    Atoms_pao = material.Atoms_pao
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    Ecut = cwf_setup.Ecut
    Ngsize = cwf_setup.Ngsize
    Wannier_Guide = cwf_setup.Wannier_Guide
    Ngrid = Calc_Ngrid(Ecut, Latvecs)
    Atoms_Cut1 = material.Atoms_Cut1
    Grid_Origin = material.Grid_Origin
    
    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)


    pao = Vector{PAO}(undef, Nspecies)
    for spe = 1:Nspecies
        pao[spe] = Read_PAO(0.0, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end
        
    ucell = UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs)
    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)

    GridVol = ucell.system_grid.GridVol
    Ngrid = ucell.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    CWF_GridN_Atom, CWF_GridListAtom, CWF_GridOrbs_Grid = CWF_UCell(cwf_setup, ucell)

    
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    CWF_TNumGrid = prod(Ngrid)*Plot_NCell
    

    println("\n")
    println("Atomic Orbital Center(Ang)")
    for p = 1:Ngsize
        @printf("  proj = %3d  %5.12f  %5.12f  %5.12f\n", p, Wannier_Guide[p][1]/Ang_to_bohr, Wannier_Guide[p][2]/Ang_to_bohr, Wannier_Guide[p][3]/Ang_to_bohr)
    end

    println("")
    if SpinPol ∈ ("off", "on")
        
        println("spin  orbital  norm   X(Ang^2)   Y(Ang^2)  Z(Ang^2)   omega(Ang^2)")
        ExpnCoef = zeros(Float64, Nfsize)
        Wannier_Orbs_Grid = zeros(Float64, CWF_TNumGrid)

        for spin = 1:spinsize, proj in CWF_Plot_Cube
            fill!(Wannier_Orbs_Grid, 0.0)
            for cell = 1:Plot_NCell
                @. ExpnCoef = CWF_ExpnCoef[spin][proj][cell]
                for atom = 1:Natom
                    
                    cwf_proj = MP[atom]
                    NO0 = Total_NumOrbs[atom]

                    _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid[atom], ExpnCoef, Wannier_Orbs_Grid)
                end
            end

            CWF_norm = check_CWF_norm(CWF_Plot_SuperCells, GridVol, Ngrid, Wannier_Orbs_Grid)
            Pos_r2, Pos_x, Pos_y, Pos_z = Calc_Pos_Omega(CWF_Plot_SuperCells, gLatvecs, GridVol, Ngrid, Wannier_Orbs_Grid)
            
            Originx = CWF_Plot_SuperCells[1]*Latvecs[1,1] + CWF_Plot_SuperCells[2]*Latvecs[2,1] + CWF_Plot_SuperCells[3]*Latvecs[3,1] - Grid_Origin[1]
            Originy = CWF_Plot_SuperCells[1]*Latvecs[1,2] + CWF_Plot_SuperCells[2]*Latvecs[2,2] + CWF_Plot_SuperCells[3]*Latvecs[3,2] - Grid_Origin[2]
            Originz = CWF_Plot_SuperCells[1]*Latvecs[1,3] + CWF_Plot_SuperCells[2]*Latvecs[2,3] + CWF_Plot_SuperCells[3]*Latvecs[3,3] - Grid_Origin[3]

            X = Pos_x - Originx
            Y = Pos_y - Originy
            Z = Pos_z - Originz
            X = X/Ang_to_bohr
            Y = Y/Ang_to_bohr
            Z = Z/Ang_to_bohr

            CWF_Omega = Pos_r2 - (Pos_x^2 + Pos_y^2 + Pos_z^2)
            omega = CWF_Omega/Ang_to_bohr/Ang_to_bohr
            @printf("%3d %3d   %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", spin, proj, CWF_norm, X, Y, Z, omega)
        end
    elseif SpinPol == "nc"

        println("orbital  norm   X(Ang^2)   Y(Ang^2)  Z(Ang^2)   omega(Ang^2)")
        CWF_norm = zeros(Float64, 2)
        CWF_Pos = zeros(Float64, 2, 3)
        CWF_R2 = zeros(Float64, 2)
        CWF_Omega = zeros(Float64, 2)
        ExpnCoef = zeros(ComplexF64, Nfsize)
        Wannier_Orbs_Grid = zeros(ComplexF64, CWF_TNumGrid)

        for proj in CWF_Plot_Cube

            for spin = 1:2
                fill!(Wannier_Orbs_Grid, 0.0)
                spin_site = ifelse(spin==1, 0, fsize)
                for cell = 1:Plot_NCell
                    @. ExpnCoef = CWF_ExpnCoef[proj][cell]
                    for atom = 1:Natom
                    
                        cwf_proj = MP[atom]
                        NO0 = Total_NumOrbs[atom]

                        _Calc_CWF_Grid8!(spin_site+cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid[atom], ExpnCoef, Wannier_Orbs_Grid)
                    end
                end

                CWF_norm[spin] = check_CWF_norm(CWF_Plot_SuperCells, GridVol, Ngrid, Wannier_Orbs_Grid)
                CWF_R2[spin], CWF_Pos[spin,1], CWF_Pos[spin,2], CWF_Pos[spin,3] = Calc_Pos_Omega(CWF_Plot_SuperCells, gLatvecs, GridVol, Ngrid, Wannier_Orbs_Grid)
            end

            Norm = sum(CWF_norm)            
            omega = sum(CWF_Omega)

            Originx = CWF_Plot_SuperCells[1]*Latvecs[1,1] + CWF_Plot_SuperCells[2]*Latvecs[2,1] + CWF_Plot_SuperCells[3]*Latvecs[3,1] - Grid_Origin[1]
            Originy = CWF_Plot_SuperCells[1]*Latvecs[1,2] + CWF_Plot_SuperCells[2]*Latvecs[2,2] + CWF_Plot_SuperCells[3]*Latvecs[3,2] - Grid_Origin[2]
            Originz = CWF_Plot_SuperCells[1]*Latvecs[1,3] + CWF_Plot_SuperCells[2]*Latvecs[2,3] + CWF_Plot_SuperCells[3]*Latvecs[3,3] - Grid_Origin[3]

            Pos_x = sum(CWF_Pos[:,1])
            Pos_y = sum(CWF_Pos[:,2])
            Pos_z = sum(CWF_Pos[:,3])
            Pos_r2 = sum(CWF_R2)
            X = Pos_x - Originx
            Y = Pos_y - Originy
            Z = Pos_z - Originz
            X = X/Ang_to_bohr
            Y = Y/Ang_to_bohr
            Z = Z/Ang_to_bohr

            CWF_Omega = Pos_r2 - (Pos_x^2 + Pos_y^2 + Pos_z^2)
            omega = CWF_Omega/Ang_to_bohr/Ang_to_bohr
            @printf("%3d   %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", proj, Norm, X, Y, Z, omega)
        end
    end
end
