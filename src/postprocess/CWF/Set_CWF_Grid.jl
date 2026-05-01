function Set_CWF_Grid(CWF_ExpnCoef, Orbs_Grid, ucell::UCell, material, 
    Guide_Symbol, Guide_Total_NumOrbs, CWF_Plot_atom, CWF_Plot_SuperCells, filename)


    Natom = material.Natom
    MP = material.MP
    atv_ijk = material.atv_ijk
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)

    Ang_to_bohr = juOpenMX.Ang_to_bohr
    Guide_Gxyz = material.Gxyz

    Latvecs = ucell.system_grid.Latvecs
    Ngrid = ucell.system_grid.Ngrid
    Grid_Origin = ucell.system_grid.Grid_Origin
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    GridVol = ucell.system_grid.GridVol
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        Plot_cell_ijk[cell] = zeros(Int32, 3)
    end
    cell = 0
    for l1 = -CWF_Plot_SuperCells[1]:CWF_Plot_SuperCells[1], l2 = -CWF_Plot_SuperCells[2]:CWF_Plot_SuperCells[2], l3 = -CWF_Plot_SuperCells[3]:CWF_Plot_SuperCells[3]
        cell += 1
        Plot_cell_ijk[cell] = [l1, l2, l3]
    end


    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
    end


    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end



    CWF_Plot_SuperCells1, CWF_Plot_SuperCells2, CWF_Plot_SuperCells3 = CWF_Plot_SuperCells
    @show CWF_TNumGrid = prod(Ngrid)*Plot_NCell


    CWF_GridN_Atom = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        CWF_GridN_Atom[cell] = zeros(Int32, Natom)
    end

    CWF_GridOrbs_Grid = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    CWF_GridListAtom = Vector{Vector{Vector{Int32}}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        CWF_GridOrbs_Grid[cell] = Vector{Vector{Int32}}(undef, Natom)
        CWF_GridListAtom[cell] = Vector{Vector{Int32}}(undef, Natom)
    end


    for cell = 1:Plot_NCell

        l1, l2, l3 = Plot_cell_ijk[cell]
        for atom = 1:Natom

            gridn = 0
            for Nc = 1:GridN_Atom[atom]

                GRc = CellListAtom[atom][Nc]+1

                m1 = l1 + atv_ijk[GRc][1]
                m2 = l2 + atv_ijk[GRc][2]
                m3 = l3 + atv_ijk[GRc][3]
                                                
                if abs(m1)<=CWF_Plot_SuperCells1 && abs(m2)<=CWF_Plot_SuperCells2 && abs(m3)<=CWF_Plot_SuperCells3
                    gridn += 1
                end
            end
            
            CWF_GridN_Atom[cell][atom] = gridn
            CWF_GridOrbs_Grid[cell][atom] = zeros(Int32, gridn)
            CWF_GridListAtom[cell][atom] = zeros(Int32, gridn)

            gridn = 0
            for Nc = 1:GridN_Atom[atom]

                GNc = GridListAtom[atom][Nc]
                GRc = CellListAtom[atom][Nc]+1

                n1 = div(GNc, Ngrid2*Ngrid3)
                n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
                n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3

                m1 = l1 + atv_ijk[GRc][1]
                m2 = l2 + atv_ijk[GRc][2]
                m3 = l3 + atv_ijk[GRc][3]
                        
                if abs(m1)<=CWF_Plot_SuperCells1 && abs(m2)<=CWF_Plot_SuperCells2 && abs(m3)<=CWF_Plot_SuperCells3

                    p1 = n1 + (m1 + CWF_Plot_SuperCells1)*Ngrid1
                    p2 = n2 + (m2 + CWF_Plot_SuperCells2)*Ngrid2
                    p3 = n3 + (m3 + CWF_Plot_SuperCells3)*Ngrid3

                    GN1 = p1*Ngrid2*(2*CWF_Plot_SuperCells2+1)*Ngrid3*(2*CWF_Plot_SuperCells3+1) + p2*Ngrid3*(2*CWF_Plot_SuperCells3+1) + p3 + 1

                    gridn += 1
        
                    CWF_GridOrbs_Grid[cell][atom][gridn] = Nc
                    CWF_GridListAtom[cell][atom][gridn] = GN1
                end
            end
        end
    end



    



    if SpinPol ∈ ("off", "on")
        
        ExpnCoef = zeros(Float64, Nfsize)
        Wannier_Orbs_Grid = zeros(Float64, CWF_TNumGrid)

        for spin = 1:spinsize, patom in CWF_Plot_atom, pst = 1:Guide_Total_NumOrbs[patom]

            fill!(Wannier_Orbs_Grid, 0.0)

            for cell = 1:Plot_NCell
                
                @. ExpnCoef = CWF_ExpnCoef[spin][patom][pst][cell]

                for atom = 1:Natom
                    
                    cwf_proj = MP[atom]
                    NO0 = Total_NumOrbs[atom]

                    _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom[cell][atom], CWF_GridOrbs_Grid[cell][atom], CWF_GridListAtom[cell][atom], Orbs_Grid[atom], ExpnCoef, Wannier_Orbs_Grid)
                end
            end

            #=
            data = open("$(filename)_CWF$(patom)_$(pst).cube", "w")
            Write_CWF_CubeInfo(data, Guide_Symbol[patom], [1,1,1], Grid_Origin, Natom, Gxyz, Latvecs, gLatvecs, Ngrid)
            Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, [1,1,1], Ngrid, Wannier_Orbs_Grid)
            close(data)
            =#

            
            data = open("$(filename)_CWF$(patom)_$(pst).cube", "w")
            Write_CWF_CubeInfo(data, Guide_Symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Guide_Gxyz, Latvecs, gLatvecs, Ngrid)
            Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, Wannier_Orbs_Grid)
            close(data)
            


            CWF_norm = check_CWF_norm_Col(CWF_Plot_SuperCells, GridVol, Ngrid, Wannier_Orbs_Grid)
            println("$patom $pst $CWF_norm")

            CWF_R2, CWF_Pos = Calc_Pos_Omega_Col(CWF_Plot_SuperCells, Latvecs, Guide_Gxyz[patom], GridVol, Grid_Origin, Ngrid, Wannier_Orbs_Grid)
            

            Originx = CWF_Plot_SuperCells[1]*Latvecs[1,1] + CWF_Plot_SuperCells[2]*Latvecs[2,1] + CWF_Plot_SuperCells[3]*Latvecs[3,1] - Grid_Origin[1]
            Originy = CWF_Plot_SuperCells[1]*Latvecs[1,2] + CWF_Plot_SuperCells[2]*Latvecs[2,2] + CWF_Plot_SuperCells[3]*Latvecs[3,2] - Grid_Origin[2]
            Originz = CWF_Plot_SuperCells[1]*Latvecs[1,3] + CWF_Plot_SuperCells[2]*Latvecs[2,3] + CWF_Plot_SuperCells[3]*Latvecs[3,3] - Grid_Origin[3]

            
            Pos_x = CWF_Pos[1]
            Pos_y = CWF_Pos[2]
            Pos_z = CWF_Pos[3]
            Pos_r2 = sum(CWF_R2)
            X = Pos_x - Originx
            Y = Pos_y - Originy
            Z = Pos_z - Originz
            X = X/Ang_to_bohr
            Y = Y/Ang_to_bohr
            Z = Z/Ang_to_bohr

            CWF_Omega = Pos_r2 - (Pos_x^2 + Pos_y^2 + Pos_z^2)
            omega = CWF_Omega/Ang_to_bohr/Ang_to_bohr
            @printf("%d %d %5.10f %5.10f %5.10f %5.10f\n", patom, pst, X, Y, Z, omega)
        end
    elseif SpinPol == "nc"
        CWF_norm = zeros(Float64, 2)
        CWF_Pos = zeros(Float64, 2, 3)
        CWF_R2 = zeros(Float64, 2)
        CWF_Omega = zeros(Float64, 2)
        ExpnCoef = zeros(ComplexF64, Nfsize)
        Wannier_Orbs_Grid = zeros(ComplexF64, CWF_TNumGrid)

        for patom in CWF_Plot_atom, pst = 1:2*Guide_Total_NumOrbs[patom]

            data = open("$(filename)_CWF$(patom)_$pst.nccube", "w")
            Write_CWF_CubeInfo(data, Guide_Symbol, CWF_Plot_SuperCells, Grid_Origin, Natom, Guide_Gxyz, Latvecs, gLatvecs, Ngrid)

            for spin = 1:2

                fill!(Wannier_Orbs_Grid, 0.0)
                spin_site = ifelse(spin==1, 0, fsize)

                for cell = 1:Plot_NCell

                    l1, l2, l3 = Plot_cell_ijk[cell]
                    @. ExpnCoef = CWF_ExpnCoef[patom][pst][cell]

                    for atom = 1:Natom, Nc = 1:GridN_Atom[atom]

                        cwf_proj = MP[atom]

                        GNc = GridListAtom[atom][Nc]
                        GRc = CellListAtom[atom][Nc]+1

                        n1 = div(GNc, Ngrid2*Ngrid3)
                        n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
                        n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3

                        m1 = l1 + atv_ijk[GRc][1]
                        m2 = l2 + atv_ijk[GRc][2]
                        m3 = l3 + atv_ijk[GRc][3]
                            
                        if abs(m1)<=CWF_Plot_SuperCells[1] && abs(m2)<=CWF_Plot_SuperCells[2] && abs(m3)<=CWF_Plot_SuperCells[3]
                            p1 = n1 + (m1 + CWF_Plot_SuperCells[1])*Ngrid1
                            p2 = n2 + (m2 + CWF_Plot_SuperCells[2])*Ngrid2
                            p3 = n3 + (m3 + CWF_Plot_SuperCells[3])*Ngrid3

                            GN1 = p1*Ngrid2*(2*CWF_Plot_SuperCells[2]+1)*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + p2*Ngrid3*(2*CWF_Plot_SuperCells[3]+1) + p3 + 1

                            temp = ComplexF64(0.0, 0.0)
                            @inbounds for ist = 1:Total_NumOrbs[atom]
                                temp += ExpnCoef[cwf_proj+ist+spin_site] * Orbs_Grid[atom][ist][Nc]
                            end
                            Wannier_Orbs_Grid[GN1] += temp
                        end
                    end
                end

                CWF_norm[spin] = check_CWF_norm_NonCol(CWF_Plot_SuperCells, GridVol, Ngrid, Wannier_Orbs_Grid)

                # CWF_Pos[spin,:], CWF_Omega[spin] = Calc_Pos_Omega_NonCol(CWF_Plot_SuperCells, Latvecs, Guide_Gxyz[patom], GridVol, Grid_Origin, Ngrid, Wannier_Orbs_Grid)
                CWF_R2[spin], CWF_Pos[spin,:] = temp_Calc_Pos_Omega_NonCol(CWF_Plot_SuperCells, Latvecs, Guide_Gxyz[patom], GridVol, Grid_Origin, Ngrid, Wannier_Orbs_Grid)
                

                Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, real(Wannier_Orbs_Grid))
                Write_Wannier_Orbs_Grid(data, CWF_Plot_SuperCells, Ngrid, imag(Wannier_Orbs_Grid))
            end

            Norm = sum(CWF_norm)
            
            
            omega = sum(CWF_Omega)
            println("$patom $pst $Norm")
            # println("$patom $pst $x  $y  $z  $omega")
            # @show CWF_R2, CWF_Pos[1,:], CWF_Pos[2,:]


            # Originx = CWF_Plot_SuperCells[1]*Latvecs[1,1] + CWF_Plot_SuperCells[2]*Latvecs[2,1] + CWF_Plot_SuperCells[3]*Latvecs[3,1] - Guide_Gxyz[patom][1] - Grid_Origin[1]
            # Originy = CWF_Plot_SuperCells[1]*Latvecs[1,2] + CWF_Plot_SuperCells[2]*Latvecs[2,2] + CWF_Plot_SuperCells[3]*Latvecs[3,2] - Guide_Gxyz[patom][2] - Grid_Origin[2]
            # Originz = CWF_Plot_SuperCells[1]*Latvecs[1,3] + CWF_Plot_SuperCells[2]*Latvecs[2,3] + CWF_Plot_SuperCells[3]*Latvecs[3,3] - Guide_Gxyz[patom][3] - Grid_Origin[3]

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
            @printf("%d %d %5.10f %5.10f %5.10f %5.10f\n", patom, pst, X, Y, Z, omega)
            # println("$patom $pst $X  $Y  $Z  $omega")
            
            close(data)
        end
    end
end


function _Calc_CWF_Grid8!(cwf_proj, NO0, CWF_GridN_Atom, CWF_GridOrbs_Grid, CWF_GridListAtom, Orbs_Grid, ExpnCoef, Wannier_Orbs_Grid)

    for Noc = 1:8:CWF_GridN_Atom-7

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

        temp0 = 0.0
        temp1 = 0.0
        temp2 = 0.0
        temp3 = 0.0
        temp4 = 0.0
        temp5 = 0.0
        temp6 = 0.0
        temp7 = 0.0
        for ist = 1:NO0
            Coef = ExpnCoef[cwf_proj+ist]
            temp0 += Coef * Orbs_Grid[ist][Nc0]
            temp1 += Coef * Orbs_Grid[ist][Nc1]
            temp2 += Coef * Orbs_Grid[ist][Nc2]
            temp3 += Coef * Orbs_Grid[ist][Nc3]
            temp4 += Coef * Orbs_Grid[ist][Nc4]
            temp5 += Coef * Orbs_Grid[ist][Nc5]
            temp6 += Coef * Orbs_Grid[ist][Nc6]
            temp7 += Coef * Orbs_Grid[ist][Nc7]
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
    for Nog = 1:rem_CWF_GridN_Atom
        Nc = CWF_GridOrbs_Grid[Nog+Nog1]
		GN = CWF_GridListAtom[Nog+Nog1]

        Sum = 0.0
        for ist = 1:NO0
            Sum += ExpnCoef[cwf_proj+ist] * Orbs_Grid[ist][Nc]
        end
        Wannier_Orbs_Grid[GN] += Sum
    end
end