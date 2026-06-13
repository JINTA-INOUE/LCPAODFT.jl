function Write_Wannier_CubeInfo(data, Atoms_Symbol, Plot_SuperCells, Grid_Origin, Natom, Gxyz, Latvecs, gtv, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    Ncell = (2*Plot_SuperCells[1]+1)*(2*Plot_SuperCells[2]+1)*(2*Plot_SuperCells[3]+1)
    
    @printf(data, " SYS1\n SYS1\n")
    @printf(data, "%5d%12.6lf%12.6lf%12.6lf\n", Natom*Ncell,Grid_Origin[1],Grid_Origin[2],Grid_Origin[3])
    @printf(data, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid1*(2*Plot_SuperCells[1]+1),gtv[1,1],gtv[1,2],gtv[1,3])
    @printf(data, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid2*(2*Plot_SuperCells[2]+1),gtv[2,1],gtv[2,2],gtv[2,3])
    @printf(data, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid3*(2*Plot_SuperCells[3]+1),gtv[3,1],gtv[3,2],gtv[3,3])
        
        
    for l = 0:2*Plot_SuperCells[1]
        for m = 0:2*Plot_SuperCells[2]
            for n = 0:2*Plot_SuperCells[3]
                for atom = 1:Natom
                    Znum = LCPAODFT.Atom_Znumber[Atoms_Symbol[atom]]
                    
                    @printf(data, "%5d%12.6lf%12.6lf%12.6lf%12.6lf\n",
                    Znum, 0, 
                    Gxyz[atom][1]+l*Latvecs[1,1]+m*Latvecs[2,1]+n*Latvecs[3,1],
                    Gxyz[atom][2]+l*Latvecs[1,2]+m*Latvecs[2,2]+n*Latvecs[3,2],
                    Gxyz[atom][3]+l*Latvecs[1,3]+m*Latvecs[2,3]+n*Latvecs[3,3]) 
                end
            end
        end
    end
end


function Write_Wannier_Orbs_Grid(data, Plot_SuperCells, Ngrid, Wannier_Orbs_Grid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    for l = 0:(2*Plot_SuperCells[1]+1)*Ngrid1-1
        for m = 0:(2*Plot_SuperCells[2]+1)*Ngrid2-1
            for n = 0:(2*Plot_SuperCells[3]+1)*Ngrid3-1
                GN = l*Ngrid2*(2*Plot_SuperCells[2]+1)*Ngrid3*(2*Plot_SuperCells[3]+1) + m*Ngrid3*(2*Plot_SuperCells[3]+1) + n + 1

                @printf(data, "%13.3E", Wannier_Orbs_Grid[GN])
                if (n+1)%6 == 0
                    @printf(data, "\n")
                end
            end
            if ((2*Plot_SuperCells[3]+1)*Ngrid3)%6 ≠ 0
                @printf(data, "\n")
            end
        end 
    end
end