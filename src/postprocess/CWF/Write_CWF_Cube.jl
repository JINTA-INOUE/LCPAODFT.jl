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


const _CWF_CUBE_FORMAT1 = Printf.Format("%13.3E\n")
const _CWF_CUBE_FORMAT2 = Printf.Format("%13.3E%13.3E\n")
const _CWF_CUBE_FORMAT3 = Printf.Format("%13.3E%13.3E%13.3E\n")
const _CWF_CUBE_FORMAT4 = Printf.Format("%13.3E%13.3E%13.3E%13.3E\n")
const _CWF_CUBE_FORMAT5 = Printf.Format("%13.3E%13.3E%13.3E%13.3E%13.3E\n")
const _CWF_CUBE_FORMAT6 = Printf.Format("%13.3E%13.3E%13.3E%13.3E%13.3E%13.3E\n")

function Write_Wannier_Orbs_Grid(
    data,
    Plot_SuperCells,
    Ngrid,
    Wannier_Orbs_Grid,
    component=identity,
)
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    output_n1 = (2*Plot_SuperCells[1]+1)*Ngrid1
    output_n2 = (2*Plot_SuperCells[2]+1)*Ngrid2
    output_n3 = (2*Plot_SuperCells[3]+1)*Ngrid3
    required_values = output_n1*output_n2*output_n3
    length(Wannier_Orbs_Grid) >= required_values || throw(DimensionMismatch(
        "Wannier grid has $(length(Wannier_Orbs_Grid)) values; $required_values are required",
    ))

    # A cube row contains six values per line. Reuse one byte buffer for all
    # rows so Printf does not allocate a temporary StringVector per grid point.
    row_buffer = Vector{UInt8}(undef, 13*output_n3+cld(output_n3, 6))
    @inbounds for l = 0:output_n1-1, m = 0:output_n2-1
        offset = (l*output_n2 + m)*output_n3
        n = 0
        pos = 1
        while n+6 <= output_n3
            first = offset+n+1
            pos = Printf.format(
                row_buffer,
                pos,
                _CWF_CUBE_FORMAT6,
                component(Wannier_Orbs_Grid[first]),
                component(Wannier_Orbs_Grid[first+1]),
                component(Wannier_Orbs_Grid[first+2]),
                component(Wannier_Orbs_Grid[first+3]),
                component(Wannier_Orbs_Grid[first+4]),
                component(Wannier_Orbs_Grid[first+5]),
            )
            n += 6
        end

        remaining = output_n3-n
        first = offset+n+1
        if remaining == 1
            pos = Printf.format(row_buffer, pos, _CWF_CUBE_FORMAT1,
                component(Wannier_Orbs_Grid[first]))
        elseif remaining == 2
            pos = Printf.format(row_buffer, pos, _CWF_CUBE_FORMAT2,
                component(Wannier_Orbs_Grid[first]), component(Wannier_Orbs_Grid[first+1]))
        elseif remaining == 3
            pos = Printf.format(row_buffer, pos, _CWF_CUBE_FORMAT3,
                component(Wannier_Orbs_Grid[first]), component(Wannier_Orbs_Grid[first+1]),
                component(Wannier_Orbs_Grid[first+2]))
        elseif remaining == 4
            pos = Printf.format(row_buffer, pos, _CWF_CUBE_FORMAT4,
                component(Wannier_Orbs_Grid[first]), component(Wannier_Orbs_Grid[first+1]),
                component(Wannier_Orbs_Grid[first+2]), component(Wannier_Orbs_Grid[first+3]))
        elseif remaining == 5
            pos = Printf.format(row_buffer, pos, _CWF_CUBE_FORMAT5,
                component(Wannier_Orbs_Grid[first]), component(Wannier_Orbs_Grid[first+1]),
                component(Wannier_Orbs_Grid[first+2]), component(Wannier_Orbs_Grid[first+3]),
                component(Wannier_Orbs_Grid[first+4]))
        end
        GC.@preserve row_buffer unsafe_write(data, pointer(row_buffer), UInt(pos-1))
    end
end
