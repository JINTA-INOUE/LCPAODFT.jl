function Print_CubeTitle(file::IOStream, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    @printf(file, " SYS1\n SYS1\n")
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Natom, Grid_Origin[1], Grid_Origin[2], Grid_Origin[3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid1, gLatvecs[1,1], gLatvecs[1,2], gLatvecs[1,3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid2, gLatvecs[2,1], gLatvecs[2,2], gLatvecs[2,3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid3, gLatvecs[3,1], gLatvecs[3,2], gLatvecs[3,3])

    for atom = 1:Natom
        Znum = Atom_Znumber[Atoms_Symbol[atom]]
        @printf(file, "%5d%12.6lf%12.6lf%12.6lf%12.6lf\n", Znum, 0.0, Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3])
    end
end


function Print_CubeTitle_psi(file::IOStream, kpts, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid, EigenValue, ChemP)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    @printf(file, "Absolute eigenvalue=%10.7f (Hartree)  Relative eigenvalue=%10.7f (Hartree)\n", EigenValue, EigenValue-ChemP)
    @printf(file, "Chemical Potential=%10.7f (Hartree)\n", ChemP)
    # @printf(file, "kpoints=%10.7f %10.7f %10.7f\n", kpts[1], kpts[2], kpts[3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Natom, Grid_Origin[1], Grid_Origin[2], Grid_Origin[3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid1, gLatvecs[1,1], gLatvecs[1,2], gLatvecs[1,3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid2, gLatvecs[2,1], gLatvecs[2,2], gLatvecs[2,3])
    @printf(file, "%5d%12.6lf%12.6lf%12.6lf\n", Ngrid3, gLatvecs[3,1], gLatvecs[3,2], gLatvecs[3,3])

    for atom = 1:Natom
        Znum = Atom_Znumber[Atoms_Symbol[atom]]
        @printf(file, "%5d%12.6lf%12.6lf%12.6lf%12.6lf\n", Znum, 0.0, Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3])
    end
end


function Print_CubeData_1DTitle(file::IOStream, kpts, EigenValue, ChemP)
    @printf(file, "# Absolute eigenvalue=%10.7f (Hartree)  Relative eigenvalue=%10.7f (Hartree)\n", EigenValue, EigenValue-ChemP)
    @printf(file, "# Chemical Potential=%10.7f (Hartree)\n", ChemP)
    @printf(file, "# kpoints=%10.7f %10.7f %10.7f\n", kpts[1], kpts[2], kpts[3])
end


function Print_CubeData(file::IOStream, data, Ngrid)
    
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)

    for ij = 0:Ngrid3:NN-1
        for k = 0:Ngrid3-1
            @printf(file, "%13.3E", data[ij+k+1])

            if iszero((k+1)%6)
                @printf(file, "\n")
            end
        end

        if Ngrid3%6 ≠ 0
            @printf(file, "\n")
        end
    end
end



function Print_CubeCData_MO(file1::IOStream, file2::IOStream, data::Vector{ComplexF64}, Ngrid)
    
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    for i = 0:Ngrid1-1, j = 0:Ngrid2-1
        for k = 0:Ngrid3-1
            GN = i*Ngrid2*Ngrid3 + j*Ngrid3 + k
            @printf(file1, "%13.3E", real(data[GN+1]))
            @printf(file2, "%13.3E", imag(data[GN+1]))

            if iszero((k+1)%6)
                @printf(file1, "\n")
                @printf(file2, "\n")
            end
        end

        if Ngrid3%6 ≠ 0
            @printf(file1, "\n")
            @printf(file2, "\n")
        end
    end
end