function Print_Density(filename, SpinPol, Atoms_Symbol, system_grid, ADensity_Grid, Density_Grid)

    println("<Print_Density> $filename.dden.cube")

    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    Grid_Origin = system_grid.Grid_Origin
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    

    Density = zeros(Float64, prod(Ngrid))
    if SpinPol == "off"
        @. Density = 2*Density_Grid[1] - 2*ADensity_Grid
    else
        @. Density = Density_Grid[1] + Density_Grid[2] - 2*ADensity_Grid
    end

    file = open(filename*".dden.cube", "w")
    Print_CubeTitle(file, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)
    Print_CubeData(file, Density, Ngrid)
    close(file)




    fill!(Density, 0.0)

    println("<Print_Density> $filename.tden.cube")
    if SpinPol == "off"
        @. Density = 2*Density_Grid[1]
    else
        @. Density = Density_Grid[1] + Density_Grid[2]
    end
    file = open(filename*".tden.cube", "w")
    Print_CubeTitle(file, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)
    Print_CubeData(file, Density, Ngrid)
    close(file)




    if SpinPol ∈ ("on", "nc")
        println("<Print_Density> $filename.den0.cube")
        file = open(filename*".tden.cube", "w")
        Print_CubeTitle(file, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)
        Print_CubeData(file, Density_Grid[1], Ngrid)
        close(file)

        println("<Print_Density> $filename.den1.cube")
        file = open(filename*".tden.cube", "w")
        Print_CubeTitle(file, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)
        Print_CubeData(file, Density_Grid[2], Ngrid)
        close(file)

        println("<Print_Density> $filename.sden.cube")
        file = open(filename*".tden.cube", "w")
        Print_CubeTitle(file, Natom, Gxyz, gLatvecs, Grid_Origin, Atoms_Symbol, Ngrid)
        Print_CubeData(file, Density_Grid[2]-Density_Grid[1], Ngrid)
        close(file)
    end
end