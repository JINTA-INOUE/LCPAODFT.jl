function Calc_dipole_moment(SpinPol::String, system_grid::System_Grid, Core_Charge, system_charge, Density_Grid)

    if SpinPol == "off"
        return Calc_dipole_moment_nospin(system_grid, Core_Charge, system_charge, Density_Grid)
    elseif SpinPol ∈ ("on", "nc")
        return Calc_dipole_moment_spin(system_grid, Core_Charge, system_charge, Density_Grid)
    else
        error("please check SpinPol")
    end
end


function Calc_dipole_moment_nospin(system_grid::System_Grid, Core_Charge, system_charge, Density_Grid)

    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    Grid_Origin = system_grid.Grid_Origin
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    
    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    Cell_Volume = abs(det(Latvecs))
    GridVol = Cell_Volume/NN

    E_dpx = 0.0
    E_dpy = 0.0
    E_dpz = 0.0

    E_dpx_BG = 0.0
    E_dpy_BG = 0.0
    E_dpz_BG = 0.0


    for i = 1:NN

        GN = i - 1
        n1 = div(GN, Ngrid2*Ngrid3)
        n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
        n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

        x = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
        y = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
        z = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

        den = 2.0*Density_Grid[1][i]

        E_dpx += den*x
        E_dpy += den*y
        E_dpz += den*z

        E_dpx_BG += x
        E_dpy_BG += y
        E_dpz_BG += z
    end

    E_dpx = E_dpx*GridVol
    E_dpy = E_dpy*GridVol
    E_dpz = E_dpz*GridVol

    cden_BG = system_charge/Cell_Volume

    E_dpx_BG = E_dpx_BG*GridVol*cden_BG
    E_dpy_BG = E_dpy_BG*GridVol*cden_BG
    E_dpz_BG = E_dpz_BG*GridVol*cden_BG


    C_dpx = 0.0
    C_dpy = 0.0
    C_dpz = 0.0

    for atom = 1:Natom
        x = Gxyz[atom][1]
        y = Gxyz[atom][2]
        z = Gxyz[atom][3]

        charge = Core_Charge[atom]
        C_dpx += charge*x
        C_dpy += charge*y
        C_dpz += charge*z
    end


    AU2Debye = 2.54174776


    dipole_moment = zeros(Float64, 4, 3)
    dipole_moment[1,1] = AU2Debye*(C_dpx - E_dpx - E_dpx_BG)
    dipole_moment[1,2] = AU2Debye*(C_dpy - E_dpy - E_dpy_BG)
    dipole_moment[1,3] = AU2Debye*(C_dpz - E_dpz - E_dpz_BG)

    dipole_moment[2,1] = AU2Debye*C_dpx
    dipole_moment[2,2] = AU2Debye*C_dpy
    dipole_moment[2,3] = AU2Debye*C_dpz

    dipole_moment[3,1] = -AU2Debye*E_dpx
    dipole_moment[3,2] = -AU2Debye*E_dpy
    dipole_moment[3,3] = -AU2Debye*E_dpz

    dipole_moment[4,1] = -AU2Debye*E_dpx_BG
    dipole_moment[4,2] = -AU2Debye*E_dpy_BG
    dipole_moment[4,3] = -AU2Debye*E_dpz_BG


    return dipole_moment
end


function Calc_dipole_moment_spin(system_grid::System_Grid, Core_Charge, system_charge, Density_Grid)


    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    Grid_Origin = system_grid.Grid_Origin
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3

    Cell_Volume = abs(det(Latvecs))
    GridVol = Cell_Volume/NN

    E_dpx = 0.0
    E_dpy = 0.0
    E_dpz = 0.0

    E_dpx_BG = 0.0
    E_dpy_BG = 0.0
    E_dpz_BG = 0.0


    for i = 1:NN

        GN = i - 1
        n1 = div(GN, Ngrid2*Ngrid3)
        n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
        n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

        x = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
        y = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
        z = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

        den = Density_Grid[1][i] + Density_Grid[2][i]

        E_dpx += den*x
        E_dpy += den*y
        E_dpz += den*z

        E_dpx_BG += x
        E_dpy_BG += y
        E_dpz_BG += z
    end

    E_dpx = E_dpx*GridVol
    E_dpy = E_dpy*GridVol
    E_dpz = E_dpz*GridVol

    cden_BG = system_charge/Cell_Volume

    E_dpx_BG = E_dpx_BG*GridVol*cden_BG
    E_dpy_BG = E_dpy_BG*GridVol*cden_BG
    E_dpz_BG = E_dpz_BG*GridVol*cden_BG


    C_dpx = 0.0
    C_dpy = 0.0
    C_dpz = 0.0

    for atom = 1:Natom
        x = Gxyz[atom][1]
        y = Gxyz[atom][2]
        z = Gxyz[atom][3]

        charge = Core_Charge[atom]
        C_dpx += charge*x
        C_dpy += charge*y
        C_dpz += charge*z
    end


    AU2Debye = 2.54174776


    dipole_moment = zeros(Float64, 4, 3)
    dipole_moment[1,1] = AU2Debye*(C_dpx - E_dpx - E_dpx_BG)
    dipole_moment[1,2] = AU2Debye*(C_dpy - E_dpy - E_dpy_BG)
    dipole_moment[1,3] = AU2Debye*(C_dpz - E_dpz - E_dpz_BG)

    dipole_moment[2,1] = AU2Debye*C_dpx
    dipole_moment[2,2] = AU2Debye*C_dpy
    dipole_moment[2,3] = AU2Debye*C_dpz

    dipole_moment[3,1] = -AU2Debye*E_dpx
    dipole_moment[3,2] = -AU2Debye*E_dpy
    dipole_moment[3,3] = -AU2Debye*E_dpz

    dipole_moment[4,1] = -AU2Debye*E_dpx_BG
    dipole_moment[4,2] = -AU2Debye*E_dpy_BG
    dipole_moment[4,3] = -AU2Debye*E_dpz_BG


    return dipole_moment
end