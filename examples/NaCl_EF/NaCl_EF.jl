using LCPAODFT

function NaCl_EF()

    Latvecs = [ 0.00000   2.81500   2.81500;
                2.81500   0.00000   2.81500;
                2.81500   2.81500   0.00000]*"Ang"
    atomorb = ["Na9.0-s3p2", "Cl7.0-s3p2d2"]
    atomsymbol = ["Na", "Cl"]
    atompos = [[0.30, 0.0, 0.0], [2.815, 2.815, 2.815]]*"Ang"

    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-5
    SCF_max = 500
    xc_type = "GGA-PBE"
    kmesh = (5,5,5)

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh)

    Geo_Opt_Max = 100
    M_GDIIS_HISTORY = 7
    OptStartDIIS = 5
    OptEveryDIIS = 10
    Geo_Opt_criterion = 0.0001
    geoopt_setup = GeoOpt_Setup(dft_setup; Geo_Opt_criterion, Geo_Opt_Max, M_GDIIS_HISTORY, OptStartDIIS, OptEveryDIIS)
    DFT(geoopt_setup)
end

NaCl_EF()
