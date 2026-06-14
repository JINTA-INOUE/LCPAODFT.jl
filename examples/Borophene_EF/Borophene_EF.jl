using LCPAODFT

function Borophene_EF()

    Latvecs = [ 5.089926331914958   0.000000000000949  -0.000000000000000;
                0.000000000000000   2.930364635816637   0.000000000000000;
                0.000000000000000   0.000000000000000  25.000000000000000]*"Ang"
    atomorb = ["B7.0-s3p3d2"]
    atomsymbol = ["B", "B", "B", "B", "B"]

    atompos = [[1.66707155244092,    1.46960461606491,   12.49984981326827],
               [3.42339959792466,    1.46960671000884,   12.49997176770496],
               [0.82635433649220,    0.00442036928587,   12.49996683427698],
               [4.26416223007023,    0.00445402128485,   12.50021499930880],
               [2.54526336896233,    0.00445461868929,   12.49999400066109]]*"Ang"

    Ecut = 200.0
    system = "Crystal"
    SCF_criterion = 1e-7
    SCF_max = 300
    xc_type = "GGA-PBE"
    kmesh = (5,10,1)
    filename = "Borophene"

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, filename)
    
    Geo_Opt_Max = 100
    M_GDIIS_HISTORY = 7
    OptStartDIIS = 5
    OptEveryDIIS = 10
    Geo_Opt_criterion = 0.0001
    geoopt_setup = GeoOpt_Setup(dft_setup; Geo_Opt_criterion, Geo_Opt_Max, M_GDIIS_HISTORY, OptStartDIIS, OptEveryDIIS)
    DFT(geoopt_setup)
end

Borophene_EF()