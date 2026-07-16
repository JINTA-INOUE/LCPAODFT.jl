using LCPAODFT

function Si()

    Latvecs = [ 2.7150  2.715  0.000;
                2.7150  0.000  2.715;
                0.0000  2.715  2.715]*"Ang"
    atomorb = ["Si7.0-s2p2d1"]
    atomsymbol = ["Si", "Si"]
    atompos = [[0.0,0.0,0.0], [0.25,0.25,0.25]]*"Frac"
    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (7,7,7)

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

Si()

kpath = [[0.5,0.75,0.25], [0.5,0.5,0.5], [0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.75,0.25], [0.375,0.75,0.375]]
kname = ["W", "L", "G", "X", "W", "K"]
# Band_kpath("Si.jld2", kpath, kname)
