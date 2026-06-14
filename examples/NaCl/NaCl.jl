using LCPAODFT

function NaCl()

    Latvecs = [ 0.00000   2.81500   2.81500;
                2.81500   0.00000   2.81500;
                2.81500   2.81500   0.00000]*"Ang"
    atomorb = ["Na9.0-s3p2d1", "Cl7.0-s2p2d1f1"]
    atomsymbol = ["Na", "Cl"]
    atompos = [[0.05, 0.0, 0.0], [2.81500, 2.81500, 2.81500]]*"Ang"

    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (9,9,9)

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

NaCl()

filepath = "NaCl.jld2"
kpath = [[0.0,0.0,0.0], [0.5,0.0,0.5], [0.5,0.25,0.75], [0.375,0.375,0.75], [0.0,0.0,0.0], [0.5,0.5,0.5], [0.625,0.25,0.625], [0.5,0.25,0.75], [0.5,0.5,0.5], [0.375,0.375,0.75]]
kname = ["G", "X", "W", "K", "G", "L", "U", "W", "L", "K"]
# Band_kpath(filepath, kpath, kname)