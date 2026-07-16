using LCPAODFT

function Cdia()

    Latvecs = [1.78  1.78  0.00;
               1.78  0.00  1.78;
               0.00  1.78  1.78]*"Ang"
    atomorb = ["C5.0-s2p2d1"]
    atomsymbol = ["C", "C"]
    atompos = [[0.0, 0.0, 0.0], [0.89, 0.89, 0.89]]*"Ang"
    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "LDA"
    kmesh = (7,7,7)
    fileout = true

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, fileout)
    DFT(dft_setup)
end

Cdia()

kpath = [[0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.75,0.25], [0.5,0.5,0.5], [0.0,0.0,0.0], [0.5,0.5,0.0]]
kname = ["G", "X", "W", "L", "G", "X"]
# Band_kpath("Cdia.jld2", kpath, kname)
