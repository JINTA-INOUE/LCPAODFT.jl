using LCPAODFT

function Al()

    Latvecs = [ 2.0250 0.0000 2.0250;
                2.0250 2.0250 0.0000;
                0.0000 2.0250 2.0250]*"Ang"
    atomorb = ["Al7.0-s2p2d1"]
    atomsymbol = ["Al"]
    atompos = [[0.0, 0.0, 0.0]]*"Frac"

    Ecut = 250.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 300
    xc_type = "GGA-PBE"
    kmesh = (9,9,9)

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

Al()

filepath = "Al.jld2"
kpath = [[0.5,0.75,0.25], [0.5,0.5,0.5], [0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.75,0.25], [0.375,0.75,0.375]]
kname = ["W", "L", "G", "X", "W", "K"]
# Band_kpath(filepath, kpath, kname)