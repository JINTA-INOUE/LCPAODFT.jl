using LCPAODFT

function GaAs_SOC()

    Latvecs = [5.367 0.000 5.367; 
               0.000 5.367 5.367; 
               5.367 5.367 0.000]*"AU"
    atomorb = ["Ga7.0-s2p2d1", "As7.0-s2p2d1"]
    atomsymbol = ["Ga", "As"]
    atompos = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]*"Frac"

    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-7
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (7,7,7)
    SpinPol = "nc"
    SO_switch = true

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; SpinPol, SO_switch, Ecut, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

GaAs_SOC()

filepath = "GaAs_SOC.jld2"
kpath = [[0.5,0.75,0.25], [0.5,0.5,0.5], [0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.75,0.25], [0.375,0.75,0.375]]
kname = ["W", "L", "G", "X", "W", "K"]
# Band_kpath(filepath, kpath, kname)