using LCPAODFT

function Bi2Se3_SOC()

    Latvecs = [ 2.0715000     1.1959811     9.5453333;
               -2.0715000     1.1959811     9.5453333;
                0.0000000    -2.3919622     9.5453333]*"Ang"
    atomorb = ["Bi8.0-s3p2d2f1", "Se7.0-s3p2d2"]
    atomsymbol = ["Bi", "Bi", "Se", "Se", "Se"]
    atompos = [[ 0.4008,  0.4008,  0.4008],
               [ 0.5992,  0.5992,  0.5992],
               [ 0.0000,  0.0000,  0.0000],
               [ 0.2117,  0.2117,  0.2117],
               [ 0.7883,  0.7883,  0.7883]]*"Frac"
    Ecut = 250.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (7,7,7)
    SpinPol = "nc"
    SO_switch = true
    filename = "Bi2Se3_SOC"

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, SpinPol, SO_switch, filename)
    DFT(dft_setup)
end

Bi2Se3_SOC()

kpath = [[0.0,0.0,0.0],[0.5,0.5,0.5],[0.5,0.5,0.0],[0.0,0.0,0.0],[0.5,0.0,0.0]]
kname = ["G", "Z", "F", "G", "L"]
# Band_kpath("Bi2Se3_SOC.jld2", kpath, kname)