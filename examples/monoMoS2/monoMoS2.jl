using LCPAODFT

function monoMoS2()

    Latvecs = [ 3.16   0.000000  0.0000;
               -1.58   2.736640  0.0000;
                0.00   0.0000000 12.345]*"Ang"
    atomorb = ["Mo7.0-s2p2d1", "S7.0-s2p2d1"]
    atomsymbol = ["Mo", "S", "S"]
    atompos = [[1.580000,  0.91221,   3.0725],
               [0.000000,  1.82443,   1.4795],
               [0.000000,  1.82443,   4.6655]]*"Ang"

    Ecut = 200.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 300
    xc_type = "GGA-PBE"
    kmesh = (9,9,1)
    filename = "monoMoS2"

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, filename)
    DFT(dft_setup)
end

monoMoS2()