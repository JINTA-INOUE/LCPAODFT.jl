using LCPAODFT

function Graphene()

    Latvecs = [0.0000000000000 2.489015870  0.0; 
               2.1555509738426 1.244507935  0.0; 
               0.0000000000000 0.000000000 10.0]*"Ang"
    atomorb = ["C6.0-s2p2d1"]
    atomsymbol = ["C", "C"]
    atompos = [[0.0, 0.0, 0.0], [2/3, 2/3, 0.0]]*"Frac"

    Ecut = 250.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 300
    xc_type = "LDA"
    kmesh = (13,13,1)

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

Graphene()

kpath = [[0.0,0.0,0.0], [0.5,0.0,0.0], [2/3,1/3,0.0], [0.0,0.0,0.0]]
kname = ["G", "M", "K", "G"]
# Band_kpath("Graphene.jld2", kpath, kname)