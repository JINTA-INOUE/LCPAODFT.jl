using LCPAODFT

function Fe_spin()

    Latvecs = [ 1.4350   1.4350  -1.4350;
               -1.4350   1.4350   1.4350;
                1.4350  -1.4350   1.4350]*"Ang"
    atomorb = ["Fe6.0H-s3p2d1"]
    atomsymbol = ["Fe"]
    atompos = [[0.0, 0.0, 0.0]]*"Ang"
    Atoms_Nspin = [[9.0, 7.0]]

    Ecut = 200.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (11,11,11)
    SpinPol = "on"

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Atoms_Nspin, Ecut, SCF_max, SCF_criterion, xc_type, kmesh, SpinPol)
    DFT(dft_setup)
end

Fe_spin()

filepath = "Fe_spin.jld2"
kpath = [[0.0,0.0,0.0], [-0.5,0.5,0.5], [0.0,0.0,0.5], [0.0,0.0,0.0], [0.25,0.25,0.25]]
kname = ["G", "H", "N", "G", "P"]
# Band_kpath(filepath, kpath, kname; Nk=50)
