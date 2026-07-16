using LCPAODFT

function monoFeCl2_spin()

    Latvecs = [  3.4749999046    0.0000000000     0.0000000000;
                -1.7374999523    3.0094381956     0.0000000000;
                 0.0000000000    0.0000000000    17.2600002289]*"Ang"
    atomorb = ["Fe6.0S-s2p3d3f1", "Cl7.0-s3p3d2"]
    atomsymbol = ["Fe", "Cl", "Cl"]
    atompos = [[0.000000000, 0.000000000, 0.000000000],
               [0.333330026, 0.666670065, 0.079930000],
               [0.666670044, 0.333330002, 0.920069993]]*"Frac"
    Atoms_Nspin = [[8.0,6.0], [3.5,3.5], [3.5,3.5]]

    Ecut = 300.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 300
    xc_type = "GGA-PBE"
    kmesh = (13,13,1)
    SpinPol = "on"

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; SpinPol, Ecut, Atoms_Nspin, SCF_max, SCF_criterion, xc_type, kmesh)
    DFT(dft_setup)
end

# monoFeCl2_spin()
