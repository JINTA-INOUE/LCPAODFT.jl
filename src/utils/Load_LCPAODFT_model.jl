function Load_LCPAODFT_model(jld2_file::String)

    data = jldopen(jld2_file, "r")
    Natom = data["Natom"]
    Nspecies = data["Nspecies"]
    Nspin = data["Nspin"]
    atom2spe = data["atom2spe"]
    Atoms_symbol = data["Atoms_symbol"]
    Atoms_Cut1 = data["Atoms_Cut1"]
    Atoms_pao = data["Atoms_pao"]
    Atoms_Core_Charge = data["Atoms_Core_Charge"]
    Init_Atoms_Nspin = data["Init_Atoms_Nspin"]
    Init_Atoms_Angle = data["Init_Atoms_Angle"]
    Atoms_Angle = data["Atoms_Angle"]
    Total_SpinS = data["Total_SpinS"]
    Latvecs = data["Latvecs"]
    Recvecs = data["Recvecs"]
    Gxyz = data["Gxyz"]
    TCpyCell = data["TCpyCell"]
    atv = data["atv"]
    atv_ijk = data["atv_ijk"]
    FNAN = data["FNAN"]
    natn = data["natn"]
    ncn = data["ncn"]
    Total_NumOrbs = data["Total_NumOrbs"]
    MP = data["MP"]
    Grid_Origin = data["Grid_Origin"]
    Ngrid = data["Ngrid"]
    SO_switch = data["SO_switch"]
    SpinPol = data["SpinPol"]
    xc_type = data["xc_type"]
    time_rev = data["time_rev"]
    E_Temp = data["E_Temp"]
    kmesh = data["kmesh"]
    SCF_criterion = data["SCF_criterion"]
    pao_file = data["pao_file"]
    pspot_file = data["pspot_file"]
    OLP = data["OLP"]
    Hks = data["Hks"]
    iHks = data["iHks"]
    DM = data["DM"]
    iDM = data["iDM"]
    ChemP = data["ChemP"]
    Eele = data["Eele"]
    Etot = data["Etot"]
    ForceAll = data["ForceAll"]
    scf_inputfile = data["scf_inputfile"]

    close(data)

    
    return LCPAO_model(
        Natom, Nspecies, Nspin, atom2spe, 
        Atoms_symbol, Atoms_Cut1, Atoms_pao, Atoms_Core_Charge, 
        Init_Atoms_Nspin, Init_Atoms_Angle, Atoms_Angle, Total_SpinS,
        Latvecs, Recvecs, Gxyz, TCpyCell, atv, atv_ijk, FNAN, natn, ncn, Total_NumOrbs, MP, 
        Grid_Origin, Ngrid, SO_switch, SpinPol, xc_type, time_rev, E_Temp, kmesh,
        SCF_criterion, pao_file, pspot_file, 
        OLP, Hks, iHks, DM, iDM, ChemP, Eele, Etot, ForceAll, scf_inputfile
    )
end
