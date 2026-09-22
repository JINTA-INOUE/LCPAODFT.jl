function Load_MLWF_model(filepath::String)

    data = jldopen(filepath, "r")

    spinsize = data["spinsize"]
    Nwann = data["Nwann"]
    Latvecs = data["Latvecs"]
    Recvecs = data["Recvecs"]
    SpinPol = data["SpinPol"]
    SO_switch = data["SO_switch"]
    kmesh = data["kmesh"]
    Dis_Energy = data["Dis_Energy"]
    NCell = data["NCell"]
    cell_list_ijk = data["cell_list_ijk"]
    Rdegens = data["Rdegens"]
    HmnR = data["HmnR"]
    ChemP = data["ChemP"]
    scf_inputfile = data["scf_inputfile"]
    mlwf_inputfile = data["mlwf_inputfile"]
    close(data)


    return MLWF_model(
        spinsize, Latvecs, Recvecs, Nwann, kmesh, 
        SpinPol, SO_switch, Dis_Energy,
        NCell, cell_list_ijk,
        Rdegens, HmnR, ChemP, 
        scf_inputfile, mlwf_inputfile)
end