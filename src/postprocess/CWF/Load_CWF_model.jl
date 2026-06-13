function Load_CWF_model(filepath::String)

    data = jldopen(filepath, "r")

    spinsize = data["spinsize"]
    gsize = data["gsize"]
    Ngsize = data["Ngsize"]
    Latvecs = data["Latvecs"]
    Recvecs = data["Recvecs"]
    SpinPol = data["SpinPol"]
    SO_switch = data["SO_switch"]
    kmesh = data["kmesh"]
    Dis_Energy = data["Dis_Energy"]
    DMfunc = data["DMfunc"]
    NCell = data["NCell"]
    cell_list = data["cell_list"]
    cell_list_ijk = data["cell_list_ijk"]
    HmnR = data["HmnR"]
    ChemP = data["ChemP"]
    weight_type = data["weight_type"]
    scf_inputfile = data["scf_inputfile"]
    cwf_inputfile = data["cwf_inputfile"]
    close(data)


    return CWF_model(
        spinsize, Latvecs, Recvecs, gsize, Ngsize, kmesh, 
        SpinPol, SO_switch, Dis_Energy, DMfunc,
        NCell, cell_list, cell_list_ijk,
        HmnR, ChemP, weight_type, 
        scf_inputfile, cwf_inputfile)
end