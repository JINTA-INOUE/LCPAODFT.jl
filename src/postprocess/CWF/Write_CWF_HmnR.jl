function Write_CWF_HmnR(cwf_setup::Union{CWF_Setup, CWF_Setup_MO}, DMfunc, NCell, cell_list, cell_list_ijk, HmnR)

    material = cwf_setup.material
    filename = cwf_setup.filename
    spinsize = cwf_setup.spinsize
    gsize = cwf_setup.gsize
    Ngsize = cwf_setup.Ngsize
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    kmesh = cwf_setup.kmesh
    SO_switch = material.SO_switch
    Dis_Energy = cwf_setup.Dis_Energy
    ChemP = material.ChemP
    SpinPol = material.SpinPol
    weight_type = cwf_setup.weight_type
    scf_inputfile = material.scf_inputfile
    
    # data = open(pwd()*"/"*PROGRAM_FILE, "r")
    # cwf_inputfile = readlines(data)
    # close(data)

    cwf_inputfile = [""]


    println("Write $(filename).CWF.jld2")
    jldopen("$(filename).CWF.jld2", "w") do file
        file["Dates"] = now()
        file["spinsize"] = spinsize
        file["gsize"] = gsize
        file["Ngsize"] = Ngsize
        file["Latvecs"] = Latvecs
        file["Recvecs"] = Recvecs
        file["SpinPol"] = SpinPol
        file["SO_switch"] = SO_switch
        file["kmesh"] = kmesh
        file["Dis_Energy"] = Dis_Energy
        file["DMfunc"] = DMfunc
        file["NCell"] = NCell
        file["cell_list"] = cell_list
        file["cell_list_ijk"] = cell_list_ijk
        file["HmnR"] = HmnR
        file["ChemP"] = ChemP
        file["weight_type"] = weight_type
        file["scf_inputfile"] = scf_inputfile
        file["cwf_inputfile"] = cwf_inputfile
    end
end