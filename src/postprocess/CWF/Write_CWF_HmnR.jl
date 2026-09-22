function _Write_CWF_metadata!(file, cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk; Ngsize=cwf_setup.Ngsize)

    material = cwf_setup.material
    spinsize = cwf_setup.spinsize
    gsize = cwf_setup.gsize
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    kmesh = cwf_setup.kmesh
    SO_switch = material.SO_switch
    Dis_Energy = cwf_setup.Dis_Energy
    ChemP = material.ChemP
    SpinPol = material.SpinPol
    weight_type = cwf_setup.weight_type
    scf_inputfile = material.scf_inputfile
    
    input_path = abspath(PROGRAM_FILE)
    cwf_inputfile = if isfile(input_path)
        open(input_path, "r") do data
            readlines(data)
        end
    else
        String[]
    end

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
    file["ChemP"] = ChemP
    file["weight_type"] = weight_type
    file["scf_inputfile"] = scf_inputfile
    file["cwf_inputfile"] = cwf_inputfile

    return nothing
end


function Write_CWF_HmnR(cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk, HmnR)

    filename = cwf_setup.filename

    println("\tWrite $(filename).CWF.jld2")
    jldopen("$(filename).CWF.jld2", "w") do file
        _Write_CWF_metadata!(file, cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk)
        file["HmnR"] = HmnR
    end
end


function Calc_Write_CWF_HmnR(
    cwf_setup,
    DMfunc,
    NCell,
    cell_list,
    cell_list_ijk,
    MinN,
    MaxN,
    Umnk,
    Enk,
    kpoints::KPoints;
    block_cells=_HmnR_block_cells(cwf_setup.Ngsize, NCell, kpoints.MPI_Nkpt),
)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    filename = cwf_setup.filename
    spinsize = cwf_setup.spinsize
    Ngsize = cwf_setup.Ngsize
    output_path = "$(filename).CWF.jld2"
    temporary_path = output_path * ".tmp"

    file = nothing
    HmnR_dataset = nothing

    if myrank == 0
        println("\tWrite $output_path (streaming)")
        # IOStream keeps the writer itself out of the memory map.  The
        # resulting contiguous dataset remains mmap-compatible when loaded.
        file = jldopen(temporary_path, "w"; iotype=IOStream)
        _Write_CWF_metadata!(file, cwf_setup, DMfunc, NCell, cell_list, cell_list_ijk)
        file["HmnR_storage"] = "contiguous-streamed-v1"
        HmnR_dataset = JLD2.create_dataset(
            file,
            "HmnR",
            ComplexF64,
            (Int(Ngsize), Int(Ngsize), Int(NCell), Int(spinsize));
            allocate=true,
        )
    end

    try
        Calc_HmnR_stream!(
            HmnR_dataset,
            spinsize,
            MinN,
            MaxN,
            NCell,
            Ngsize,
            cell_list_ijk,
            Umnk,
            Enk,
            kpoints;
            block_cells,
        )

        if myrank == 0
            close(file)
            file = nothing
            mv(temporary_path, output_path; force=true)
        end
    catch
        if myrank == 0 && file !== nothing
            close(file)
        end
        rethrow()
    end

    return nothing
end

