function Write_CIFfile(filename::String, Natom, Latvecs, Gxyz_frac, Atoms_symbol)

    if !(Natom == length(Gxyz_frac) == length(Atoms_symbol))
        error("please check input")
    end

    Cell_Volume = abs(det(Latvecs))

    length_a = norm(Latvecs[1,:])
    length_b = norm(Latvecs[2,:])
    length_c = norm(Latvecs[3,:])
    
    t1 = dot(Latvecs[2,:], Latvecs[3,:])
    alpha = ifelse(abs(t1) < 1e-14, 90.0, acos(t1/(length_b*length_c)))/pi*180
    
    t1 = dot(Latvecs[3,:], Latvecs[1,:])
    beta = ifelse(abs(t1) < 1e-14, 90.0, acos(t1/(length_c*length_a)))/pi*180
    
    t1 = dot(Latvecs[1,:], Latvecs[2,:])
    gamma = ifelse(abs(t1) < 1e-14, 90.0, acos(t1/(length_a*length_b)))/pi*180
    

    data = open("$filename.cif", "w")

    @printf(data, "data_%s\n", data)
    @printf(data, "_audit_creation_date              %s\n", today())
    @printf(data, "_audit_creation_method            'LCPAODFT.jl'\n")

    @printf(data, "_symmetry_space_group_name_H-M    'P1'\n")
    @printf(data, "_symmetry_Int_Tables_number       1\n")
    @printf(data, "_symmetry_cell_setting            triclinic\n")

    @printf(data, "loop_\n")
    @printf(data, "_symmetry_equiv_pos_as_xyz\n")
    @printf(data, "  x,y,z\n")

    @printf(data, "_cell_length_a%26.8f\n", length_a/Ang_to_bohr)
    @printf(data, "_cell_length_b%26.8f\n", length_b/Ang_to_bohr)
    @printf(data, "_cell_length_c%26.8f\n", length_c/Ang_to_bohr)

    @printf(data, "_cell_angle_alpha%24.8f\n", alpha)
    @printf(data, "_cell_angle_beta %24.8f\n", beta)
    @printf(data, "_cell_angle_gamma%24.8f\n", gamma)

    @printf(data, "loop_\n")
    @printf(data, "_atom_site_label\n")
    @printf(data, "_atom_site_type_symbol\n")
    @printf(data, "_atom_site_fract_x\n")
    @printf(data, "_atom_site_fract_y\n")
    @printf(data, "_atom_site_fract_z\n")
    @printf(data, "_atom_site_Uiso_or_equiv\n")
    @printf(data, "_atom_site_adp_type\n")
    @printf(data, "_atom_site_occupancy\n")


    for atom = 1:Natom
        @printf(data, "%s%-6d%-3s%10.8f %10.8f %10.8f %10.8f  Uiso   1.00\n",
                Atoms_symbol[atom],
                atom, 
	            Atoms_symbol[atom],
	            Gxyz_frac[atom][1], Gxyz_frac[atom][2], Gxyz_frac[atom][3],
                0.0)
    end
    close(data)
end


function Write_xyzfile(filename::String, Natom, Gxyz, Atoms_symbol)

    if !(Natom == length(Gxyz) == length(Atoms_symbol))
        error("please check input")
    end

    data = open("$filename.xyz", "w")
    for atom = 1:Natom
        @printf(data, "%4s  %18.14f %18.14f %18.14f\n", Atoms_symbol[atom], Gxyz[atom][1]/Ang_to_bohr, Gxyz[atom][2]/Ang_to_bohr, Gxyz[atom][3]/Ang_to_bohr)
    end
    close(data)
end


function Write_xyzfile(filename::String, Natom, Gxyz, ForceAll, Atoms_symbol)

    if !(Natom == length(Gxyz) == length(Atoms_symbol))
        error("please check input")
    end

    data = open("$filename.xyz", "w")
    for atom = 1:Natom
        @printf(data, "%4s  %8.8f  %8.8f  %8.8f  %18.15f %18.15f %18.15f\n", 
        Atoms_symbol[atom], 
        Gxyz[atom][1]/Ang_to_bohr, Gxyz[atom][2]/Ang_to_bohr, Gxyz[atom][3]/Ang_to_bohr,
        ForceAll[atom,1], ForceAll[atom,2], ForceAll[atom,3])
    end
    close(data)
end


#=
function Write_outFile()

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    nthreads = Threads.nthreads()

    _inputfile = @__FILE__
    data = open(_input_file, "r")
    inputfile = readlines(data)
    close(data)

    File = open(filename*".out", "w")
    println(File, "===========================================================")
    println(File, "")
    println(File, "  This calculation was performed by LCPAODFT Ver. 1.0.0")
    println(File, "  using $nprocs MPI processes and $nthreads OpenMP threads.")
    println(File, "")
    println(File, "  $(now())")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    for i = 1:length(inputfile)
        println(File, inputfile[i])
    end
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "  Required cutoff energy (Ryd) for 3D-grids = $Ecut")
    println(File, "  Num. of grids of a-, b-, and c-axes = $(Ngrid[1]) $(Ngrid[2]) $(Ngrid[3])")
    println(File, "")
    println(File, "  Cell_Volume = $Cell_Volume (Bohr^3)")
    println(File, "  GridVol     = $GridVol (Bohr^3)")
    println(File, "  Cell vectors (bohr) of the grid cell (gtv)")
    println(File, "    gtv_a = $(gLatvecs[1,1]), $(gLatvecs[1,2]), $(gLatvecs[1,3])")
    println(File, "    gtv_b = $(gLatvecs[2,1]), $(gLatvecs[2,2]), $(gLatvecs[2,3])")
    println(File, "    gtv_c = $(gLatvecs[3,1]), $(gLatvecs[3,2]), $(gLatvecs[3,3])")
    println(File, "   |gtv_a| = $(norm(gLatvecs[1,:]))")
    println(File, "   |gtv_b| = $(norm(gLatvecs[2,:]))")
    println(File, "   |gtv_c| = $(norm(gLatvecs[3,:]))")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       SCF history")
    println(File, "")
    for iter = 1:Conv_iter
        @printf(File, "    SCF = %3d  NormRD = %5.12f Uele = %5.12f\n", iter, NormRD[iter], HisEele[iter])
    end
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       Total Energy (Hartree, eV)")
    println(File, "")
    @printf(File, "  Uele.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Ukin.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  UH0.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  UH1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Una.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Unl.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc0.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Ucore.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Utot.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    @printf(File, "  Uxc1.       %5.12f     %5.12f\n", energy.Eele, energy.Eele*eV2Hartree)
    println(File, "")
    println(File, "  Chemical potential %5.12f     %5.12f", energy.ChemP, energy.ChemP*eV2Hartree)
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       Mulliken populations")
    println(File, "")
    println(File, "  Total spin moment (muB)")
    println(File, "")
    println(File, "")
    println(File, "")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       Dipole moment (Debye)")
    println(File, "")
    println(File, "  Absolute D")
    println(File, "             Dx                Dy                Dz")
    println(File, "  Total")
    println(File, "  Core")
    println(File, "  Electron")
    println(File, "  Back ground")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       xyz-coordinates (Ang) and forces (Hartree/Bohr)  ")
    println(File, "")
    println(File, "<coordinates.forces")
    println(File, "  $Natom")
    for atom = 1:Natom
        @printf(File, "    \n", atom, Atoms_symbol[atom], Gxyz[atom][1], Gxyz[atom][2], Gxyz[atom][3], Force[atom,1], Force[atom,2], Force[atom,3])
    end
    println(File, "coordinates.forces>")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       Fractional coordinates")
    println(File, "")
    println(File, "")
    println(File, "")
    println(File, "===========================================================")
    println(File, "")
    println(File, "                       Computational Time (second)")
    println(File, "")

    close(File)
end
=#


function WriteFile!(
    filename::String,
    SpinPol::String, Ngrid, ADensity_Grid, Density_Grid)

    # Density_Grid overwritten in utils/OutData.jl/WriteFile!
    if SpinPol == "off"
        @. Density_Grid[1] = Density_Grid[1] - ADensity_Grid
    else
        @. Density_Grid[1] = Density_Grid[1] - ADensity_Grid
        @. Density_Grid[2] = Density_Grid[2] - ADensity_Grid
    end
    
    println("Write $filename")
    Generate_rhoFile(filename, Ngrid, Density_Grid)
end


function WriteFile(
    filename::String,
    system_grid::System_Grid, 
    energy::Energy, force::Force,
    DM, Hks, iHks)

    Natom = system_grid.Natom
    CpyCell = system_grid.CpyCell
    TCpyCell = (2*CpyCell + 1)^3 - 1
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    Total_NumOrbs = system_grid.Total_NumOrbs
    atv_ijk = system_grid.atv_ijk
    Total_Hsize = system_grid.Total_Hsize
    ChemP = energy.ChemP
    Eele = energy.Eele
    Etot = energy.Etot
    ForceAll = force.ForceAll
    
    println("Write $filename")
    jldopen("$filename", "w") do file
        file["Natom"] = Natom
        file["TCpyCell"] = TCpyCell
        file["FNAN"] = FNAN
        file["natn"] = natn
        file["ncn"] = ncn
        file["Total_NumOrbs"] = Total_NumOrbs
        file["atv_ijk"] = atv_ijk
        file["Total_Hsize"] = Total_Hsize
        file["Hks"] = Hks
        file["iHks"] = iHks
        file["ChemP"] = ChemP
        file["Eele"] = Eele
        file["Etot"] = Etot
        file["ForceAll"] = ForceAll
        file["Dates"] = now()
    end
end


function WriteFile(
    mulliken_charge::Mulliken_Charge, 
    dft_setup::DFT_Setup, 
    system_grid::System_Grid, 
    dipole_moment,
    energy::Energy, force::Force,
    DM, iDM, OLP, Hks, iHks)

    filename = dft_setup.filename
    Nspin = dft_setup.Nspin
    Natom = dft_setup.Natom
    Nspecies = dft_setup.Nspecies
    atom2spe = dft_setup.atom2spe
    Latvecs = dft_setup.Latvecs
    Recvecs = dft_setup.Recvecs
    Init_Atoms_Nspin = dft_setup.Init_Atoms_Nspin
    Init_Atoms_Angle = dft_setup.Init_Atoms_Angle
    Atoms_symbol = dft_setup.Atoms_symbol
    Atoms_pao = dft_setup.Atoms_pao
    Atoms_Core_Charge = mulliken_charge.Atoms_Core_Charge
    Atoms_Angle = mulliken_charge.Angle_Spin
    Total_SpinS = mulliken_charge.Total_SpinS
    CpyCell = system_grid.CpyCell
    TCpyCell = (2*CpyCell + 1)^3 - 1
    Gxyz = system_grid.Gxyz
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    atv = system_grid.atv
    atv_ijk = system_grid.atv_ijk
    ncn = system_grid.ncn
    Grid_Origin = system_grid.Grid_Origin
    MP = system_grid.MP
    Total_NumOrbs = system_grid.Total_NumOrbs
    Atoms_Cut1 = system_grid.Atom_Cut1
    xc_type = dft_setup.xc_type
    SpinPol = dft_setup.SpinPol
    SO_switch = dft_setup.SO_switch
    xc_type = dft_setup.xc_type
    Ngrid = dft_setup.Ngrid
    kmesh = dft_setup.kmesh
    E_Temp = dft_setup.E_Temp
    SCF_criterion = dft_setup.SCF_criterion
    time_rev = dft_setup.time_rev
    pao_file = dft_setup.pao_file
    pspot_file = dft_setup.pspot_file
    ChemP = energy.ChemP
    Eele = energy.Eele
    Etot = energy.Etot
    ForceAll = force.ForceAll
    filename = dft_setup.filename

    data = open(pwd()*"/"*PROGRAM_FILE, "r")
    scf_inputfile = readlines(data)
    close(data)
    
    println("Write $filename.jld2 LCPAO_model")
    jldopen("$filename.jld2", "w") do file
        file["Natom"] = Natom
        file["Nspecies"] = Nspecies
        file["Nspin"] = Nspin
        file["atom2spe"] = atom2spe
        file["Atoms_symbol"] = Atoms_symbol
        file["Atoms_Cut1"] = Atoms_Cut1
        file["Atoms_pao"] = Atoms_pao
        file["Atoms_Core_Charge"] = Atoms_Core_Charge
        file["Init_Atoms_Nspin"] = Init_Atoms_Nspin
        file["Init_Atoms_Angle"] = Init_Atoms_Angle
        file["dipole_moment"] = dipole_moment
        file["Atoms_Angle"] = Atoms_Angle
        file["Total_SpinS"] = Total_SpinS
        file["Latvecs"] = Latvecs
        file["Recvecs"] = Recvecs
        file["Gxyz"] = Gxyz
        file["Grid_Origin"] = Grid_Origin
        file["TCpyCell"] = TCpyCell
        file["atv"] = atv
        file["atv_ijk"] = atv_ijk
        file["FNAN"] = FNAN
        file["natn"] = natn
        file["ncn"] = ncn
        file["Total_NumOrbs"] = Total_NumOrbs
        file["MP"] = MP
        file["Ngrid"] = Ngrid
        file["SO_switch"] = SO_switch
        file["SpinPol"] = SpinPol
        file["xc_type"] = xc_type
        file["time_rev"] = time_rev
        file["E_Temp"] = E_Temp
        file["kmesh"] = kmesh
        file["SCF_criterion"] = SCF_criterion
        file["pao_file"] = pao_file
        file["pspot_file"] = pspot_file
        file["OLP"] = OLP
        file["Hks"] = Hks
        file["iHks"] = iHks
        file["DM"] = DM
        file["iDM"] = iDM
        file["ChemP"] = ChemP
        file["Eele"] = Eele
        file["Etot"] = Etot
        file["ForceAll"] = ForceAll
        file["Dates"] = now()
        file["scf_inputfile"] = scf_inputfile
    end
end