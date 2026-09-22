function Print_Cube(
    filepath::String, 
    Ecut::AbstractFloat, 
    _mode::Union{Vector{String}, String}; 
    kpts=nothing,
    Nk_hwf=nothing,
    pflag=nothing)

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    MKL.set_num_threads(1)

    if nprocs > 1
        error("please run serial.")
    end


    filename, ext = splitext(basename(filepath))
    if ext ≠ ".jld2"
        error("not support filepath.")
    end


    mode_type = typeof(_mode)
    if mode_type == String
        mode = [_mode]
    elseif mode_type == Vector{String}
        mode = _mode
    else
        error("please check mode.")
    end


    Nmode = length(mode)
    mode = lowercase.(mode)
    
    for i = 1:Nmode
        if mode[i] ∉ ("rho", "psi", "hwf")
            error("please check mode.")
        end

        if mode[i] == "psi" && isnothing(kpts)
            error("please set kpoints for psi")
        end

        if mode[i] == "hwf" && isnothing(kpts)
            error("please set kpoints for hwf")
        end

        if mode[i] == "hwf" && isnothing(Nk_hwf)
            error("please set kpoints for hwf")
        end

        if mode[i] == "hwf" && isnothing(pflag)
            error("please set kpoints for hwf")
        end
    end


    if !isnothing(pflag)
        if length(pflag) ≠ 3
            error("please check pflag.")
        end

        if sum(pflag) > 1
            error("please check pflag.")
        end
    end


    if !isnothing(kpts)
        if typeof(kpts) == Vector{Float64}
            kpts2 = [kpts]
        elseif typeof(kpts) == Vector{Vector{Float64}}
            kpts2 = kpts
        else
            error("please check kpts")
        end
    else
        kpts2 = nothing
    end





    material = Load_LCPAODFT_model(filepath)
    Nspin = material.Nspin
    TCpyCell = material.TCpyCell
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Atoms_symbol = material.Atoms_symbol
    Atoms_Cut1 = material.Atoms_Cut1
    Atoms_pao = material.Atoms_pao
    Init_Atoms_Nspin = material.Init_Atoms_Nspin
    Init_Atoms_Angle = material.Init_Atoms_Angle
    Latvecs = material.Latvecs
    FNAN = material.FNAN
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    Grid_Origin = material.Grid_Origin
    SpinPol = material.SpinPol
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    DM = material.DM
    Ngrid = Calc_Ngrid(Ecut, Latvecs)


    if !isnothing(pflag)
        system = Check_system(FNAN, ncn, atv_ijk)
        if system == "bluk"
            # error("not support Print_Cube(HWF) for bulk")
        end
    end

    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)

    pao = Vector{PAO}(undef, Nspecies)
    pspot = Vector{Pspot}(undef, Nspecies)
    for spe = 1:Nspecies
        pspot[spe] = Read_VPS(Spe_symbol[spe], Spe_extra[spe], xc_type, SO_switch)
        pao[spe] = Read_PAO(pspot[spe].Spe_Core_Charge, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end



    ucell = UCell(Nspin, TCpyCell, Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin, Total_NumOrbs)
    system_grid = ucell.system_grid


    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)
    ADensity_Grid, _, Density_Grid = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, pao, pspot, ucell)



    for i = 1:Nmode
        if mode[i] == "rho"
            DM_1D = Set_DM_Vec2DM(DM, system_grid)
            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            else
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM_1D, Density_Grid)
                diagonalize_nc_density!(Density_Grid)
            end
            
            Print_Density(filename, SpinPol, Atoms_symbol, system_grid, ADensity_Grid, Density_Grid)
        elseif mode[i] == "psi"
            Enk, Bulk_HOMO, HOMOs_Coef = Calc_psi_HOMO(material, kpts2)
            Print_psi(filename, kpts2, material, ucell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)
        elseif mode[i] == "hwf"
            OLPpos = Set_OLPpos(Orbs_Grid, ucell)
            _, Umjk = Generate_HWF(material, OLPpos, Nk_hwf, pflag)
            HOMOs_Coef = Calc_hwfs_HOMO(material, Umjk)
            Print_hwfs(filename, zeros(3), material, ucell, Orbs_Grid, HOMOs_Coef)
        end
    end

end


export Print_Cube_HWF
function Print_Cube_HWF(
    filepath::String, 
    Ecut::AbstractFloat,
    _mode::Union{Vector{String}, String}; 
    kpts=nothing,
    Nk_hwf=nothing,
    pflag=nothing)

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    if nprocs > 1
        # error("please run serial.")
    end


    filename, ext = splitext(basename(filepath))
    if ext ≠ ".jld2"
        error("not support filepath.")
    end


    mode_type = typeof(_mode)
    if mode_type == String
        mode = [_mode]
    elseif mode_type == Vector{String}
        mode = _mode
    else
        error("please check mode.")
    end


    Nmode = length(mode)
    mode = lowercase.(mode)
    
    #=
    for i = 1:Nmode
        if mode[i] ∉ ("rho", "psi", "hwf")
            error("please check mode.")
        end

        if mode[i] == "psi" && isnothing(kpts)
            error("please set kpoints for psi")
        end

        if mode[i] == "hwf" && isnothing(kpts)
            error("please set kpoints for hwf")
        end

        if mode[i] == "hwf" && isnothing(Nk_hwf)
            error("please set kpoints for hwf")
        end

        if mode[i] == "hwf" && isnothing(pflag)
            error("please set kpoints for hwf")
        end
    end
    =#


    if !isnothing(pflag)
        if length(pflag) ≠ 3
            error("please check pflag.")
        end

        if sum(pflag) > 1
            error("please check pflag.")
        end
    end


    if !isnothing(kpts)
        if typeof(kpts) == Vector{Float64}
            kpts2 = [kpts]
        elseif typeof(kpts) == Vector{Vector{Float64}}
            kpts2 = kpts
        else
            error("please check kpts")
        end
    else
        kpts2 = nothing
    end





    material = Load_LCPAODFT_model(filepath)
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Atoms_symbol = material.Atoms_symbol
    Atoms_Cut1 = material.Atoms_Cut1
    Atoms_pao = material.Atoms_pao
    Init_Atoms_Nspin = material.Init_Atoms_Nspin
    Init_Atoms_Angle = material.Init_Atoms_Angle
    Latvecs = material.Latvecs
    FNAN = material.FNAN
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    Grid_Origin = material.Grid_Origin
    SpinPol = material.SpinPol
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    DM = material.DM
    Ngrid = Calc_Ngrid(Ecut, Latvecs)


    if !isnothing(pflag)
        system = Check_system(FNAN, ncn, atv_ijk)
        if system == "bluk"
            # error("not support Print_Cube(HWF) for bulk")
        end
    end

    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)

    pao = Vector{PAO}(undef, Nspecies)
    pspot = Vector{Pspot}(undef, Nspecies)
    for spe = 1:Nspecies
        pspot[spe] = Read_VPS(Spe_symbol[spe], Spe_extra[spe], xc_type, SO_switch)
        pao[spe] = Read_PAO(pspot[spe].Spe_Core_Charge, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end



    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs)
    system_grid = ucell.system_grid


    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)



    for i = 1:Nmode
        if mode[i] == "rho"
            DM_1D = Set_DM_Vec2DM(DM, system_grid)

            ADensity_Grid, _, Density_Grid = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, pao, pspot, ucell)
            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            else
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM_1D, Density_Grid)
                diagonalize_nc_density!(Density_Grid)
            end
            
            Print_Density(filename, SpinPol, Atoms_symbol, system_grid, ADensity_Grid, Density_Grid)
        elseif mode[i] == "psi"
            Enk, Bulk_HOMO, HOMOs_Coef = Calc_psi_HOMO(material, kpts2)
            Print_psi(filename, kpts2, material, ucell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)
        elseif mode[i] == "hwf"
            kpoints = KPoints((11,11,1), false, 0.0)
            OLPpos = Set_OLPpos(Orbs_Grid, ucell)
            OLPpos2 = Set_OLPpos2(Orbs_Grid, ucell)
            Cnk, Wannier_Center, Umjk = Generate_HWF_slab(material, OLPpos, kpoints, 3)
            # Cnk, Wannier_Center, Umjk = MPI_Generate_HWF_slab(material, OLPpos, kpoints, 3)
            Calc_Spread_HWF_NonCol(material, kpoints, Wannier_Center, OLPpos2, Cnk, Umjk, 3)
            # HOMOs_Coef = Calc_hwfs_HOMO(material, Cnk, Umjk)
            # Print_hwfs(filename, [zeros(3)], material, ucell, Orbs_Grid, HOMOs_Coef)
        end
    end

    @show LCPAODFT.timer
end


export temp_Print_Cube
function temp_Print_Cube(
    filepath::String, 
    Ecut,
    _mode::Union{Vector{String}, String};
    Ngrid = nothing,
    kpts = nothing)

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    LCPAODFT.reset_timer!(LCPAODFT.timer)

    if nprocs > 1
        error("please run serial.")
    end


    filename, ext = splitext(basename(filepath))
    if ext ≠ ".jld2"
        error("not support filepath.")
    end


    mode_type = typeof(_mode)
    if mode_type == String
        mode = [_mode]
    elseif mode_type == Vector{String}
        mode = _mode
    else
        error("please check mode.")
    end


    Nmode = length(mode)
    mode = lowercase.(mode)
    
    for i = 1:Nmode
        if mode[i] ∉ ("rho", "psi", "psi_1d")
            error("please check mode.")
        end

        if mode[i] == "psi" && isnothing(kpts)
            error("please set kpoints for psi")
        end
    end


    if !isnothing(kpts)
        if typeof(kpts) == Vector{Float64}
            kpts2 = [kpts]
        elseif typeof(kpts) == Vector{Vector{Float64}}
            kpts2 = kpts
        else
            error("please check kpts")
        end
    else
        kpts2 = nothing
    end
    @show kpts2


    if !isnothing(Ngrid)
        if length(Ngrid) ≠ 3
            error("please check Ngrid")
        end

        for i = 1:3
            if Ngrid[i] <= 0
                error("please check Ngrid")
            end
        end
    end





    material = Load_LCPAODFT_model(filepath)
    Print_LCPAO_model(filepath, material)
    
    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Atoms_symbol = material.Atoms_symbol
    Atoms_Cut1 = material.Atoms_Cut1
    Atoms_pao = material.Atoms_pao
    Init_Atoms_Nspin = material.Init_Atoms_Nspin
    Init_Atoms_Angle = material.Init_Atoms_Angle
    Latvecs = material.Latvecs
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    Grid_Origin = material.Grid_Origin
    SpinPol = material.SpinPol
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    DM = material.DM

    if isnothing(Ngrid)
        Ngrid = Calc_Ngrid(Ecut, Latvecs)
    end
    @show Ngrid



    Spe_symbol, Spe_cutoff, Spe_orb, Spe_extra = Get_Atoms_data(Atoms_pao)

    pao = Vector{PAO}(undef, Nspecies)
    pspot = Vector{Pspot}(undef, Nspecies)
    for spe = 1:Nspecies
        pspot[spe] = Read_VPS(Spe_symbol[spe], Spe_extra[spe], xc_type, SO_switch)
        pao[spe] = Read_PAO(pspot[spe].Spe_Core_Charge, Spe_symbol[spe], Spe_cutoff[spe], Spe_orb[spe], Spe_extra[spe])
    end



    ucell = UCell(Latvecs, Natom, atom2spe, Gxyz, Atoms_Cut1, Ngrid, Grid_Origin; Total_NumOrbs)
    system_grid = ucell.system_grid


    Orbs_Grid = Set_Orbitals_Grid(pao, ucell)
    ADensity_Grid, _, Density_Grid = Set_AdenPCC_Grid(SpinPol, Init_Atoms_Nspin, Init_Atoms_Angle, pao, pspot, ucell)



    for i = 1:Nmode
        if mode[i] == "rho"
            DM_1D = Set_DM_Vec2DM(DM, system_grid)
            if SpinPol == "off"
                Set_Density_Grid_nonpol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            elseif SpinPol == "on"
                Set_Density_Grid_pol!(ucell, Orbs_Grid, DM_1D, Density_Grid)
            else
                Set_Density_Grid_nc!(ucell, Orbs_Grid, DM_1D, Density_Grid)
                diagonalize_nc_density!(Density_Grid)
            end
            
            Print_Density(filename, SpinPol, Atoms_symbol, system_grid, ADensity_Grid, Density_Grid)
        elseif mode[i] == "psi"
            Enk, Bulk_HOMO, HOMOs_Coef = Calc_psi_HOMO(material, kpts2)
            Print_psi(filename, kpts2, material, ucell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)
        elseif mode[i] == "psi_1d"
            Enk, Bulk_HOMO, HOMOs_Coef = Calc_psi_HOMO(material, kpts2)
            Print_psi_1D(filename, kpts2, material, ucell, Orbs_Grid, Enk, Bulk_HOMO, HOMOs_Coef)
        end
    end

end
