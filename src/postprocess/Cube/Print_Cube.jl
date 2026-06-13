function Print_Cube(
    filepath::String, 
    Ecut::AbstractFloat, 
    _mode::Union{Vector{String}, String}; 
    kpts=nothing)

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

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
        if mode[i] ∉ ("rho", "psi")
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
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    Grid_Origin = material.Grid_Origin
    SpinPol = material.SpinPol
    SO_switch = material.SO_switch
    xc_type = material.xc_type
    DM = material.DM
    Ngrid = Calc_Ngrid(Ecut, Latvecs)



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
        end
    end 
end