function memory_usage(cal_force::Bool, SpinPol::String, pao::Vector{PAO}, pspot::Vector{Pspot}, ucell::UCell, Orbs_Grid, electron, dft_mixing, xc_func, Ham::Hamiltonian)
    
    ram_pao = memory_usage(pao)
    ram_pspot = memory_usage(pspot)
    ram_ucell = memory_usage(ucell)
    ram_dft_mixing = memory_usage(dft_mixing)
    ram_xc_func = memory_usage(xc_func)
    ram_Ham = memory_usage_Hamiltonian(Ham)
    ram_Orbs_Grid = memory_usage_Orbs_Grid(cal_force, Orbs_Grid)
    ram_electron = memory_usage(electron)
    ram_Den_Pot = memory_usage_Density_Potentials(ucell.system_grid.Ngrid, SpinPol)
    ram_Hks_DM = memory_usage_Hks_DM(SpinPol, ucell.system_grid.MPI_Hsize, ucell.system_grid.Total_Hsize)
    ram_EDM = memory_usage_EDM(SpinPol, cal_force, ucell.system_grid.Total_Hsize)

    ram_total = ram_pao + ram_pspot + ram_ucell + ram_dft_mixing + ram_xc_func + ram_Ham + ram_Orbs_Grid + ram_Den_Pot + ram_Hks_DM + ram_EDM

    println("<memory_usage>")
    @printf("\tram_pao : %6.3f MB\n", ram_pao)
    @printf("\tram_pspot : %6.3f MB\n", ram_pspot)
    @printf("\tram_ucell : %6.3f MB\n", ram_ucell)
    @printf("\tram_dft_mixing : %6.3f MB\n", ram_dft_mixing)
    @printf("\tram_xc_func : %6.3f MB\n", ram_xc_func)
    @printf("\tram_Ham : %6.3f MB\n", ram_Ham)
    @printf("\tram_Orbs_Grid : %6.3f MB\n", ram_Orbs_Grid)
    @printf("\tram_electron : %6.3f MB\n", ram_electron)
    @printf("\tram_Den_Pot : %6.3f MB\n", ram_Den_Pot)
    @printf("\tram_Hks_DM : %6.3f MB\n", ram_Hks_DM)
    @printf("\tram_EDM : %6.3f MB\n", ram_EDM)
    @printf("\tram_total : %6.3f MB\n", ram_total)
end


function memory_usage(pao::Vector{PAO})

    Nspecies = length(pao)

    MB = 1024^2
    ram = 0.0
    for spe = 1:Nspecies
        ram += Base.summarysize(pao[spe].Spe_PAO_XV)
        ram += Base.summarysize(pao[spe].Spe_PAO_RV)
        ram += Base.summarysize(pao[spe].Spe_Atomic_Den)
        ram += Base.summarysize(pao[spe].Spe_PAO_RWF)
        ram += Base.summarysize(pao[spe].Spe_RF_Bessel)
    end

    return ram/MB
end


function memory_usage(pspot::Vector{Pspot})

    Nspecies = length(pspot)

    MB = 1024^2
    ram = 0.0
    for spe = 1:Nspecies
        ram += Base.summarysize(pspot[spe].Spe_VPS_XV)
        ram += Base.summarysize(pspot[spe].Spe_VPS_RV)
        ram += Base.summarysize(pspot[spe].Spe_Vcore)
        ram += Base.summarysize(pspot[spe].Spe_Atomic_PCC)
        ram += Base.summarysize(pspot[spe].Spe_VNL)
    end

    return ram/MB
end


function memory_usage(ucell::UCell)

    MB = 1024^2
    system_grid = ucell.system_grid
    ram = Base.summarysize(system_grid.atv)
    ram += Base.summarysize(system_grid.Dis)
    ram += Base.summarysize(ucell.GridListAtom)
    ram += Base.summarysize(ucell.CellListAtom)
    ram += Base.summarysize(ucell.MPI_GListTAtoms1)
    ram += Base.summarysize(ucell.MPI_GListTAtoms2)
    ram += Base.summarysize(ucell.density_scratch)
    ram += Base.summarysize(ucell.density_matrix_scratch)
    ram += Base.summarysize(ucell.density_orbital_scratch)
    ram += Base.summarysize(ucell.hamiltonian_orbital_scratch)
    ram += Base.summarysize(ucell.hamiltonian_product_scratch)

    return ram/MB
end


function memory_usage_Hamiltonian(Ham::Hamiltonian)

    MB = 1024^2
    ram = Base.summarysize(Ham.OLP)
    ram += Base.summarysize(Ham.MPI_Hkin)
    ram += Base.summarysize(Ham.MPI_HNL)
    ram += Base.summarysize(Ham.MPI_iHNL)
    ram += Base.summarysize(Ham.MPI_HVNA)
    ram += Base.summarysize(Ham.MPI_NLPforce)
    ram += Base.summarysize(Ham.MPI_DS_VNAforce)
    ram += Base.summarysize(Ham.MPI_HVNA2force)
    ram += Base.summarysize(Ham.MPI_HVNA3force)

    return ram/MB
end


function memory_usage(electron::AbstractBloch)

    MB = 1024^2
    ram = Base.summarysize(electron.FF)
    ram += Base.summarysize(electron.Enk)
    ram += Base.summarysize(electron.Cnk)

    return ram/MB
end


function memory_usage(dft_mixing::RMM_DIISH_Mixing)

    MB = 1024^2
    ram = Base.summarysize(dft_mixing.Density_xyz)
    ram += Base.summarysize(dft_mixing.HisH)
    ram += Base.summarysize(dft_mixing.ResH)

    return ram/MB
end


function memory_usage(xc_func::XC_GGA_PBE)
    MB = 1024^2
    ram = Base.summarysize(xc_func.Diff_Coef)
    ram += Base.summarysize(xc_func.Vxc_Grid)
    ram += Base.summarysize(xc_func.dDensity_Grid)
    ram += Base.summarysize(xc_func.dEXC_dGD)

    return ram/MB
end


function memory_usage_Orbs_Grid(cal_force, Orbs_Grid)
    MB = 1024^2
    ram = Base.summarysize(Orbs_Grid)
    if cal_force
        ram += 3*Base.summarysize(Orbs_Grid)
    end

    return ram/MB
end


function memory_usage_Density_Potentials(Ngrid, SpinPol)

    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end

    MB = 1024^2
    ram = (Nspin+1)*prod(Ngrid)     # dVHart_Grid, Vpot_Grid
    ram += (Nspin+2)*prod(Ngrid)    # ADensity_Grid, PCCDensity_Grid, Density_Grid

    return sizeof(Float64)*ram/MB
end


function memory_usage_Hks_DM(SpinPol, MPI_Hsize, Total_Hsize)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end

    MB = 1024^2
    ram = 2*Nspin*Total_Hsize     # DM, Hks
    if SpinPol == "nc"
        ram += 5*Total_Hsize    # iDM, iHks
    end
    ram += Nspin*MPI_Hsize[myrank+1]    # MPI_Hks


    return sizeof(Float64)*ram/MB
end


function memory_usage_EDM(SpinPol, cal_force, Total_Hsize)

    Nspin_EDM = ifelse(SpinPol=="off", 1, 2)
    EDMsize = ifelse(cal_force, Total_Hsize, 1)
    MB = 1024^2
    ram = Nspin_EDM*EDMsize*sizeof(Float64)

    return ram/MB
end
