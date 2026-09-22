@timeit timer "Set_MLWF_ExpnCoef" function Set_MLWF_ExpnCoef(mlwf_setup::MLWF_Setup, kpoints::KPoints, chk::Vector{w90_Chk}, BANDNUM, Nk, Cnk)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    filename = mlwf_setup.filename
    WANNUM = mlwf_setup.WANNUM
    MLWF_Plot_SuperCells = mlwf_setup.MLWF_Plot_SuperCells
    Plot_NCell = prod(2*MLWF_Plot_SuperCells.+1)
    write_coef = mlwf_setup.write_coef

    AllNkpt = kpoints.AllNkpt
    

    Umnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umnk[spin] = Vector{Matrix{ComplexF64}}(undef, AllNkpt)
        for ik = 1:AllNkpt
            Umnk[spin][ik] = zeros(ComplexF64, BANDNUM, WANNUM)
        end
    end

    for spin = 1:spinsize
        Umnk[spin] = gauge_matrices(chk[spin])
    end

    if SpinPol ∈ ("off", "on")
        coeffsize = ifelse(SpinPol=="off", 1, 2)
        MLWF_ExpnCoef = Vector{Vector{Vector{Vector{Float64}}}}(undef, coeffsize)
        for spin = 1:coeffsize
            MLWF_ExpnCoef[spin] = Vector{Vector{Vector{Float64}}}(undef, WANNUM)
            for pst = 1:WANNUM
                MLWF_ExpnCoef[spin][pst] = Vector{Vector{Float64}}(undef, Plot_NCell)
                for cell = 1:Plot_NCell
                    MLWF_ExpnCoef[spin][pst][cell] = zeros(Float64, Nfsize)
                end
            end
        end
    elseif SpinPol == "nc"
        MLWF_ExpnCoef = Vector{Vector{Vector{ComplexF64}}}(undef, WANNUM)
        for pst = 1:WANNUM
            MLWF_ExpnCoef[pst] = Vector{Vector{ComplexF64}}(undef, Plot_NCell)
            for cell = 1:Plot_NCell
                MLWF_ExpnCoef[pst][cell] = zeros(ComplexF64, Nfsize)
            end
        end
    end    


    if SpinPol ∈ ("off", "on")
        Set_MLWF_ExpnCoef_Col!(mlwf_setup, kpoints, BANDNUM, Nk, Cnk, Umnk, MLWF_ExpnCoef)
    else
        Set_MLWF_ExpnCoef_NonCol!(mlwf_setup, kpoints, BANDNUM, Nk, Cnk, Umnk, MLWF_ExpnCoef)
    end


    if write_coef && myrank == 0
        println("Write $filename.ExpnCoef.jld2")
        jldopen("$filename.ExpnCoef.jld2", "w") do file
            file["Dates"] = now()
            file["SpinPol"] = SpinPol
            file["spinsize"] = spinsize
            file["Nfsize"] = Nfsize
            file["WANNUM"] = WANNUM
            file["WF_Plot_SuperCells"] = MLWF_Plot_SuperCells
            file["WF_ExpnCoef"] = MLWF_ExpnCoef
        end
    end


    return MLWF_ExpnCoef
end


function Set_MLWF_ExpnCoef_Col!(mlwf_setup::MLWF_Setup, kpoints::KPoints, BANDNUM, Nk, Cnk, Umnk, MLWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nfsize = sum(material.Total_NumOrbs)
    WANNUM = mlwf_setup.WANNUM
    MLWF_Plot_SuperCells = mlwf_setup.MLWF_Plot_SuperCells
    Plot_NCell = prod(2*MLWF_Plot_SuperCells.+1)
    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    Plot_cell_ijk = CWF_plot_cells(MLWF_Plot_SuperCells)

    phases = zeros(ComplexF64, Plot_NCell)
    band_coef = zeros(ComplexF64, Nfsize, WANNUM)

    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        @inbounds for cell = 1:Plot_NCell
            l, m, n = Plot_cell_ijk[cell]
            phases[cell] = cispi(2*(ka*l + kb*m + kc*n))/AllNkpt
        end

        MinN = Nk[spin,ik,1]
        _Cnk = Cnk[spin][ik]
        _Umnk = Umnk[spin][ik]

        @inbounds for proj = 1:WANNUM, ist = 1:Nfsize
            temp = ComplexF64(0.0, 0.0)
            for mu = 1:BANDNUM
                temp += _Umnk[mu,proj]*_Cnk[ist,mu+MinN]
            end
            band_coef[ist,proj] = temp
        end

        @inbounds for proj = 1:WANNUM, cell = 1:Plot_NCell, ist = 1:Nfsize
            MLWF_ExpnCoef[spin][proj][cell][ist] += real(band_coef[ist,proj]*phases[cell])
        end
    end

    for spin = 1:spinsize, proj = 1:WANNUM, cell = 1:Plot_NCell
        MPI.Allreduce!(MLWF_ExpnCoef[spin][proj][cell], MPI.SUM, comm)
    end
end


function Set_MLWF_ExpnCoef_NonCol!(mlwf_setup::MLWF_Setup, kpoints::KPoints, BANDNUM, Nk, Cnk, Umnk, MLWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = mlwf_setup.material
    Nfsize = 2*sum(material.Total_NumOrbs)
    WANNUM = mlwf_setup.WANNUM
    MLWF_Plot_SuperCells = mlwf_setup.MLWF_Plot_SuperCells
    Plot_NCell = prod(2*MLWF_Plot_SuperCells.+1)
    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    Plot_cell_ijk = CWF_plot_cells(MLWF_Plot_SuperCells)

    phases = zeros(ComplexF64, Plot_NCell)
    band_coef = zeros(ComplexF64, Nfsize, WANNUM)


    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        @inbounds for cell = 1:Plot_NCell
            l, m, n = Plot_cell_ijk[cell]
            phases[cell] = cispi(2*(ka*l + kb*m + kc*n))/AllNkpt
        end
        MinN = Nk[1,ik,1]
        _Cnk = Cnk[1][ik]
        _Umnk = Umnk[1][ik]
        @inbounds for proj = 1:WANNUM, ist = 1:Nfsize
            temp = ComplexF64(0.0, 0.0)
            @inbounds for mu = 1:BANDNUM
                temp += _Umnk[mu,proj]*_Cnk[ist,mu+MinN]
            end
            band_coef[ist,proj] = temp
        end
        @inbounds for proj = 1:WANNUM, cell = 1:Plot_NCell, ist = 1:Nfsize
            MLWF_ExpnCoef[proj][cell][ist] += band_coef[ist,proj]*phases[cell]
        end
    end


    for proj = 1:WANNUM, cell = 1:Plot_NCell
        MPI.Allreduce!(MLWF_ExpnCoef[proj][cell], MPI.SUM, comm)
    end
end