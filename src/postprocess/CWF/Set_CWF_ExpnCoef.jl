@timeit timer "Set_CWF_ExpnCoef" function Set_CWF_ExpnCoef(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, MinN, MaxN, Cnk, Umnk)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol=="off", 1, 2)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)

    filename = cwf_setup.filename
    spinsize = cwf_setup.spinsize
    gsize = cwf_setup.gsize
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    write_coef = cwf_setup.write_coef


    

    if SpinPol ∈ ("off", "on")
        CWF_ExpnCoef = Vector{Vector{Vector{Vector{Float64}}}}(undef, spinsize)
        for spin = 1:spinsize
            CWF_ExpnCoef[spin] = Vector{Vector{Vector{Float64}}}(undef, gsize)
            for pst = 1:gsize
                CWF_ExpnCoef[spin][pst] = Vector{Vector{Float64}}(undef, Plot_NCell)
                for cell = 1:Plot_NCell
                    CWF_ExpnCoef[spin][pst][cell] = zeros(Float64, Nfsize)
                end
            end
        end
    elseif SpinPol == "nc"
        CWF_ExpnCoef = Vector{Vector{Vector{ComplexF64}}}(undef, 2*gsize)
        for pst = 1:2*gsize
            CWF_ExpnCoef[pst] = Vector{Vector{ComplexF64}}(undef, Plot_NCell)
            for cell = 1:Plot_NCell
                CWF_ExpnCoef[pst][cell] = zeros(ComplexF64, Nfsize)
            end
        end
    end    


    if SpinPol ∈ ("off", "on")
        Set_CWF_ExpnCoef_Col!(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk, CWF_ExpnCoef)
    else
        Set_CWF_ExpnCoef_NonCol!(cwf_setup, kpoints, MinN, MaxN, Cnk, Umnk, CWF_ExpnCoef)
    end


    if write_coef && myrank == 0
        println("Write $(filename).ExpnCoef.jld2")
        jldopen("$(filename).ExpnCoef.jld2", "w") do file
            file["Dates"] = now()
            file["SpinPol"] = SpinPol
            file["spinsize"] = spinsize
            file["Nfsize"] = Nfsize
            file["gsize"] = gsize
            file["CWF_Plot_SuperCells"] = CWF_Plot_SuperCells
            file["CWF_ExpnCoef"] = CWF_ExpnCoef
        end
    end


    return CWF_ExpnCoef
end


function Set_CWF_ExpnCoef_Col!(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, MinN, MaxN, Cnk, Umnk, CWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="off", 1, 2)
    Nfsize = sum(material.Total_NumOrbs)
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    Ngsize = cwf_setup.Ngsize
    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    BANDNUM = MaxN - MinN + 1

    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        Plot_cell_ijk[cell] = zeros(Int32, 3)
    end

    cell = 0
    for l1 = -CWF_Plot_SuperCells[1]:CWF_Plot_SuperCells[1], l2 = -CWF_Plot_SuperCells[2]:CWF_Plot_SuperCells[2], l3 = -CWF_Plot_SuperCells[3]:CWF_Plot_SuperCells[3]
        cell += 1
        Plot_cell_ijk[cell] = [l1, l2, l3]
    end



    Umnk_tmp = zeros(ComplexF64, BANDNUM, Ngsize)
    Cnk_tmp = zeros(ComplexF64, Nfsize, Nfsize)

    for spin = 1:spinsize, proj = 1:Ngsize, cell = 1:Plot_NCell
        l, m, n = Plot_cell_ijk[cell]
        for ist = 1:Nfsize
            Sum = ComplexF64(0.0, 0.0)
            for ik = 1:MPI_Nkpt
                @. Cnk_tmp = Cnk[spin][ik]
                @. Umnk_tmp = Umnk[spin][ik]
                kRn = MPI_kpts[ik][1]*l + MPI_kpts[ik][2]*m + MPI_kpts[ik][3]*n
                ex = cispi(2*kRn)/AllNkpt
                temp = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:BANDNUM
                    temp += Umnk_tmp[μ,proj]*Cnk_tmp[ist,μ+MinN-1]
                end
                Sum += temp*ex
            end

            CWF_ExpnCoef[spin][proj][cell][ist] = real(Sum)
        end
    end

    for spin = 1:spinsize, proj = 1:Ngsize, cell = 1:Plot_NCell
        MPI.Allreduce!(CWF_ExpnCoef[spin][proj][cell], MPI.SUM, comm)
    end
end


function Set_CWF_ExpnCoef_NonCol!(cwf_setup::Union{CWF_Setup,CWF_Setup_MO}, kpoints::KPoints, MinN, MaxN, Cnk, Umnk, CWF_ExpnCoef)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    material = cwf_setup.material
    SpinPol = material.SpinPol
    Nfsize = 2*sum(material.Total_NumOrbs)
    CWF_Plot_SuperCells = cwf_setup.CWF_Plot_SuperCells
    Ngsize = cwf_setup.Ngsize
    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    BANDNUM = MaxN - MinN + 1

    Plot_NCell = prod(2*CWF_Plot_SuperCells.+1)
    Plot_cell_ijk = Vector{Vector{Int32}}(undef, Plot_NCell)
    for cell = 1:Plot_NCell
        Plot_cell_ijk[cell] = zeros(Int32, 3)
    end

    cell = 0
    for l1 = -CWF_Plot_SuperCells[1]:CWF_Plot_SuperCells[1], l2 = -CWF_Plot_SuperCells[2]:CWF_Plot_SuperCells[2], l3 = -CWF_Plot_SuperCells[3]:CWF_Plot_SuperCells[3]
        cell += 1
        Plot_cell_ijk[cell] = [l1, l2, l3]
    end

    Umnk_tmp = zeros(ComplexF64, Nfsize, Ngsize)
    Cnk_tmp = zeros(ComplexF64, Nfsize, Nfsize)

    for proj = 1:Ngsize, cell = 1:Plot_NCell
        l, m, n = Plot_cell_ijk[cell]
        for ist = 1:Nfsize
            Sum = ComplexF64(0.0, 0.0)
            for ik = 1:MPI_Nkpt
                @. Cnk_tmp = Cnk[1][ik]
                @. Umnk_tmp = Umnk[1][ik]
                kRn = MPI_kpts[ik][1]*l + MPI_kpts[ik][2]*m + MPI_kpts[ik][3]*n
                ex = cispi(2*kRn)/Nkpt
                temp = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:Nfsize
                    temp += Umnk_tmp[μ,proj]*Cnk_tmp[ist,μ+MinN-1]
                end
                Sum += temp*ex
            end

            CWF_ExpnCoef[proj][cell][ist] = Sum
        end
    end


    for proj = 1:Ngsize, cell = 1:Plot_NCell
        MPI.Allreduce!(CWF_ExpnCoef[proj][cell], MPI.SUM, comm)
    end
end
