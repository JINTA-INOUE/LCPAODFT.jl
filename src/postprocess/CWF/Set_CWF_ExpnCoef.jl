function Set_CWF_ExpnCoef(GNatom, Guide_Total_NumOrbs, CWF_Plot_SuperCells, Cnk, Umnk, material, kmesh)

    println("Set_CWF_ExpnCoef CWF_SOC on site")

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)

    kpoints = Set_KPoints(kmesh, false)

    Nkpt = kpoints.Nkpt
    kpts = kpoints.kpts


    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    elseif SpinPol == "nc"
        error("please check")
    end


    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end

    gsize = sum(Guide_Total_NumOrbs)

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



    if SpinPol ∈ ("off", "on")
        CWF_ExpnCoef = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, spinsize)
        for spin = 1:spinsize
            CWF_ExpnCoef[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, GNatom)
            for atom = 1:GNatom
                CWF_ExpnCoef[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, Guide_Total_NumOrbs[atom])
                for pst = 1:Guide_Total_NumOrbs[atom]
                    CWF_ExpnCoef[spin][atom][pst] = Vector{Vector{Float64}}(undef, Plot_NCell)
                    for cell = 1:Plot_NCell
                        CWF_ExpnCoef[spin][atom][pst][cell] = zeros(Float64, Nfsize)
                    end
                end
            end
        end
    elseif SpinPol == "nc"
        CWF_ExpnCoef = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, GNatom)
        for atom = 1:GNatom
            CWF_ExpnCoef[atom] = Vector{Vector{Vector{ComplexF64}}}(undef, 2*Guide_Total_NumOrbs[atom])
            for pst = 1:2*Guide_Total_NumOrbs[atom]
                CWF_ExpnCoef[atom][pst] = Vector{Vector{ComplexF64}}(undef, Plot_NCell)
                for cell = 1:Plot_NCell
                    CWF_ExpnCoef[atom][pst][cell] = zeros(ComplexF64, Nfsize)
                end
            end
        end
    end



    # CWF_ExpnCoef = ∑_{kμ}e^{ikR}u_{kμ,p}c_{kμ,iα}
    if SpinPol ∈ ("off", "on")
        for spin = 1:spinsize
            proj = 0
            for atom = 1:GNatom, pst = 1:Guide_Total_NumOrbs[atom]
                proj += 1
                for cell = 1:Plot_NCell
                    l, m, n = Plot_cell_ijk[cell]
                    for ist = 1:Nfsize
                        Sum = ComplexF64(0.0, 0.0)
                        for ik = 1:Nkpt
                            kRn = kpts[ik][1]*l + kpts[ik][2]*m + kpts[ik][3]*n
                            ex = cispi(2*kRn)/Nkpt
                            temp = ComplexF64(0.0, 0.0)
                            @inbounds for μ = 1:Nfsize
                                temp += Umnk[spin][ik][μ,proj]*Cnk[spin][ik][ist,μ]
                            end
                            Sum += temp * ex
                        end

                        CWF_ExpnCoef[spin][atom][pst][cell][ist] = real(Sum)
                    end
                end
            end
        end
    else SpinPol == "nc"
        error("please check")
        #=
        for spin = 1:2
            proj = 0
            spin_site = ifelse(spin==1, 0, gsize)
            for atom = 1:GNatom, pst = 1:Guide_Total_NumOrbs[atom]
                proj += 1
                Anum = ifelse(spin==1, 0, Guide_Total_NumOrbs[atom])
                for cell = 1:Plot_NCell
                    l, m, n = Plot_cell_ijk[cell]
                    for ist = 1:Nfsize
                        Sum = ComplexF64(0.0, 0.0)
                        for ik = 1:Nkpt
                            kRn = kpts[ik][1]*l + kpts[ik][2]*m + kpts[ik][3]*n
                            ex = cispi(2*kRn)/Nkpt
                            temp = ComplexF64(0.0, 0.0)
                            @inbounds for μ = 1:Nfsize
                                temp += Umnk[1][ik][μ,spin_site+proj]*Cnk[1][ik][ist,μ]
                            end
                            Sum += temp * ex
                        end

                        CWF_ExpnCoef[atom][Anum+pst][cell][ist] = Sum
                    end
                end
            end
        end
        =#
    end


    @show [0,0,0]
    for atom = 1:GNatom, pst = 1:Guide_Total_NumOrbs[atom], cell = 1:Plot_NCell
        if Plot_cell_ijk[cell] == [0,0,0]
            @show atom, pst
            for ist = 1:Nfsize
                @show ist, CWF_ExpnCoef[1][atom][pst][cell][ist]
            end
        end
    end
    #=
    @show [1,0,0]
    for cell = 1:Plot_NCell
        if Plot_cell_ijk[cell] == [1,0,0]
            @show CWF_ExpnCoef[1][1][1][cell]
        end
    end
    @show [0,0,-1]
    for cell = 1:Plot_NCell
        if Plot_cell_ijk[cell] == [0,0,-1]
            @show CWF_ExpnCoef[1][1][1][cell]
        end
    end
    @show [0,1,0]
    for cell = 1:Plot_NCell
        if Plot_cell_ijk[cell] == [0,1,0]
            @show CWF_ExpnCoef[1][1][1][cell]
        end
    end
    =#

    return CWF_ExpnCoef
end