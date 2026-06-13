function Get_cell_list(kmesh)

    for i = 1:3
        if iseven(kmesh[1])
            error("please check kmesh")
        end
    end

    cell_list = Vector{UnitRange{Int64}}(undef, 3)
    for i = 1:3
        MinCell = 0
        MaxCell = 0
        if kmesh[i] < 0
            println(kmesh)
            error("error : kmesh mush positive values")
        elseif iszero(kmesh[i])
            MinCell = 0
            MaxCell = 0
        elseif iseven(kmesh[i])
            MinCell = -div(kmesh[i]-1, 2)
            MaxCell = div(kmesh[i], 2)
        elseif isodd(kmesh[i])
            MinCell = -div(kmesh[i], 2)
            MaxCell = div(kmesh[i], 2)
        end
        cell_list[i] = MinCell:MaxCell
    end

    NCell = 0
    for i in cell_list[1], j in cell_list[2], k in cell_list[3]
        NCell += 1
    end

    cell = 0
    cell_list_ijk = Vector{Vector{Int32}}(undef, NCell)
    for i in cell_list[1], j in cell_list[2], k in cell_list[3]
        cell += 1
        cell_list_ijk[cell] = [i, j, k]
    end
    
    return NCell, cell_list, cell_list_ijk
end


function CWF_weight(ϵnk, ChemP, Dis_Energy; δ=1e-12)

    ϵ0 = Dis_Energy[1]
    ϵ1 = Dis_Energy[2]
    kbT0 = Dis_Energy[3]
    kbT1 = Dis_Energy[4]

    b0 = 1/kbT0
    b1 = 1/kbT1
    e0 = ϵ0 + ChemP
    e1 = ϵ1 + ChemP

    x0 = b0*(e0 - ϵnk)
    x1 = b1*(ϵnk - e1)

    weight = (1/(exp(x0)+1))*(1/(exp(x1)+1)) + δ

    return weight
end


function CWF_weight2(ϵnk, ChemP, Dis_Energy)

    out0 = Dis_Energy[1]
    in0 = Dis_Energy[2]
    in1 = Dis_Energy[3]
    out1 = Dis_Energy[4]
    enk = ϵnk - ChemP

    c0 = 1.0
    c1 = 0.0
    c2 = -6.0
    c3 = 8.0
    c4 = -3.0

    if enk <= out0
        return 0.0
    elseif out0 < enk <= in0
        y = (in0-enk)/(in0-out0)
        return c0 + c1*y + c2*y*y + c3*y*y*y + c4*y*y*y*y
    elseif in0 < enk <= in1
        return 1.0
    elseif in1 < enk <= out1
        y = (enk-in1)/(out1-in1)
        return c0 + c1*y + c2*y*y + c3*y*y*y + c4*y*y*y*y
    elseif out1 < enk
        return 0.0
    end
end


function Find_MinN_MaxN(material, Dis_Energy, Enk)

    out0 = Dis_Energy[1]
    out1 = Dis_Energy[4]
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)

    MinN = 0
    MaxN = 0
    for μ = Nfsize:-1:1
        ene = Enk[1,μ,1] - ChemP
        if ene < out1
            MaxN = μ
            break
        end
    end

    for μ = 1:Nfsize
        ene = Enk[1,μ,1] - ChemP
        if out0 < ene
            MinN = μ
            break
        end
    end

    
    return MinN, MaxN
end

