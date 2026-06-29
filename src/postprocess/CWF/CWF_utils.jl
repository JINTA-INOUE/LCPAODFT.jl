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

    ϵ0, ϵ1, kbT0, kbT1 = Dis_Energy
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

    out0, in0, in1, out1 = Dis_Energy
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


function Calc_MLWF_Enk(MLWF_kpts::Vector{Float64}, material::LCPAO_model)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    Hks = material.Hks
    iHks = material.iHks
    OLP = material.OLP

    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MLWF_kpts)
        MLWF_Enk = eigvals(Hermitian(H), Hermitian(S))
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        HS_matrix_NC!(H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MLWF_kpts)
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MLWF_kpts)
        @. S[1:fsize, 1:fsize] = tmpH
        @. S[fsize+1:end, fsize+1:end] = tmpH
        MLWF_Enk = eigvals(Hermitian(H), Hermitian(S))
    end


    return MLWF_Enk
end


function Find_MinN_MaxN(Dis_Energy, Nfsize, MLWF_Enk, ChemP)

    out0 = Dis_Energy[1]
    out1 = Dis_Energy[4]

    MinN = 0
    MaxN = 0
    for μ = Nfsize:-1:1
        ene = MLWF_Enk[μ] - ChemP
        if ene < out1
            MaxN = μ
            break
        end
    end

    for μ = 1:Nfsize
        ene = MLWF_Enk[μ] - ChemP
        if out0 < ene
            MinN = μ
            break
        end
    end

    
    return MinN, MaxN
end


function Calc_BANDNUM_KS_state(weight_type::String, MLWF_kpts, Dis_Energy, material::LCPAO_model)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    ChemP = material.ChemP

    if weight_type == "poly"
        MLWF_Enk = Calc_MLWF_Enk(MLWF_kpts, material)
        MinN, MaxN = Find_MinN_MaxN(Dis_Energy, Nfsize, MLWF_Enk, ChemP)
        BANDNUM = MaxN - MinN + 1
    else
        MinN = 1
        MaxN = Nfsize
        BANDNUM = Nfsize
    end


    return MinN, MaxN, BANDNUM
end


function Check_CWF_Band(weight_type::String, MinN, MaxN, BANDNUM, Ngsize)

    if weight_type == "poly"

        println("\tMinN = $MinN  MaxN = $MaxN  BANDNUM = $BANDNUM  Num_CWFs = $Ngsize")

        if MaxN < MinN
            println("Could not find any state to be included.")
            error("Parameters for CWF disentangling must be improper.")
        end

        if BANDNUM < Ngsize
            println("(MaxN-MinN+1) should be larger than TNum_CWFs.")
            error("Parameters for CWF disentangling must be improper.")
        end
    end
end