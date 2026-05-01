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


function weight_CWF(ϵnk, kbT0, kbT1, ϵ0, ϵ1; δ=1.0e-10, ecut=10.0)
    x0 = (ϵ0 - ϵnk)/kbT0
    x1 = (ϵnk - ϵ1)/kbT1

    tmp0 = ifelse(x0 > ecut, exp(-x0), 1/(1+exp(x0)))
    tmp1 = ifelse(x1 > ecut, exp(-x1), 1/(1+exp(x1)))
    tmp2 = ifelse(x0+x1 > 2*ecut, -exp(2*ecut), 1-exp(x0 + x1))

    return tmp2*tmp0*tmp1 + δ
end


function Calc_Enk_Cnk(material::LCPAO_model, kpoints::KPoints)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)

    Hks = material.Hks
    iHks = material.iHks
    OLP = material.OLP


    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    else SpinPol == "on"
        spinsize = 2
    end



    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    if SpinPol ∈ ("off", "on")
        Enk = zeros(Float64, spinsize, fsize, MPI_Nkpt)
        for spin = 1:spinsize
            Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
            for ik = 1:MPI_Nkpt
                Cnk[spin][ik] = zeros(ComplexF64, fsize, fsize)
            end
        end
    elseif SpinPol == "nc"
        Enk = zeros(Float64, 1, 2*fsize, MPI_Nkpt)
        Cnk[1] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Cnk[1][ik] = zeros(ComplexF64, 2*fsize, 2*fsize)
        end
    end



    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        @inbounds for spin = 1:spinsize, ik = 1:MPI_Nkpt
            HS_matrix!(S, material, OLP, MPI_kpts[ik])
            HS_matrix!(H, material, Hks[spin], MPI_kpts[ik])
            Enk[spin,:,ik], Cnk[spin][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:MPI_Nkpt
            HS_matrix_NC!(tmpH, H, Hks, iHks, material, MPI_kpts[ik])
            HS_matrix!(tmpH, material, OLP, MPI_kpts[ik])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk[1,:,ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    end

    return Enk, Cnk
end


function Calc_Amnk(kBT, ε, Guide_Orbs, Ngsize, Enk, Cnk, material, kpoints)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    MP = material.MP
    Total_NumOrbs = material.Total_NumOrbs
    SpinPol = material.SpinPol
    fsize = sum(Total_NumOrbs)
    ChemP = material.ChemP
    OLP = material.OLP


    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
    end

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end

    Guide_Total_NumOrbs = zeros(Int32, Natom)
    for atom = 1:Natom
        Guide_Total_NumOrbs[atom] = length(Guide_Orbs[atom])
    end
    gsize = sum(Guide_Total_NumOrbs)



    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, Nfsize)
            for μ = 1:Nfsize
                Amnk[spin][ik][μ] = zeros(ComplexF64, Ngsize)
            end
        end
    end


    # Calculation Akμ,p
    if SpinPol == "off"
        for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            pst = 0
            for patom = 1:Natom, proj in Guide_Orbs[patom]
                
                pst += 1
                Sum = ComplexF64(0.0, 0.0)
    
                for Rn = 1:FNAN[patom]+1
    
                    jatom = natn[patom][Rn]
                    cell = ncn[patom][Rn]+1
                    kRn = dot(MPI_kpts[ik], atv_ijk[cell])
                    ex = cispi(-2*kRn)
                    Bnum = MP[jatom]
                    tmp1 = ComplexF64(0.0, 0.0)
                    for jβ = 1:Total_NumOrbs[jatom]
                        tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ])*OLP[patom][Rn][proj][jβ]
                    end
                    Sum += tmp1*ex
                end
    
                weight = weight_CWF((Enk[1,μ,ik]-ChemP)*Hartree2eV,kBT[1],kBT[2], ε[1], ε[2])
                Amnk[1][ik][μ][pst] = weight * Sum
            end
        end
    elseif SpinPol == "on"
        for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            pst = 0
            for patom = 1:Natom, proj in Guide_Orbs[patom]
                
                pst += 1
                Sum1 = ComplexF64(0.0, 0.0)
                Sum2 = ComplexF64(0.0, 0.0)
    
                for Rn = 1:FNAN[patom]+1
    
                    jatom = natn[patom][Rn]
                    cell = ncn[patom][Rn]+1
                    kRn = dot(MPI_kpts[ik], atv_ijk[cell])
                    ex = cispi(-2*kRn)
                    Bnum = MP[jatom]
                    tmp1 = ComplexF64(0.0, 0.0)
                    tmp2 = ComplexF64(0.0, 0.0)
                    for jβ = 1:Total_NumOrbs[jatom]
                        tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ])*OLP[patom][Rn][proj][jβ]
                        tmp2 += conj(Cnk[2][ik][Bnum+jβ,μ])*OLP[patom][Rn][proj][jβ]
                    end
                    Sum1 += tmp1*ex
                    Sum2 += tmp2*ex
                end
    
                weight1 = weight_CWF((Enk[1,μ,ik]-ChemP)*Hartree2eV,kBT[1], kBT[2], ε[1], ε[2])
                weight2 = weight_CWF((Enk[2,μ,ik]-ChemP)*Hartree2eV,kBT[1], kBT[2], ε[1], ε[2])
                Amnk[1][ik][μ][pst] = weight1 * Sum1
                Amnk[2][ik][μ][pst] = weight2 * Sum2
            end
        end
    elseif SpinPol == "nc"
        for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            pst = 0
            for patom = 1:Natom, proj in Guide_Orbs[patom]
                
                pst += 1
                Sum1 = ComplexF64(0.0, 0.0)
                Sum2 = ComplexF64(0.0, 0.0)
    
                for Rn = 1:FNAN[patom]+1
    
                    jatom = natn[patom][Rn]
                    cell = ncn[patom][Rn]+1
                    kRn = dot(MPI_kpts[ik], atv_ijk[cell])
                    ex = cispi(-2*kRn)
                    Bnum = MP[jatom]
                    tmp1 = ComplexF64(0.0, 0.0)
                    tmp2 = ComplexF64(0.0, 0.0)
                    for jβ = 1:Total_NumOrbs[jatom]
                        tmp1 += conj(Cnk[1][ik][Bnum+jβ,μ])*OLP[patom][Rn][proj][jβ]
                        tmp2 += conj(Cnk[1][ik][fsize+Bnum+jβ,μ])*OLP[patom][Rn][proj][jβ]
                    end
                    Sum1 += tmp1*ex
                    Sum2 += tmp2*ex
                end
    
                weight = weight_CWF((Enk[1,μ,ik]-ChemP)*Hartree2eV,kBT[1],kBT[2], ε[1], ε[2])
                Amnk[1][ik][μ][pst] = weight * Sum1
                Amnk[1][ik][μ][gsize+pst] = weight * Sum2
            end
        end
    end


    return Amnk
end



function Calc_Σmk_Umnk(Amnk, Ngsize, material, kpoints)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MPI_Nkpt = kpoints.MPI_Nkpt

    if SpinPol ∈ ("off", "on")
        Nfsize = fsize
    elseif SpinPol == "nc"
        Nfsize = 2*fsize
    end

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end
    
    
    Amatrix = zeros(ComplexF64, Nfsize, Ngsize)
    Umnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    Σmk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        Σmk[spin] = Vector{Vector{Float64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Umnk[spin][ik] = zeros(ComplexF64, Nfsize, Ngsize)
            Σmk[spin][ik] = zeros(Float64, Ngsize)
        end
    end
    

    # Calculate Uk using Uk = W*V† where A = WΣV†
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        for μ = 1:Nfsize, p = 1:Ngsize
            Amatrix[μ,p] = Amnk[spin][ik][μ][p]
        end

        W, Σmk[spin][ik], V = svd(Amatrix)
        Umnk[spin][ik] = W * V'
    end


    return Σmk, Umnk
end


function Calc_DM(spinsize, Ngsize, Σmk, kpoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Nkpt = kpoints.Nkpt
    MPI_Nkpt = kpoints.MPI_Nkpt


    DM = 0.0
    for spin = 1:spinsize, ik = 1:MPI_Nkpt, p = 1:Ngsize
        if Σmk[spin][ik][p] < 0.0
            error("singler values is negative")
        end
        DM += (Σmk[spin][ik][p]-1)^2
    end
    DM = MPI.Allreduce(DM, MPI.SUM, comm)
    DM = DM/Nkpt/Ngsize/spinsize

    myrank == 0 && println("DM per CWF is $(DM)")


    return DM
end


function Calc_HmnR!(HmnR, spinsize, Nfsize, NCell, Ngsize, Nkpt, MPI_Nkpt, MPI_kpts, cell_list_ijk, Umnk, Enk)

    for spin = 1:spinsize, cell = 1:NCell, pst = 1:Ngsize, qst = 1:Ngsize
                
        Sum = ComplexF64(0.0, 0.0)
        for ik = 1:MPI_Nkpt
            tmp = ComplexF64(0.0, 0.0)
            for μ = 1:Nfsize
                tmp += Enk[spin,μ,ik] * conj(Umnk[spin][ik][μ,pst]) * Umnk[spin][ik][μ,qst]
            end

            kRn = dot(MPI_kpts[ik], cell_list_ijk[cell])
            ex = cispi(-2*kRn)
            Sum += tmp*ex
        end

        HmnR[spin,cell,pst,qst] = Sum/Nkpt
    end

    HmnR = MPI.Allreduce(HmnR, MPI.SUM, comm)
end