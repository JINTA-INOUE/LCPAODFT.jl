function Calc_psi_HOMO(material::LCPAO_model, kpts::Vector{Vector{Float64}})

    Natom = material.Natom
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    E_Temp = material.E_Temp
    ChemP = material.ChemP
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks

    fsize = sum(Total_NumOrbs)
    Nkpt = length(kpts)

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

    if SpinPol ∈ ("on", "nc")
        Spindeg = 2
    elseif SpinPol == "off"
        Spindeg = 1
    end



    Enk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Enk[spin] = Vector{Vector{Float64}}(undef, Nkpt)
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, Nkpt)
        for k = 1:Nkpt
            Enk[spin][k] = zeros(Float64, Nfsize)
            Cnk[spin][k] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end


    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for spin = 1:spinsize, ik = 1:Nkpt
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])

            Enk[spin][ik], Cnk[spin][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    else
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:Nkpt
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[ik])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk[1][ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
        end
    end



    Beta = 1/E_Temp/kb*eV2Hartree
    max_x = 30.0
    Bulk_HOMO = Vector{Vector{Int32}}(undef, Nkpt)
    for ik = 1:Nkpt
        Bulk_HOMO[ik] = zeros(Int32, 2)
    end
    for spin = 1:spinsize, ik = 1:Nkpt
        for μ = 1:Nfsize
            x = (Enk[spin][ik][μ] - ChemP)*Beta
            x = ifelse(x <= -max_x, -max_x, x)
            x = ifelse(x >=  max_x,  max_x, x)

            FermiF = 1/(1 + exp(x))
            if FermiF > 0.5
                Bulk_HOMO[ik][spin] = μ
            end
        end
    end

    if SpinPol == "nc"
        for ik = 1:Nkpt
            Bulk_HOMO[ik][2] = deepcopy(Bulk_HOMO[ik][1])
        end
    end




    NHOMO = 1    
    HOMOs_Coef = Vector{Vector{Vector{Vector{Vector{ComplexF64}}}}}(undef, Nkpt)
    for ik = 1:Nkpt
        HOMOs_Coef[ik] = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Spindeg)
        for spin = 1:Spindeg
            HOMOs_Coef[ik][spin] = Vector{Vector{Vector{ComplexF64}}}(undef, NHOMO)
            for μ = 1:NHOMO
                HOMOs_Coef[ik][spin][μ] = Vector{Vector{ComplexF64}}(undef, Natom)
                for atom = 1:Natom
                    HOMOs_Coef[ik][spin][μ][atom] = zeros(ComplexF64, Total_NumOrbs[atom])
                end
            end
        end
    end


    if SpinPol ∈ ("off", "on")
        for ik = 1:Nkpt, spin = 1:Spindeg
            μ1 = Bulk_HOMO[ik][spin]
            for atom = 1:Natom
                Anum = MP[atom]
                for ist = 1:Total_NumOrbs[atom]
                    HOMOs_Coef[ik][spin][1][atom][ist] = Cnk[spin][ik][Anum+ist,μ1]
                end
            end
        end
    else
        for ik = 1:Nkpt
            μ1 = Bulk_HOMO[ik][1]
            for atom = 1:Natom
                Anum = MP[atom]
                for ist = 1:Total_NumOrbs[atom]
                    HOMOs_Coef[ik][1][1][atom][ist] = Cnk[1][ik][Anum+ist,μ1]
                    HOMOs_Coef[ik][2][1][atom][ist] = Cnk[1][ik][Anum+ist+fsize,μ1]
                end
            end
        end
    end


    return Enk, Bulk_HOMO, HOMOs_Coef
end


function Calc_hwfs_HOMO(material::LCPAO_model, Umjk)

    Natom = material.Natom
    SpinPol = material.SpinPol
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    Atoms_Core_Charge = material.Atoms_Core_Charge
    Valence_Electrons = sum(Atoms_Core_Charge)
    Total_SpinS = 0.0
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks
    fsize = sum(Total_NumOrbs)

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

    if SpinPol ∈ ("on", "nc")
        Spindeg = 2
    elseif SpinPol == "off"
        Spindeg = 1
    end

    if SpinPol == "off"
        Nocc = Int64(div(Valence_Electrons,2))
    elseif SpinPol == "on"
        Nocc = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
    elseif SpinPol == "nc"
        Nocc = Int64(Valence_Electrons)
    end



    Enk = zeros(Float64, spinsize, Nfsize)
    Cnk = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = zeros(ComplexF64, Nfsize, Nfsize)
    end


    kpts = zeros(Float64, 3)
    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for spin = 1:spinsize
            HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
            Enk[spin,:], Cnk[spin] = eigen(Hermitian(H), Hermitian(S))
        end
    else
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
        @. S[1:fsize, 1:fsize] = tmpH
        @. S[fsize+1:end, fsize+1:end] = tmpH
        Enk[1,:], Cnk[1] = eigen(Hermitian(H), Hermitian(S))
    end


    # Calc C̃j,iα(k) =  ∑_{μ} Ujμ(k)Cμ,iα(k) (|hjk> = ∑_{μ} Ujμ(k)|ψμk>)
    Umjk2 = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Umjk2[spin] = zeros(ComplexF64, Nocc, Nocc)
    end

    for spin = 1:spinsize, j = 1:Nocc, μ = 1:Nocc
        Umjk2[spin][μ,j] = Umjk[spin][j,μ]
    end


    Cnk2 = Vector{Matrix{ComplexF64}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk2[spin] = zeros(ComplexF64, Nfsize, Nocc)
    end

    for spin = 1:spinsize, ist = 1:Nfsize, j = 1:Nocc
        Sum = 0.0 + im*0.0
        for μ = 1:Nocc
            Sum += Cnk[spin][ist,μ]*Umjk2[spin][μ,j]
        end
        Cnk2[spin][ist,j] = Sum
    end




    NHOMO = Nocc
    HOMOs_Coef = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Spindeg)
    for spin = 1:Spindeg
        HOMOs_Coef[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, NHOMO)
        for μ = 1:NHOMO
            HOMOs_Coef[spin][μ] = Vector{Vector{ComplexF64}}(undef, Natom)
            for atom = 1:Natom
                HOMOs_Coef[spin][μ][atom] = zeros(ComplexF64, Total_NumOrbs[atom])
            end
        end
    end


    if SpinPol ∈ ("off", "on")
        for spin = 1:Spindeg, μ = 1:NHOMO
            for atom = 1:Natom
                NO0 = Total_NumOrbs[atom]
                Anum = MP[atom]
                for ist = 1:NO0
                    HOMOs_Coef[spin][μ][atom][ist] = Cnk2[spin][Anum+ist,μ]
                end
            end
        end
    else
        for μ = 1:NHOMO
            for atom = 1:Natom
                NO0 = Total_NumOrbs[atom]
                Anum = MP[atom]
                for ist = 1:NO0
                    HOMOs_Coef[1][μ][atom][ist] = Cnk2[1][Anum+ist,μ]
                    HOMOs_Coef[2][μ][atom][ist] = Cnk2[1][Anum+ist+fsize,μ]
                end
            end
        end
    end


    return HOMOs_Coef
end