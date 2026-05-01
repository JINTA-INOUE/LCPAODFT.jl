function Mixing_H!(SCF_iter, Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    Start_Pulay_SCF = dft_options.Start_Pulay_SCF
    if SCF_iter <= Start_Pulay_SCF-1
        Simple_Mixing_H!(SCF_iter, Hks, dft_options, dft_mixing)
    else
        Pulay_Mixing_H!(SCF_iter, Hks, dft_options, dft_mixing)
    end

    if SCF_iter == 1
        dft_mixing.NormRD[1] = 1.0
    end
end


function Simple_Mixing_H!(SCF_iter, Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    Natom = dft_mixing.Natom
    Nspin = dft_mixing.Nspin
    Total_Hsize = dft_mixing.Total_Hsize
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay
    dim = ifelse(SCF_iter<Num_Mixing_Pulay, SCF_iter, Num_Mixing_Pulay)

    ResH = dft_mixing.ResH
    HisH = dft_mixing.HisH

    # shift Residual H
    for m = dim:-1:2, spin = 1:Nspin, hst = 1:Total_Hsize
        ResH[m][hst,spin] = ResH[m-1][hst,spin]
    end



    # calculate the current Residual H
    for spin = 1:Nspin, hst = 1:Total_Hsize
        ResH[1][hst,spin] = Hks[spin][hst] - HisH[1][hst,spin]
    end


    Norm = H_norm(Nspin, Total_Hsize, Hks, HisH[1])
    Norm = Norm/Natom
    update_NormRD!(Norm, dft_mixing)


    # shift Historical H
    for m = dim+1:-1:2, spin = 1:Nspin, hst = 1:Total_Hsize
        HisH[m][hst,spin] = HisH[m-1][hst,spin]
    end


    if SCF_iter ≠ 1
        get_Mixing_weight!(dft_options, dft_mixing)
        Mixing_weight = dft_options.Mixing_weight
        Hmix!(Nspin, Total_Hsize, Mixing_weight, Hks, HisH[2])
    end


    for spin = 1:Nspin, hst = 1:Total_Hsize
        HisH[1][hst,spin] = Hks[spin][hst]
    end

    dft_mixing.ResH = ResH
    dft_mixing.HisH = HisH
end


function Pulay_Mixing_H!(SCF_iter, Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    Natom = dft_mixing.Natom
    Nspin = dft_mixing.Nspin
    FNAN = dft_mixing.FNAN
    natn = dft_mixing.natn
    Total_NumOrbs = dft_mixing.Total_NumOrbs
    Total_Hsize = dft_mixing.Total_Hsize
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay
    dim = ifelse(SCF_iter<=Num_Mixing_Pulay, SCF_iter-1, Num_Mixing_Pulay)

    ChemP = dft_mixing.ChemP

    ResH = dft_mixing.ResH
    HisH = dft_mixing.HisH


    metric = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        metric[atom] = zeros(Float64, Total_NumOrbs[atom])
    end
    get_metric!(Natom, FNAN, natn, Total_NumOrbs, metric, HisH[1], ChemP)


    # shift Residual H
    for m = dim+1:-1:2, spin = 1:Nspin, hst = 1:Total_Hsize
        ResH[m][hst,spin] = ResH[m-1][hst,spin]
    end


    # calculate the current Residual H
    for spin = 1:Nspin, hst = 1:Total_Hsize
        ResH[1][hst,spin] = Hks[spin][hst] - HisH[1][hst,spin]
    end



    A = zeros(Float64, Num_Mixing_Pulay+3, Num_Mixing_Pulay+3)

    for m = 1:dim, n = m:dim
        Sum = 0.0
        for spin = 1:Nspin
            hst = 0
            for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
                hst += 1
                Sum += metric[atom][ist]*ResH[m][hst,spin]*ResH[n][hst,spin]
            end
            A[m,n] = Sum
            A[n,m] = A[m,n]
        end
    end

    Norm = A[1,1]/Natom
    update_NormRD!(Norm, dft_mixing)



    for m = 2:dim+1
        A[m-1,dim+1] = -1.0
        A[dim+1,m-1] = -1.0
    end
    A[dim+1,dim+1] = 0.0

    invA = inv(A[1:dim+1,1:dim+1])
    coes = zeros(Float64, dim)
    for m = 1:dim
        coes[m] = -invA[m,dim+1]
    end

    # NaN Inf check
    flag = false
    for m = 1:dim
        flag = flag || isnan(coes[m]) || isinf(coes[m])
    end

    if flag
        coes = zeros(Float64, dim)
        coes[1] = 0.05
        coes[2] = 0.95
    end
    
    

    # calculation of optimum Residual H
    for spin = 1:Nspin, hst = 1:Total_Hsize
        r = 0.0
        for m = 1:dim
            r += ResH[m][hst,spin]*coes[m]
        end
        ResH[dim+1][hst,spin] = r
    end

    

    # mixing Hamiltonian
    if Norm >= 1e-1
        alpha = 0.5
    elseif 1e-2 < Norm < 1e-1
        alpha = 0.6
    elseif 1e-3 < Norm < 1e-2
        alpha = 0.7
    elseif 1e-4 < Norm < 1e-3
        alpha = 0.8
    else
        alpha = 1.0
    end


    for spin = 1:Nspin, hst = 1:Total_Hsize
        r = 0.0
        h = 0.0
        for m = 1:dim
            r += ResH[m][hst,spin]*coes[m]
            h += HisH[m][hst,spin]*coes[m]
        end
        Hks[spin][hst] = h + alpha*r
    end



    # shift Historical H
    for m = dim:-1:2, spin = 1:Nspin, hst = 1:Total_Hsize
        HisH[m][hst,spin] = HisH[m-1][hst,spin]
    end

    for spin = 1:Nspin, hst = 1:Total_Hsize
        HisH[1][hst,spin] = Hks[spin][hst]
    end


    dft_mixing.ResH = ResH
    dft_mixing.HisH = HisH
end


function get_metric!(Natom, FNAN, natn, Total_NumOrbs, metric, HisH, ChemP)

    hks_sum = 0
    for atom = 1:Natom
        for ist = 1:Total_NumOrbs[atom]
            d = abs(HisH[hks_sum+ist,1] - ChemP)
            metric[atom][ist] = 5/(d^2 + 5)
        end

        for Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            hks_sum += 1
        end
    end
end


function H_norm(Nspin, Total_Hsize, Hks, HisH)
    Norm = 0.0
    for spin = 1:Nspin, hst = 1:Total_Hsize
        tmp = Hks[spin][hst] - HisH[hst,spin]
        Norm += tmp^2
    end

    return Norm
end


function Hmix!(Nspin, Total_Hsize, weight, Hks, HisH)
    weight2 = 1.0 - weight
    for spin = 1:Nspin, hst = 1:Total_Hsize
        Hks[spin][hst] = weight2*HisH[hst,spin] + weight*Hks[spin][hst]
    end
end
