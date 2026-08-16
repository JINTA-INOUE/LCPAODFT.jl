@timeit timer "Mixing_H!"  function Mixing_H!(SCF_iter, MPI_Hks, Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    Nspin = dft_mixing.Nspin
    MPI_Hsize = dft_mixing.MPI_Hsize
    Start_Pulay_SCF = dft_options.Start_Pulay_SCF
    if SCF_iter <= Start_Pulay_SCF-1
        Simple_Mixing_H!(SCF_iter, MPI_Hks, dft_options, dft_mixing)
    else
        Pulay_Mixing_H!(SCF_iter, MPI_Hks, Hks, dft_options, dft_mixing)
    end

    for spin = 1:Nspin
        MPI.Allgatherv!(MPI_Hks[spin], VBuffer(Hks[spin], MPI_Hsize), comm)
    end

    if SCF_iter == 1
        dft_mixing.NormRD[1] = 1.0
    end
end


function Simple_Mixing_H!(SCF_iter, MPI_Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = dft_mixing.Natom
    Nspin = dft_mixing.Nspin
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay
    dim = ifelse(SCF_iter<Num_Mixing_Pulay, SCF_iter, Num_Mixing_Pulay)
    ResH = dft_mixing.ResH
    HisH = dft_mixing.HisH
    MPI_Hsize = dft_mixing.MPI_Hsize
    myHsize = MPI_Hsize[myrank+1]


    # shift Residual H
    for m = dim:-1:2, spin = 1:Nspin, hst = 1:myHsize
        ResH[m][hst,spin] = ResH[m-1][hst,spin]
    end



    # calculate the current Residual H
    for spin = 1:Nspin, hst = 1:myHsize
        ResH[1][hst,spin] = MPI_Hks[spin][hst] - HisH[1][hst,spin]
    end


    Norm = H_norm(Nspin, myHsize, MPI_Hks, HisH[1])
    Norm = MPI.Allreduce(Norm, MPI.SUM, comm)
    Norm = Norm/Natom
    update_NormRD!(Norm, dft_mixing)


    # shift Historical H
    for m = dim+1:-1:2, spin = 1:Nspin, hst = 1:myHsize
        HisH[m][hst,spin] = HisH[m-1][hst,spin]
    end


    if SCF_iter ≠ 1
        get_Mixing_weight!(dft_options, dft_mixing)
        Mixing_weight = dft_options.Mixing_weight
        Hmix!(Nspin, myHsize, Mixing_weight, MPI_Hks, HisH[2])
    end


    for spin = 1:Nspin, hst = 1:myHsize
        HisH[1][hst,spin] = MPI_Hks[spin][hst]
    end

    dft_mixing.ResH = ResH
    dft_mixing.HisH = HisH
end


function Pulay_Mixing_H!(SCF_iter, MPI_Hks, Hks, dft_options::DFT_Options, dft_mixing::Ham_Mixing)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = dft_mixing.Natom
    Nspin = dft_mixing.Nspin
    FNAN = dft_mixing.FNAN
    natn = dft_mixing.natn
    Total_NumOrbs = dft_mixing.Total_NumOrbs
    MPI_size = dft_mixing.MPI_size
    MPI_atom = dft_mixing.MPI_atom
    MPI_natn = dft_mixing.MPI_natn
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay
    dim = ifelse(SCF_iter<=Num_Mixing_Pulay, SCF_iter-1, Num_Mixing_Pulay)
    ChemP = dft_mixing.ChemP
    ResH = dft_mixing.ResH
    HisH = dft_mixing.HisH
    MPI_Hsize = dft_mixing.MPI_Hsize
    myHsize = MPI_Hsize[myrank+1]


    # shift Residual H
    for m = dim+1:-1:2, spin = 1:Nspin, hst = 1:myHsize
        ResH[m][hst,spin] = ResH[m-1][hst,spin]
    end

    # calculate the current Residual H
    for spin = 1:Nspin, hst = 1:myHsize
        ResH[1][hst,spin] = MPI_Hks[spin][hst] - HisH[1][hst,spin]
    end

    metric = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        metric[atom] = zeros(Float64, Total_NumOrbs[atom])
    end
    MPI.Allgatherv!(HisH[1], VBuffer(Hks[1], MPI_Hsize), comm)
    get_metric!(Natom, FNAN, natn, Total_NumOrbs, metric, Hks[1], ChemP)



    A = zeros(Float64, Num_Mixing_Pulay+3, Num_Mixing_Pulay+3)

    for m = 1:dim, n = m:dim
        MPI_Sum = 0.0
        ResHm = ResH[m]
        ResHn = ResH[n]
        for spin = 1:Nspin
            hst = 0
            for loop = 1:MPI_size
                atom = MPI_atom[loop]
                jatom = MPI_natn[loop]
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                for ist = 1:NO0, jst = 1:NO1
                    hst += 1
                    MPI_Sum += metric[atom][ist]*ResHm[hst,spin]*ResHn[hst,spin]
                end
            end
        end
        Sum = MPI.Allreduce(MPI_Sum, MPI.SUM, comm)
        A[m,n] = Sum
        A[n,m] = A[m,n]
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
        fill!(coes, 0.0)
        coes[1] = 0.05
        coes[2] = 0.95
    end

    

    # mixing Hamiltonian
    if Norm >= 1e-1
        alpha = 0.5
    elseif 1e-2 <= Norm < 1e-1
        alpha = 0.6
    elseif 1e-3 <= Norm < 1e-2
        alpha = 0.7
    elseif 1e-4 <= Norm < 1e-3
        alpha = 0.8
    else
        alpha = 1.0
    end


    for spin = 1:Nspin, hst = 1:myHsize
        r = 0.0
        h = 0.0
        for m = 1:dim
            r += ResH[m][hst,spin]*coes[m]
            h += HisH[m][hst,spin]*coes[m]
        end
        ResH[dim+1][hst,spin] = r
        MPI_Hks[spin][hst] = h + alpha*r
    end



    # shift Historical H
    for m = dim:-1:2, spin = 1:Nspin, hst = 1:myHsize
        HisH[m][hst,spin] = HisH[m-1][hst,spin]
    end

    for spin = 1:Nspin, hst = 1:myHsize
        HisH[1][hst,spin] = MPI_Hks[spin][hst]
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


function H_norm(Nspin, myHsize, MPI_Hks, HisH)
    Norm = 0.0
    for spin = 1:Nspin, hst = 1:myHsize
        tmp = MPI_Hks[spin][hst] - HisH[hst,spin]
        Norm += tmp^2
    end

    return Norm
end


function Hmix!(Nspin, myHsize, weight, MPI_Hks, HisH)
    weight2 = 1.0 - weight
    for spin = 1:Nspin, hst = 1:myHsize
        MPI_Hks[spin][hst] = weight2*HisH[hst,spin] + weight*MPI_Hks[spin][hst]
    end
end

#=
function Pulay_H_inv!(dim, IA)
    val, vec = eigen(Symmetric(IA))
    for i = 1:dim
        val[i] = 1/(val[i] + 1.0e-13)
    end

    for i = 1:dim, j = 1:dim
        Sum = 0.0
        for k = 1:dim
            Sum += vec[i,k]*val[k]*vec[j,k]
        end
        IA[i,j] = Sum
    end
end
=#