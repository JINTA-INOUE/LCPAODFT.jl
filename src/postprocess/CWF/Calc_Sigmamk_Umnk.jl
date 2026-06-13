function Calc_Smk_Umnk(Amnk, Ngsize, material::LCPAO_model, kpoints::KPoints)

    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MPI_Nkpt = kpoints.MPI_Nkpt
    Nks_size = length(Amnk[1][1])

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
    
    
    Amatrix = zeros(ComplexF64, Nks_size, Ngsize)
    Umnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    Σmk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        Σmk[spin] = Vector{Vector{Float64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Umnk[spin][ik] = zeros(ComplexF64, Nks_size, Ngsize)
            Σmk[spin][ik] = zeros(Float64, Ngsize)
        end
    end
    

    # Calculate Uk using Uk = W*V† where A = WΣV†
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        for μ = 1:Nks_size, p = 1:Ngsize
            Amatrix[μ,p] = Amnk[spin][ik][μ][p]
        end

        W, Σmk[spin][ik], V = svd(Amatrix)
        Umnk[spin][ik] = W * V'
    end


    return Σmk, Umnk
end