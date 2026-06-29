@timeit timer "Calc_Smk_Umnk" function Calc_Smk_Umnk(BANDNUM, Amnk, Ngsize, material::LCPAO_model, kpoints::KPoints)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    MPI_Nkpt = kpoints.MPI_Nkpt
    
    Amatrix = zeros(ComplexF64, BANDNUM, Ngsize)
    Umnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    Σmk = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        Umnk[spin] = Vector{Matrix{ComplexF64}}(undef, MPI_Nkpt)
        Σmk[spin] = Vector{Vector{Float64}}(undef, MPI_Nkpt)
        for ik = 1:MPI_Nkpt
            Umnk[spin][ik] = zeros(ComplexF64, BANDNUM, Ngsize)
            Σmk[spin][ik] = zeros(Float64, Ngsize)
        end
    end
    

    # Calculate Uk using Uk = W*V† where A = WΣV†
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        for μ = 1:BANDNUM, p = 1:Ngsize
            Amatrix[μ,p] = Amnk[spin][ik][μ][p]
        end

        W, Σmk[spin][ik], V = svd(Amatrix)
        Umnk[spin][ik] = W * V'
    end


    return Σmk, Umnk
end