include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_FT_ProductPAO_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    
    MaxL = pao.Spe_MaxL_Basis
    Mu = pao.Spe_Num_Basis

    Spe_ProductRF_Bessel = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, MaxL+1)
    for GL1 = 0:MaxL
        Spe_ProductRF_Bessel[GL1+1] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Mu[GL1+1])
        for Mul1 = 1:Mu[GL1+1]
            Spe_ProductRF_Bessel[GL1+1][Mul1] = Vector{Vector{Vector{Vector{Float64}}}}(undef, MaxL+1)
            for GL2 = 0:MaxL
                Spe_ProductRF_Bessel[GL1+1][Mul1][GL2+1] = Vector{Vector{Vector{Float64}}}(undef, Mu[GL2+1])
                    
                Lmax = ifelse(GL1 <= GL2, 2*GL2, 1)
                num = ifelse(GL1 <= GL2, GL_Mesh, 1)
                    
                for Mul2 = 1:Mu[GL2+1]
                    Spe_ProductRF_Bessel[GL1+1][Mul1][GL2+1][Mul2] = Vector{Vector{Float64}}(undef, Lmax+1)
                    for l = 0:Lmax
                        Spe_ProductRF_Bessel[GL1+1][Mul1][GL2+1][Mul2][l+1] = zeros(Float64, num)
                    end
                end
            end
        end
    end

    @code_warntype FT_ProductPAO!(pao, pspot.Spe_VPS_RV[begin], Spe_ProductRF_Bessel)
end

check_FT_ProductPAO_code_warntype()