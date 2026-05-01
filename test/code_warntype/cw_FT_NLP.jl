include("../../src/LCPAODFT.jl")
using .LCPAODFT
using Printf

function check_FT_NLP_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pao = Read_PAO(1.0, "H", 5.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    pspot = Read_VPS("H", "", "LDA", false; verbosity)
    Spe_Num_RVPS = pspot.Spe_Num_RVPS
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    VPS_j_dependency = pspot.VPS_j_dependency
    NLRF_Bessel = Vector{Vector{Vector{Float64}}}(undef, VPS_j_dependency+1)
    for so = 1:VPS_j_dependency+1
        NLRF_Bessel[so] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
        for l = 1:Spe_Num_RVPS
            NLRF_Bessel[so][l] = zeros(Float64, NkGrid+1)
        end
    end
    FT_NLP!(Spe_Atom_Cut1, NLRF_Bessel, pspot)
end

check_FT_NLP_code_warntype()