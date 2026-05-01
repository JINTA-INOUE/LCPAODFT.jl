include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_FT_ProExpn_VNA_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    
    maxL = 2 + BufferL_ProVNA
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Projector_VNA = Vector{Vector{Vector{Float64}}}(undef, maxL+1)
    for L = 0:maxL
        Projector_VNA[L+1] = Vector{Vector{Float64}}(undef, maxM)
        for m = 1:maxM
            Projector_VNA[L+1][m] = zeros(Float64, Spe_Num_Mesh_VPS)
        end
    end

    VNA_proj_ene = Vector{Vector{Float64}}(undef, maxL+1)
    Spe_VNA_Bessel = Vector{Vector{Vector{Float64}}}(undef, maxL+1)
    for L = 0:maxL
        VNA_proj_ene[L+1] = zeros(Float64, maxM)
        Spe_VNA_Bessel[L+1] = Vector{Vector{Float64}}(undef, maxM)
        for m = 1:maxM
            Spe_VNA_Bessel[L+1][m] = zeros(Float64, GL_Mesh)
        end
    end
    @code_warntype FT_ProExpn_VNA!(pao, pspot, maxL, Projector_VNA, Spe_VNA_Bessel)
end

check_FT_ProExpn_VNA_code_warntype()