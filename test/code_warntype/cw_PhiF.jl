include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_PhiF_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    pspot = Read_VPS("Si", "", "LDA", false; verbosity)
    maxL = 8
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    
    VNA_proj_ene = Vector{Vector{Float64}}(undef, maxL+1)
    for L = 0:maxL
        VNA_proj_ene[L+1] = zeros(Float64, maxM)
    end
        

    Projector_VNA = Vector{Vector{Vector{Float64}}}(undef, maxL+1)
    for L = 0:maxL
        Projector_VNA[L+1] = Vector{Vector{Float64}}(undef, maxM)
        for m = 1:maxM
            Projector_VNA[L+1][m] = zeros(Float64, Spe_Num_Mesh_VPS)
        end
    end
    Calc_ProExpn_VNA!(pao, pspot, maxL, Projector_VNA, VNA_proj_ene)


    @show typeof(Spe_Num_Mesh_VPS)
    @code_warntype PhiF(Spe_Num_Mesh_VPS, 1.0, Spe_VPS_RV, Projector_VNA[1][1])
end

check_PhiF_code_warntype()