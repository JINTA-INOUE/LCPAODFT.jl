include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_Calc_Bessel_Pro00_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    
    Spe_MaxL_Basis = pao.Spe_MaxL_Basis
    Spe_Num_Basis = pao.Spe_Num_Basis
    
    Bessel_Pro00 = Vector{Vector{Vector{Float64}}}(undef, Spe_MaxL_Basis+1)
    for l = 0:Spe_MaxL_Basis
        Bessel_Pro00[l+1] = Vector{Vector{Float64}}(undef, Spe_Num_Basis[l+1])
        for p = 1:Spe_Num_Basis[l+1]
            Bessel_Pro00[l+1][p] = zeros(Float64, GL_Mesh)
        end
    end

    @code_warntype Calc_Bessel_Pro00!(pao, Bessel_Pro00)
end

check_Calc_Bessel_Pro00_code_warntype()