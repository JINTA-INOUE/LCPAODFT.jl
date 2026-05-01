include("../../src/LCPAODFT.jl")
using .LCPAODFT


function check_RF_BesselF_code_warntype()
    verbosity = 0
    pao = Read_PAO(4.0, "Si", 7.0, "s2p2d1", ""; verbosity)
    Spe_RF_Bessel = pao.Spe_RF_Bessel

    kmax = Nkmax
    dk = kmax/NkGrid
    NormK = zeros(Float64, NkGrid)
    for ik = 1:NkGrid
        NormK[ik] = (ik-1)*dk
    end

    @code_warntype RF_BesselF(1.0, NormK, Spe_RF_Bessel[1][1])
end

check_RF_BesselF_code_warntype()