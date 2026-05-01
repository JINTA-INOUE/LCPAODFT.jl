include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "LCPAODFT.jl" begin
    include("test_Gauss_Legendre.jl")
    include("test_Get_Atoms_data.jl")
    include("test_Set_Comp2Real.jl")
    include("test_Ylm.jl")
    include("test_DFT_Setup.jl")
    include("test_KPoints.jl")
    include("test_Read_PAO.jl")
    include("test_Read_VPS.jl")
    include("test_UCell.jl")
    include("test_Set_AdenPCC_Grid.jl")
    include("test_Set_Orbs_Grid.jl")
    include("test_Set_XC_Grid.jl")
end