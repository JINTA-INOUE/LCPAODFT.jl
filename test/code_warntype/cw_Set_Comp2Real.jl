include("../../src/utils/Set_Comp2Real.jl")

function check_Set_Comp2Real_code_warntype()
    
    C1 = zeros(ComplexF64, 3, 3)
    C2 = zeros(ComplexF64, 13, 13)
    C3 = zeros(ComplexF64, 18, 18)
    Spe_MaxL_Basis = 2
    Spe_Num_Basis = Int32.([2,2,1])
    Spe_VPS_List = Int32.([0, 0, 1, 1, 2, 2])

    # @code_warntype Set_Comp2Real(1)
    # @code_warntype Set_Comp2Real!(C1,1)
    # @code_warntype Set_Comp2Real(1)
    # @code_warntype Calc_Comp2Real(1)
    # @code_warntype Set_Comp2Real!(C2, Spe_MaxL_Basis, Spe_Num_Basis)
    # @code_warntype Set_NLComp2Real!(C3, Spe_VPS_List)
    @code_warntype Set_VNAComp2Real(10)
end

check_Set_Comp2Real_code_warntype()