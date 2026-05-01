include("../../src/KPoints.jl")

function check_KPoints_code_warntype()
    kmesh = (3,3,3)
    Shift_K_Point = 1e-8
    @code_warntype _Gen_KPoints_Gcenter_TRS(kmesh, Shift_K_Point)
    @code_warntype _Gen_KPoints_Gcenter_noTRS(kmesh, Shift_K_Point)
    @code_warntype _Gen_KPoints_MP_TRS(kmesh, Shift_K_Point)
    @code_warntype _Gen_KPoints_MP_noTRS(kmesh, Shift_K_Point)
end

check_KPoints_code_warntype()