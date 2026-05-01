include("../../src/utils/Generation_ATV.jl")
include("../../src/utils/Find_CGrid.jl")


function check_Find_ln_code_warntype()
    Ngrid = (3,3,3)
    @code_warntype Find_ln(Int32(Ngrid[1]), 1)
    @code_warntype Find_ln(Int64(Ngrid[1]), 1)
end


function check_Find_CGrid_code_warntype()
    CpyCell = 3
    Latvecs = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]
    Ngrid = Int32.((3,3,3))

    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int32, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid[1]
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid[2]
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid[3]

    Cxyz = zeros(Float64, 3)
    NOC = zeros(Int32, 4)
    
    Grid_Origin = [0.0, 0.0, 0.0]

    Generation_ATV!(CpyCell, Latvecs, atv, atv_ijk, ratv)
    @code_warntype Find_CGrids!(NOC, Cxyz, CpyCell, Ngrid, 1, 2, 3, atv, ratv, gLatvecs, Grid_Origin)
end

# check_Find_ln_code_warntype()
check_Find_CGrid_code_warntype()