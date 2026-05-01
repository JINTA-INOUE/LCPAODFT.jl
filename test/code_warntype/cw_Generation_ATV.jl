include("../../src/utils/Generation_ATV.jl")


function check_Generation_ATV_code_warntype()
    CpyCell = 3
    Latvecs = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int32, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    # @code_warntype Generation_ATV!(CpyCell, Latvecs, atv, atv_ijk, ratv)
    # @code_warntype Generation_ATV!(CpyCell, Latvecs, atv)
    # @code_warntype Generation_ATV_ijk!(CpyCell, atv_ijk)
    # @code_warntype Generation_RATV!(CpyCell, ratv)
end

check_Generation_ATV_code_warntype()