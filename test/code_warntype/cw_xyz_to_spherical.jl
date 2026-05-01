include("../../src/utils/xyz_to_spherical.jl")


function check_xyz_to_spherical_code_warntype()
    x = rand(Float64)
    y = rand(Float64)
    z = rand(Float64)
    tau = [x, y, z]
    @code_warntype xyz_to_spherical(x, y, z)
    @code_warntype xyz_to_spherical(tau)
end


check_xyz_to_spherical_code_warntype()