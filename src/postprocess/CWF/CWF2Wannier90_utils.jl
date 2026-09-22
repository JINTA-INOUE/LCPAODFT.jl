"""Return the absolute Wannier90 seed path while preserving a user-supplied directory."""
_cwf_seed_path(filename::AbstractString) = abspath(filename)

"""Return the private directory used for rank-local CWF2MLWF intermediates."""
_cwf_work_dir(filename::AbstractString) = _cwf_seed_path(filename)*"_work_cwf"

function _cwf_work_file(filename::AbstractString, suffix::AbstractString)
    seed = _cwf_seed_path(filename)
    return joinpath(_cwf_work_dir(filename), basename(seed)*"_"*suffix*".jld2")
end

function Write_Variable_work_file(filename::AbstractString, myrank, variable, name::AbstractString)
    work_file = _cwf_work_file(filename, "$(name)$(myrank)")
    jldopen(work_file, "w") do file
        file[name] = variable
    end
end
