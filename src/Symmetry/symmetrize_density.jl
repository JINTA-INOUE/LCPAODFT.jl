struct DensitySymmetryPlan{Index<:Integer,Forward,Backward}
    Ngrid::NTuple{3,Int}
    Grid_xyz::Array{ComplexF64,3}
    Grid_xyz_sym::Array{ComplexF64,3}
    forward::Forward
    backward::Backward
    targets::Matrix{Index}
    phase1::Matrix{ComplexF64}
    phase2::Matrix{ComplexF64}
    phase3::Matrix{ComplexF64}
    scale::Float64
end

function DensitySymmetryPlan(
    Ngrid::NTuple{3,<:Integer},
    rotations::AbstractVector{<:AbstractMatrix{<:Integer}},
    translations::AbstractVector{<:AbstractVector{<:Real}},
    origin::AbstractVector{<:Real};
    flags=FFTW.ESTIMATE,
    tolerance::Real=1e-8,
)
    Base.require_one_based_indexing(rotations, translations, origin)
    grid = (Int(Ngrid[1]), Int(Ngrid[2]), Int(Ngrid[3]))
    all(>(0), grid) || throw(ArgumentError("grid dimensions must be positive"))
    ng = Base.checked_mul(Base.checked_mul(grid[1], grid[2]), grid[3])
    nops = length(rotations)
    nops > 0 || throw(ArgumentError("at least one symmetry operation is required"))
    length(translations) == nops || throw(DimensionMismatch("rotations and translations must have the same length"))
    length(origin) == 3 || throw(DimensionMismatch("fractional origin must have three components"))
    all(isfinite, origin) || throw(ArgumentError("origin must be finite"))
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError("tolerance must be finite and positive"))

    # Finish validating the small inputs before allocating the large index table.
    effective = Matrix{Float64}(undef, 3, nops)
    for s = 1:nops
        rotation, translation = rotations[s], translations[s]
        Base.require_one_based_indexing(rotation, translation)
        size(rotation) == (3, 3) || throw(DimensionMismatch("each rotation must be 3 × 3"))
        length(translation) == 3 || throw(DimensionMismatch("each translation must have three components"))
        all(isfinite, translation) || throw(ArgumentError("translations must be finite"))
        for j = 1:3
            rotated_origin = 0.0
            for i = 1:3
                rotated_origin += Float64(origin[i]) * rotation[i, j]
                # h = S*g must be well defined modulo the FFT grid dimensions.
                scaled_rotation = Base.checked_mul(Int(rotation[j, i]), grid[i])
                mod(scaled_rotation, grid[j]) == 0 || throw(ArgumentError(
                    "grid $grid is incompatible with rotation $s",
                ))
            end
            effective[j, s] = Float64(translation[j]) + Float64(origin[j]) - rotated_origin
            scaled = effective[j, s] * grid[j]
            isfinite(scaled) || throw(ArgumentError("effective translation must be finite"))
            abs(scaled - round(scaled)) <= tolerance || throw(ArgumentError(
                "grid $grid with origin $origin is not mapped onto itself by operation $s",
            ))
        end
    end

    Index = ng <= typemax(Int32) ? Int32 : Int
    targets = Matrix{Index}(undef, ng, nops)
    phase1 = Matrix{ComplexF64}(undef, grid[1], nops)
    phase2 = Matrix{ComplexF64}(undef, grid[2], nops)
    phase3 = Matrix{ComplexF64}(undef, grid[3], nops)
    phases = (phase1, phase2, phase3)
    seen = falses(ng)
    linear = LinearIndices(grid)
    for s = 1:nops
        fill!(seen, false)
        for axis = 1:3, i = 1:grid[axis]
            g = _signed_fft_index(i, grid[axis])
            phases[axis][i, s] = cis(-2pi * effective[axis, s] * g)
        end
        rotation = rotations[s]
        for (ig, xyz) = enumerate(CartesianIndices(grid))
            target, _, _, _ = reciprocal_grid_target(rotation, xyz, grid)
            target_index = linear[target]
            seen[target_index] && throw(ArgumentError(
                "operation $s does not produce a one-to-one FFT grid map",
            ))
            seen[target_index] = true
            targets[ig, s] = target_index
        end
    end

    Grid_xyz = zeros(ComplexF64, grid)
    Grid_xyz_sym = zeros(ComplexF64, grid)
    # MEASURE/PATIENT may overwrite these scratch arrays during planning.
    forward = FFTW.plan_fft!(Grid_xyz; flags=flags)
    backward = FFTW.plan_ifft!(Grid_xyz_sym; flags=flags)
    return DensitySymmetryPlan(
        grid, Grid_xyz, Grid_xyz_sym, forward, backward, targets,
        phase1, phase2, phase3, 1.0 / nops,
    )

end

function symmetrize_Grid!(
    Ngrid, 
    Rots_direct, translations,
    Grid_Origin,
    Grid_xyz::AbstractArray{<:ComplexF64,3},
    Grid_xyz_sym::AbstractArray{<:ComplexF64,3},
    Grid_Vec::Vector{Float64})

    Set_Density2Density_xyz!(Grid_Vec, Grid_xyz)
    symmetrize_density_reciprocal!(Ngrid, Rots_direct, translations, Grid_Origin, Grid_xyz, Grid_xyz_sym)
    Set_Density_xyz2Density!(Grid_xyz_sym, Grid_Vec)
end


function symmetrize_density_reciprocal!(Ngrid, Rots_direct, translations, Grid_Origin, Grid_xyz::AbstractArray{<:ComplexF64,3}, Grid_xyz_sym::AbstractArray{<:ComplexF64,3})
    
    n_operations = length(Rots_direct)
    
    fill!(Grid_xyz_sym, 0.0)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    effective = zeros(Float64, 3)
    in1 = zeros(ComplexF64, Ngrid1)
    in2 = zeros(ComplexF64, Ngrid2)
    in3 = zeros(ComplexF64, Ngrid3)
    FFT_Grid!(Ngrid, in1, in2, in3, Grid_xyz)

    for is = 1:n_operations
        Rot_direct = Rots_direct[is]
        real_grid_target!(Grid_Origin, Rot_direct, translations[is], effective)
        @inbounds for xyz in CartesianIndices(Ngrid)
            xyz_id_rot, freq1, freq2, freq3 = reciprocal_grid_target(Rot_direct, xyz, Ngrid)
            phase = effective[1]*freq1 + effective[2]*freq2 + effective[3]*freq3
            Grid_xyz_sym[xyz_id_rot] += Grid_xyz[xyz]*cis(-2*pi*phase)
        end
    end
    @. Grid_xyz_sym = Grid_xyz_sym/n_operations

    iFFT_Grid!(Ngrid, in1, in2, in3, Grid_xyz_sym)
end


@inline function _signed_fft_index(index::Int, size::Int)
    zero_based = index - 1
    if zero_based <= fld(size, 2)
        return zero_based
    else
        return zero_based - size
    end
end


@inline function real_grid_target!(Grid_Origin, Rot_direct, translation, effective)
    for j = 1:3
        rotated_origin = 0.0
        for i = 1:3
            rotated_origin += Grid_Origin[i]*Rot_direct[i,j]
        end

        effective[j] = translation[j] + Grid_Origin[j] - rotated_origin
    end
end


@inline function reciprocal_grid_target(Rot_direct, xyz, Ngrid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    
    freq1 = _signed_fft_index(xyz[1], Ngrid1)
    freq2 = _signed_fft_index(xyz[2], Ngrid2)
    freq3 = _signed_fft_index(xyz[3], Ngrid3)

    r1 = Rot_direct[1,1]*freq1 + Rot_direct[1,2]*freq2 + Rot_direct[1,3]*freq3
    r2 = Rot_direct[2,1]*freq1 + Rot_direct[2,2]*freq2 + Rot_direct[2,3]*freq3
    r3 = Rot_direct[3,1]*freq1 + Rot_direct[3,2]*freq2 + Rot_direct[3,3]*freq3

    i_rot = mod(r1, Ngrid1)+1
    j_rot = mod(r2, Ngrid2)+1
    k_rot = mod(r3, Ngrid3)+1


    return CartesianIndex(i_rot,j_rot,k_rot), freq1, freq2, freq3
end


function Set_Density2Density_xyz!(Grid_Vec::Vector{Float64}, Grid_xyz::AbstractArray{<:ComplexF64,3})
    Ngrid1, Ngrid2, Ngrid3 = size(Grid_xyz)
    @inbounds for i = 1:Ngrid1, j = 1:Ngrid2, k = 1:Ngrid3
        idx = (i-1)*Ngrid2*Ngrid3 + (j-1)*Ngrid3 + k
        Grid_xyz[i,j,k] = Grid_Vec[idx]
    end
end


function Set_Density_xyz2Density!(Grid_xyz::AbstractArray{<:ComplexF64,3}, Grid_Vec::Vector{Float64})
    Ngrid1, Ngrid2, Ngrid3 = size(Grid_xyz)
    counts = 0
    @inbounds for i = 1:Ngrid1, j = 1:Ngrid2, k = 1:Ngrid3
        counts += 1
        Grid_Vec[counts] = real(Grid_xyz[i,j,k])
    end
end