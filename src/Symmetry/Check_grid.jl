_nint(x::Real) = round(Int, x, RoundNearestTiesAway)

function check_grid_compatibility(
    rotations::AbstractVector{<:AbstractMatrix{<:Integer}},
    translations::AbstractVector{<:AbstractVector{<:Real}},
    Ngrid::NTuple{3,<:Integer};
    check_translation::Bool=true,
    tolerance::Real=1e-5)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    operation_count = length(rotations)
    incompatible = zeros(Int32, operation_count)
    incompatible_count = 0

    for isym = 1:operation_count
        rotation = rotations[isym]    
        bad_rotation = any(mod(rotation[i,j]*Ngrid[j],Ngrid[i]) ≠ 0 for i = 1:3 for j = 1:3 if i ≠ j)

        bad_translation = false
        if check_translation
            for axis = 1:3
                scaled = translations[isym][axis] * Ngrid[axis]
                if abs(scaled - _nint(scaled)) / Ngrid[axis] > tolerance
                    bad_translation = true
                    break
                end
            end
        end

        if bad_rotation || bad_translation
            incompatible_count += 1
            incompatible[incompatible_count] = isym
        end
    end
    _trim_vectors!(incompatible_count, incompatible)

    return isempty(incompatible), incompatible
end


function _next_symmetry_compatible_fft_size(size::Integer, factor::Integer)

    size > 0 || throw(ArgumentError("FFT grid dimensions must be positive"))
    factor > 0 || throw(ArgumentError("FFT symmetry factors must be positive"))

    candidate = Int(size)
    while true
        if mod(candidate, factor) == 0
            remainder = candidate
            for prime = (2, 3, 5, 7)
                while mod(remainder, prime) == 0
                    remainder = div(remainder, prime)
                end
            end
            remainder == 1 && return candidate
        end
        candidate = Base.checked_add(candidate, 1)
    end
end


function _trim_vectors!(new_length::Int, vectors::Vector...)
    new_length >= 0 || throw(ArgumentError("new vector length must be nonnegative"))
    for vector = vectors
        new_length <= length(vector) || throw(ArgumentError("preallocated vectors may only be trimmed, not expanded"))
        resize!(vector, new_length)
    end
    return nothing
end


function _fft_factors_from_translations(translations::AbstractVector{<:AbstractVector{<:Real}}; tolerance::Real=1e-5)
    factors = ones(Int, 3)
    for translation = translations, axis = 1:3
        value = translation[axis]
        abs(value) <= tolerance && continue
        denominator = _nint(1.0 / abs(value))
        factors[axis] = lcm(factors[axis], denominator)
    end
    return Tuple(factors)
end


function symmetry_compatible_grid(symmetry, Ngrid::NTuple{3,<:Integer})

    fft_factors = _fft_factors_from_translations(symmetry.translations)
    compatible = ntuple(axis -> _next_symmetry_compatible_fft_size(Ngrid[axis], fft_factors[axis]), 3)
    grid_ok, incompatible = check_grid_compatibility(symmetry.rotations, symmetry.translations, compatible)
    grid_ok || throw(ArgumentError(
        "FFT grid $compatible is not closed under symmetry operations $incompatible; " *
        "choose Ngrid dimensions compatible with rotations and fractional translations"))

    return compatible
end
