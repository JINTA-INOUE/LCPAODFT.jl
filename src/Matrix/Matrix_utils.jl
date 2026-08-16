function _packed_allreduce_matrices!(arrays, comm)
    MPI.Comm_size(comm) == 1 && return arrays
    total_length = sum(length, arrays; init=0)
    iszero(total_length) && return arrays
    packed = Vector{Float64}(undef, total_length)
    offset = 1
    @inbounds for array in arrays
        copyto!(packed, offset, array, 1, length(array))
        offset += length(array)
    end
    MPI.Allreduce!(packed, MPI.SUM, comm)
    offset = 1
    @inbounds for array in arrays
        copyto!(array, 1, packed, offset, length(array))
        offset += length(array)
    end
    return arrays
end