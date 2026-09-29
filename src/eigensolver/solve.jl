"""
    eigen!(ws, H; generalized=true)

Solve H*C = S*C*Λ using the S previously supplied to prepare_overlap!, or
H*C = C*Λ with generalized=false. H is copied and is not modified.
Returns LinearAlgebra.GeneralizedEigen (or Eigen for a standard problem), with
ascending values and column eigenvectors.

The returned arrays ALIAS workspace storage: a later solve overwrites them.
Copy F.values and F.vectors when retaining a result or warm-start subspace.
No BLAS backend/thread settings or LinearAlgebra methods are changed.
"""
function eigen!(ws::EigenWorkspace{T}, H::AbstractMatrix{T}; generalized::Bool=true) where T
    _check_backend(ws)
    _check_matrix(H, length(ws.values))
    _check_no_alias(ws, H)
    if generalized && !ws.overlap_ready
        isempty(ws.overlap_factor) && throw(ArgumentError("workspace has store_overlap=false; use eigen!(ws, H, factor::OverlapFactor)"))
        throw(ArgumentError("call prepare_overlap!(ws, S) before a generalized solve"))
    end
    return _solve!(ws,H,generalized ? ws.overlap_factor : nothing)
end

"""
    eigen!(ws, H, factor::OverlapFactor)

Solve H*C = S(k)*C*Λ with a separately cached factor. The factor is read-only
during this call; it is neither copied nor retained by ws. Reuse one ws for the
local k-point loop, and copy each result into its destination before the next call.
The workspace and factor must have matching element types and uplo. For an
OverlapFactor their dimensions match; for SpinOverlapFactor(S0) the workspace
and H have twice the orbital dimension. H retains the full spin coupling.
"""
eigen!(ws::EigenWorkspace{T}, H::AbstractMatrix{T}, factor::AbstractOverlapFactor{T}) where T = _solve_external!(ws, H, factor, true)

"""
    eigvals!(ws, H, factor::OverlapFactor)

Eigenvalues-only generalized solve using a cached S factor and existing scratch.
The result aliases ws.values; H and factor are unchanged. No eigenvectors or
triangular backtransform are computed. A subsequent eigen!/eigvals! overwrites it.
"""
eigvals!(ws::EigenWorkspace{T}, H::AbstractMatrix{T}, factor::OverlapFactor{T}) where T = _solve_external!(ws, H, factor, false)

function _solve_external!(ws::EigenWorkspace{T}, H::AbstractMatrix{T}, factor::AbstractOverlapFactor{T}, vectors::Bool) where T
    _check_backend(ws)
    _check_backend(factor)
    factor.ready || throw(ArgumentError("call prepare_overlap!(factor, S) before solving"))
    n = length(ws.values)
    (_problem_size(factor) == n && size(factor.matrix,1) == size(factor.matrix,2)) ||
        throw(DimensionMismatch("overlap factor dimension does not match workspace"))
    factor.uplo == ws.uplo || throw(ArgumentError("overlap factor and workspace must use the same uplo"))
    _check_no_alias(ws,factor.matrix)
    Base.mightalias(H,factor.matrix) && throw(ArgumentError("H must not alias cached factor storage"))
    _check_matrix(H,n)
    _check_no_alias(ws,H)
    return _solve!(ws,H,_overlap_storage(factor); vectors)
end

function _backtransform!(ws::EigenWorkspace, overlap::Matrix)
    factor = ws.uplo == 'U' ? LA.UpperTriangular(overlap) : adjoint(LA.LowerTriangular(overlap))
    LA.ldiv!(factor, ws.matrix)
    return nothing
end

# Standard, dense-overlap and spin-overlap problems share the same solve steps.
# Only standardization and backtransformation depend on the overlap storage.
function _solve!(ws::EigenWorkspace, H, overlap; vectors::Bool=true)
    copyto!(ws.matrix, H)
    if overlap !== nothing
        _reduce!(ws,overlap)
    end
    _evd!(ws, false, vectors)
    vectors || return ws.values
    if overlap !== nothing
        _backtransform!(ws, overlap)
    end
    return overlap === nothing ? LA.Eigen(ws.values, ws.matrix) : LA.GeneralizedEigen(ws.values, ws.matrix)
end

"""
    LCPAODEigen.eigen(H[, S])

Non-destructive convenience interface for Float64 / ComplexF64 dense Hermitian
problems. Creates an owned workspace and uses the dedicated solver.
For SCF reuse, create EigenWorkspace and call prepare_overlap! / eigen! instead.
This is a separate function, not an extension of LinearAlgebra.eigen.
"""
function eigen(H::AbstractMatrix{T}, S::AbstractMatrix{T}) where {T<:SupportedScalar}
    n = size(H, 1)
    n >= 1 || throw(ArgumentError("empty matrices are not supported"))
    _check_matrix(H, n)
    ws = EigenWorkspace(T, n)
    prepare_overlap!(ws, S)
    return _solve!(ws, H, ws.overlap_factor)
end

function eigen(H::AbstractMatrix{T}) where {T<:SupportedScalar}
    n = size(H, 1)
    n >= 1 || throw(ArgumentError("empty matrices are not supported"))
    _check_matrix(H, n)
    ws = EigenWorkspace(T, n;store_overlap=false)
    return _solve!(ws, H, nothing)
end

"""
    diagnostics(H, S, F)
    diagnostics(H, F)

Return normalized Frobenius residual and S-orthogonality (ordinary orthogonality
without S). Allocates diagnostic matrices; use for verification, not every hot loop.
"""
function diagnostics(H, S, F::Union{LA.Eigen,LA.GeneralizedEigen})
    C, w = F.vectors, F.values
    SC = isnothing(S) ? C : S*C
    nr = LA.norm(H*C - SC*LA.Diagonal(w))
    ns = isnothing(S) ? 1.0 : LA.norm(S)
    denominator = (LA.norm(H) + maximum(abs, w)*ns)*LA.norm(C)
    residual = denominator == 0 ? (nr == 0 ? 0.0 : Inf) : nr/denominator
    orthogonality = LA.norm(C'*SC - LA.I)/sqrt(length(w))
    return (; residual, orthogonality)
end
diagnostics(H, F::LA.Eigen) = diagnostics(H, nothing, F)
