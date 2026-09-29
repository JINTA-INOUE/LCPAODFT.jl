"""
    EigenWorkspace(T, n; uplo=:U, workspace=:optimized, store_overlap=true)

Reusable dense all-eigenpair workspace for Float64 or ComplexF64 (n ≥ 1).
Owns LAPACK scratch arrays and optionally a Cholesky factor of S. Create one per
concurrent solve; it is not thread-safe. Choose the BLAS backend before creation.
Use prepare_overlap! again when geometry, basis or k point changes.
`uplo` selects the internal LAPACK triangle (:U or :L), independently of input
wrappers. `workspace=:optimized` reserves scratch for blocked eigenvector recovery;
`:lapack` retains the native eigensolver query for controlled comparisons.
With `store_overlap=false`, no internal N×N overlap buffer is allocated. Pass an
OverlapFactor to eigen!(ws,H,factor), reusing one workspace across local k points.
"""
mutable struct EigenWorkspace{T<:SupportedScalar}
    matrix::Matrix{T}
    overlap_factor::Matrix{T}
    values::Vector{Float64}
    work::Vector{T}
    rwork::Vector{Float64}
    iwork::Vector{BlasInt}
    info::Base.RefValue{BlasInt}
    overlap_ready::Bool
    backend::LA.BLAS.LBTConfig
    uplo::Char
end

function EigenWorkspace(::Type{T}, n::Integer; uplo::Symbol=:U,
                        workspace::Symbol=:optimized,
                        store_overlap::Bool=true) where {T<:SupportedScalar}
    n >= 1 || throw(ArgumentError("matrix dimension must be positive"))
    triangle = _triangle_char(uplo)
    workspace in (:optimized, :lapack) || throw(ArgumentError("workspace must be :optimized or :lapack"))
    overlap_storage = store_overlap ? zeros(T,n,n) : Matrix{T}(undef,0,0)
    ws = EigenWorkspace(zeros(T, n, n), overlap_storage, zeros(n),
                        zeros(T, 1), zeros(1), zeros(BlasInt, 1),
                        Ref{BlasInt}(0), false, LA.BLAS.get_config(),
                        triangle)
    _query_workspace!(ws, workspace)
    return ws
end

function _check_backend(ws::EigenWorkspace)
    # get_config is cached by LBT; changing forwards invalidates that snapshot.
    LA.BLAS.get_config() === ws.backend || throw(ArgumentError("BLAS/LAPACK configuration changed; create a new EigenWorkspace"))
end

function _check_matrix(A::AbstractMatrix, n::Integer)
    Base.require_one_based_indexing(A)
    size(A) == (n, n) || throw(DimensionMismatch("expected a $n × $n matrix"))
    all(isfinite, A) || throw(ArgumentError("matrix contains non-finite values"))
    LA.ishermitian(A) || throw(ArgumentError("matrix must be symmetric/Hermitian; wrap an authoritative triangle explicitly if needed"))
    return nothing
end

function _check_no_alias(ws::EigenWorkspace, A)
    for storage in (ws.matrix,ws.overlap_factor,ws.values,ws.work,ws.rwork,ws.iwork)
        Base.mightalias(A,storage) && throw(ArgumentError("input must not alias workspace storage"))
    end
end

"""
    prepare_overlap!(ws, S)

Copy and factor S without modifying the input. S must be Hermitian positive
definite. Failure invalidates any previous cached factor; no silent regularization
or basis truncation is performed. Cached S is not automatically synchronized with
later edits to the caller's S: call this function explicitly after any change.
"""
function prepare_overlap!(ws::EigenWorkspace{T}, S::AbstractMatrix{T}) where T
    ws.overlap_ready = false
    _check_backend(ws)
    isempty(ws.overlap_factor) && throw(ArgumentError("workspace has store_overlap=false; prepare an OverlapFactor and pass it to eigen!"))
    _check_matrix(S, length(ws.values))
    _check_no_alias(ws, S)
    copyto!(ws.overlap_factor, S)
    _, info = LA.LAPACK.potrf!(ws.uplo, ws.overlap_factor)
    info == 0 || throw(LA.PosDefException(info))
    ws.overlap_ready = true
    return ws
end
