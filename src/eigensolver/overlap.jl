abstract type AbstractOverlapFactor{T<:SupportedScalar} end

_triangle_char(uplo::Symbol) = uplo === :U ? 'U' : uplo === :L ? 'L' : throw(ArgumentError("uplo must be :U or :L"))

"""
    OverlapFactor(T, n; uplo=:U)
    OverlapFactor(S; uplo=:U)

Own only the Cholesky factor of one S(k), without eigensolver scratch arrays.
The (T,n) constructor requires prepare_overlap! before use; the S constructor
copies and factors S. Geometry/basis/k-point changes require preparing it again.
Use one factor per local k point and one EigenWorkspace per concurrent solve.
Solves only read this cache. Do not modify/reprepare it during a solve, or send
workspace/cache objects between MPI processes: create them in their owning rank.
"""
mutable struct OverlapFactor{T<:SupportedScalar} <: AbstractOverlapFactor{T}
    matrix::Matrix{T}
    uplo::Char
    ready::Bool
    backend::LA.BLAS.LBTConfig
end

function OverlapFactor(::Type{T}, n::Integer; uplo::Symbol=:U) where {T<:SupportedScalar}
    n >= 1 || throw(ArgumentError("matrix dimension must be positive"))
    triangle = _triangle_char(uplo)
    return OverlapFactor(Matrix{T}(undef,n,n), triangle, false, LA.BLAS.get_config())
end

function OverlapFactor(S::AbstractMatrix{T}; uplo::Symbol=:U) where {T<:SupportedScalar}
    return prepare_overlap!(OverlapFactor(T,size(S,1);uplo=uplo), S)
end

function _check_backend(factor::AbstractOverlapFactor)
    LA.BLAS.get_config() === factor.backend || throw(ArgumentError("BLAS/LAPACK configuration changed; create a new overlap factor"))
end

"""
    prepare_overlap!(factor::OverlapFactor, S)
    prepare_overlap!(factor::SpinOverlapFactor, S0)

Replace this k point's cached Cholesky factor, preserving the input S. Failure
invalidates the previous factor. The same factor storage is reused.
For SpinOverlapFactor supply only the orbital overlap S0, not diag(S0,S0).
"""
function prepare_overlap!(factor::AbstractOverlapFactor{T}, S::AbstractMatrix{T}) where T
    factor.ready = false
    _check_backend(factor)
    Base.mightalias(factor.matrix,S) && throw(ArgumentError("input must not alias cached factor storage"))
    _check_matrix(S,size(factor.matrix,1))
    copyto!(factor.matrix,S)
    _, info = LA.LAPACK.potrf!(factor.uplo,factor.matrix)
    info == 0 || throw(LA.PosDefException(info))
    factor.ready = true
    return factor
end

"""
    SpinOverlapFactor(S0; uplo=:U)
    SpinOverlapFactor(n; uplo=:U)

Complex noncollinear overlap S = diag(S0,S0), stored as ONE n×n Cholesky
factor. S0 must be ComplexF64 Hermitian positive definite. The integer constructor
is unprepared; call prepare_overlap!(factor,S0) before use. Pair with a 2n×2n
EigenWorkspace and eigen!(ws,H,factor). H remains a coupled complex spinor problem.
The factor's matrix owns n² complex elements; no expanded overlap is retained.
Geometry/basis/k-point changes require re-preparation. Backend and ownership rules
are the same as for OverlapFactor.
"""
mutable struct SpinOverlapFactor <: AbstractOverlapFactor{ComplexF64}
    matrix::Matrix{ComplexF64}
    uplo::Char
    ready::Bool
    backend::LA.BLAS.LBTConfig
end

function SpinOverlapFactor(n::Integer; uplo::Symbol=:U)
    n >= 1 || throw(ArgumentError("orbital dimension must be positive"))
    return SpinOverlapFactor(Matrix{ComplexF64}(undef,n,n), _triangle_char(uplo), false, LA.BLAS.get_config())
end
SpinOverlapFactor(S0::AbstractMatrix{ComplexF64}; uplo::Symbol=:U) = prepare_overlap!(SpinOverlapFactor(size(S0,1); uplo), S0)

_problem_size(factor::OverlapFactor) = size(factor.matrix,1)
_problem_size(factor::SpinOverlapFactor) = 2size(factor.matrix,1)
_overlap_storage(factor::OverlapFactor) = factor.matrix
_overlap_storage(factor::SpinOverlapFactor) = factor
