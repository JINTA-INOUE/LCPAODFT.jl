# Native calls are isolated here. Link through Julia's libblastrampoline, using
# its BlasInt and symbol convention. Never load or replace a BLAS backend here.
# LAPACK routines: xSYGST/xHEGST and xSYEVD/xHEEVD.

function _check_info(info::Integer)
    info < 0 && throw(ArgumentError("LAPACK rejected argument $(-info)"))
    info > 0 && throw(LA.LAPACKException(info))
    return nothing
end

for (T, routine) in ((Float64, :dsygst_), (ComplexF64, :zhegst_))
    @eval function _reduce_matrix!(ws::EigenWorkspace{$T}, A::StridedMatrix{$T}, U::Matrix{$T})
        n = size(A, 1)
        ccall((@blasfunc($routine), liblbt), Cvoid,
              (Ref{BlasInt}, Ref{UInt8}, Ref{BlasInt}, Ptr{$T}, Ref{BlasInt},
               Ptr{$T}, Ref{BlasInt}, Ref{BlasInt}, Clong),
              1, ws.uplo, n, A, stride(A,2), U, stride(U,2), ws.info, 1)
        _check_info(ws.info[])
        return nothing
    end
end

_reduce!(ws::EigenWorkspace, U::Matrix) = _reduce_matrix!(ws,ws.matrix,U)

function _evd!(ws::EigenWorkspace{Float64}, query::Bool=false, vectors::Bool=true)
    n = length(ws.values)
    lw = query ? -1 : length(ws.work)
    liw = query ? -1 : length(ws.iwork)
    ccall((@blasfunc(dsyevd_), liblbt), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
           Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt},
           Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
          vectors ? 'V' : 'N', ws.uplo, n, ws.matrix, n, ws.values, ws.work, lw,
          ws.iwork, liw, ws.info, 1, 1)
    _check_info(ws.info[])
    return nothing
end

function _evd!(ws::EigenWorkspace{ComplexF64}, query::Bool=false, vectors::Bool=true)
    n = length(ws.values)
    lw = query ? -1 : length(ws.work)
    lrw = query ? -1 : length(ws.rwork)
    liw = query ? -1 : length(ws.iwork)
    ccall((@blasfunc(zheevd_), liblbt), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
           Ptr{Float64}, Ptr{ComplexF64}, Ref{BlasInt}, Ptr{Float64},
           Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
          vectors ? 'V' : 'N', ws.uplo, n, ws.matrix, n, ws.values, ws.work, lw,
          ws.rwork, lrw, ws.iwork, liw, ws.info, 1, 1)
    _check_info(ws.info[])
    return nothing
end

# ZHEEVD's native query can leave only N elements for ZUNMTR after the
# N-element TAU and N²-element eigenvector regions. Query the recovery routine
# as well as its QR/QL kernel: some LAPACK versions need extra block storage in
# that kernel beyond the N*NB advertised by ZUNMTR. Query calls do not read A/TAU.
function _recovery_lwork(ws::EigenWorkspace{ComplexF64})
    n = length(ws.values)
    n <= 1 && return 1
    answer = zeros(ComplexF64, 1)
    ccall((@blasfunc(zunmtr_), liblbt), Cvoid,
          (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
           Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
           Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong, Clong),
          'L', ws.uplo, 'N', n, n, ws.matrix, n, ws.work, ws.matrix, n,
          answer, -1, ws.info, 1, 1, 1)
    _check_info(ws.info[])
    result = ceil(Int, real(answer[1]))
    result = max(result, _recovery_kernel_lwork(ws, answer))
    return result
end

for (triangle, routine) in (('U', :zunmql_), ('L', :zunmqr_))
    @eval function _recovery_kernel_lwork(ws::EigenWorkspace{ComplexF64}, answer, ::Val{$(QuoteNode(triangle))})
        n = length(ws.values)
        ccall((@blasfunc($routine), liblbt), Cvoid,
              (Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
               Ptr{ComplexF64}, Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
               Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt}, Ref{BlasInt}, Clong, Clong),
              'L', 'N', n-1, n, n-1, ws.matrix, n, ws.work, ws.matrix, n,
              answer, -1, ws.info, 1, 1)
        _check_info(ws.info[])
        return ceil(Int, real(answer[1]))
    end
end
_recovery_kernel_lwork(ws, answer) = _recovery_kernel_lwork(ws, answer, Val(ws.uplo))

function _query_workspace!(ws::EigenWorkspace{T}, workspace::Symbol) where T
    _evd!(ws, true)
    lw = max(1,ceil(Int,real(ws.work[1])))
    liw = max(1,Int(ws.iwork[1]))
    lrw = T === ComplexF64 ? max(1,ceil(Int,ws.rwork[1])) : 1
    if T === ComplexF64
        n = length(ws.values)
        if workspace === :optimized && n > 1
            # EVD stores N² vectors in WORK after the N-element TAU region.
            offset = n*n+n
            lw = max(lw,offset+_recovery_lwork(ws))
        end
    end
    # Query all sizes before allocating. A second resize! would copy the large
    # WORK buffer and reserve excess Vector capacity (not merely the extra LWORK).
    ws.work = Vector{T}(undef,lw)
    ws.rwork = Vector{Float64}(undef,lrw)
    ws.iwork = Vector{BlasInt}(undef,liw)
    return ws
end
