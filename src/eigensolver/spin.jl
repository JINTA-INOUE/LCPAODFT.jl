# S0 = U†U (upper) or LL† (lower). Reduce the two Hermitian diagonal
# blocks with HEGST and the authoritative off-diagonal block with TRSM.
# Views use the full spinor matrix's leading dimension, not the orbital size.
function _reduce!(ws::EigenWorkspace{ComplexF64}, factor::SpinOverlapFactor)
    A, U = ws.matrix, factor.matrix
    n = size(U,1)
    top, bottom = 1:n, n+1:2n
    _reduce_matrix!(ws, view(A,top,top), U)
    _reduce_matrix!(ws, view(A,bottom,bottom), U)
    if ws.uplo == 'U'
        B = view(A,top,bottom)
        LA.BLAS.trsm!('L','U','C','N',one(ComplexF64),U,B)
        LA.BLAS.trsm!('R','U','N','N',one(ComplexF64),U,B)
    else
        B = view(A,bottom,top)
        LA.BLAS.trsm!('L','L','N','N',one(ComplexF64),U,B)
        LA.BLAS.trsm!('R','L','C','N',one(ComplexF64),U,B)
    end
    return nothing
end

function _backtransform!(ws::EigenWorkspace{ComplexF64}, factor::SpinOverlapFactor)
    U = factor.matrix
    n = size(U,1)
    trans = ws.uplo == 'U' ? 'N' : 'C'
    for rows in (1:n, n+1:2n)
        LA.BLAS.trsm!('L',ws.uplo,trans,'N',one(ComplexF64),U,view(ws.matrix,rows,:))
    end
    return nothing
end
