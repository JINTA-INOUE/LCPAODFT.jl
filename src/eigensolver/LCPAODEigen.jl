module LCPAODEigen

# Do not import or extend LinearAlgebra.eigen/eigen!.
import LinearAlgebra as LA
using LinearAlgebra.BLAS: @blasfunc

export EigenWorkspace
export OverlapFactor
export SpinOverlapFactor
export prepare_overlap!
export eigen
export eigen!
export eigvals!
export diagnostics

const SupportedScalar = Union{Float64,ComplexF64}
const BlasInt = LA.BlasInt
const liblbt = LA.BLAS.libblastrampoline

include("overlap.jl")
include("workspace.jl")
include("lapack.jl")
include("solve.jl")
include("spin.jl")

end
