# Rank-local state, owned by a single KSsolve_SCF! invocation. Only the overlap
# factors scale with local k-point count; LAPACK scratch is shared across spins/k.
struct SCFEigenCache{T,F<:LCPAODEigen.AbstractOverlapFactor{T}}
    workspace::Union{Nothing,LCPAODEigen.EigenWorkspace{T}}
    factors::Vector{F}
    buffers::Vector{Matrix{T}}
end

_scf_structured(H::AbstractMatrix{Float64}) = Symmetric(H, :U)
_scf_structured(H::AbstractMatrix{ComplexF64}) = Hermitian(H, :U)
_scf_matrix_args(g) = (g.Natom, g.Total_NumOrbs, g.MP, g.FNAN, g.natn)
_scf_matrix_args(g, k) = (_scf_matrix_args(g)..., g.ncn, g.atv_ijk, k)

"""
    prepare_scf_eigensolver(OLP, grid, electron, kpoints; workspace, cal_force)

Prepare one geometry's rank-local S factors before the SCF loop. Empty crystal
ranks own no dense storage. Call again after changing geometry, basis or k mesh.
"""
function prepare_scf_eigensolver(OLP, grid, electron, kpoints; workspace=:optimized, cal_force=false)

    workspace in (:optimized, :lapack) || throw(ArgumentError("eigen_workspace must be :optimized or :lapack"))
    crystal = electron isa CrystalBloch
    nc = electron.SpinPol == "nc"
    T = crystal || nc ? ComplexF64 : Float64
    nk = crystal ? kpoints.MPI_Nkpt : 1
    n = Int(electron.Nfsize)
    block_spin = nc && electron.cal_mode == 1
    F = block_spin ? LCPAODEigen.SpinOverlapFactor : LCPAODEigen.OverlapFactor{T}
    factors = F[]
    buffers = Matrix{T}[]
    nk == 0 && return SCFEigenCache{T,F}(nothing, factors, buffers)
    ws = LCPAODEigen.EigenWorkspace(T, n; workspace, store_overlap=false)
    for ik in 1:nk
        args = crystal ? _scf_matrix_args(grid, kpoints.MPI_kpts[ik]) : _scf_matrix_args(grid)
        S = ws.matrix # Temporary assembly storage, not retained by factor.
        if block_spin
            # Assemble only S0, using existing rank scratch as a temporary view.
            S0 = view(S, 1:n÷2, 1:n÷2)
            HS_matrix!(S0, OLP, args...)
            push!(factors, LCPAODEigen.SpinOverlapFactor(Hermitian(S0, :U)))
            continue
        elseif nc
            m = n ÷ 2
            fill!(S, 0)
            HS_matrix!(view(S, 1:m, 1:m), OLP, args...)
            copyto!(view(S, m+1:n, m+1:n), view(S, 1:m, 1:m))
        else
            HS_matrix!(S, OLP, args...)
        end
        push!(factors, LCPAODEigen.OverlapFactor(_scf_structured(S)))
    end
    # Mode 1 uses Cnk as H input. Mode 2 retains only rank-local H/DM/EDM
    # buffers, reused across SCF iterations. Cluster has one reusable H buffer.
    nbuffers = crystal ? (electron.cal_mode == 2 ? 2 + Int(cal_force) : 0) : 1
    for _ in 1:nbuffers
        push!(buffers, zeros(T, n, n))
    end
    return SCFEigenCache{T,F}(ws, factors, buffers)
end

_scf_buffer(cache::SCFEigenCache{T}, i) where T = cache.workspace === nothing ? Matrix{T}(undef, 0, 0) : cache.buffers[i]

function _scf_eigen!(cache, H, ik)
    return LCPAODEigen.eigen!(cache.workspace, _scf_structured(H), cache.factors[ik])
end
function _scf_eigvals!(cache, H, ik)
    return LCPAODEigen.eigvals!(cache.workspace, _scf_structured(H), cache.factors[ik])
end

function _crystal_eigen_cached!(cache, Hks, iHks, electron, kpoints, grid)
    fill!(electron.Enk, 0)
    offset = kpoints.MPkpts[MPI.Comm_rank(MPI.COMM_WORLD)+1]
    for ik in 1:kpoints.MPI_Nkpt, spin in 1:electron.spinsize
        H = electron.Cnk[spin][ik]
        args = _scf_matrix_args(grid, kpoints.MPI_kpts[ik])
        if electron.SpinPol == "nc"
            HS_matrix_NC!(H, Hks, iHks, args...)
        else
            HS_matrix!(H, Hks[spin], args...)
        end
        F = _scf_eigen!(cache, H, ik)
        copyto!(view(electron.Enk, :, offset+ik, spin), F.values)
        # The next solve overwrites workspace vectors; DM needs every local Cnk.
        copyto!(H, F.vectors)
    end
    MPI.Allreduce!(electron.Enk, MPI.SUM, MPI.COMM_WORLD)
    return nothing
end

function _cluster_eigen_cached!(cache, Hks, iHks, electron, grid)
    H = cache.buffers[1]
    args = _scf_matrix_args(grid)
    for spin in 1:electron.spinsize
        if electron.SpinPol == "nc"
            HS_matrix_NC!(H, Hks, iHks, args...)
        else
            HS_matrix!(H, Hks[spin], args...)
        end
        F = _scf_eigen!(cache, H, 1)
        copyto!(view(electron.Enk, :, spin), F.values)
        if size(electron.Cnk[spin]) != size(F.vectors)
            electron.Cnk[spin] = Matrix{ComplexF64}(F.vectors)
        else
            copyto!(electron.Cnk[spin], F.vectors)
        end
    end
    return nothing
end
