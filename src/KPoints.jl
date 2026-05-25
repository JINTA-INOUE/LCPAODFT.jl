"""
The type for describing Bloch wave vector of electronic states for MPI
"""
struct KPoints
    AllNkpt::Int32
    Nkpt::Int32
    MPI_Nkpt::Int32
    kmesh::Tuple{Int32,Int32,Int32}
    MPI_kpts::Vector{Vector{Float64}}
    MPI_kweight::Vector{Int32}
    All_kweight::Vector{Int32}
    MPI_krange::Vector{UnitRange{Int32}}
    MPkpts::Vector{Int32}
    KP_flag::String
    crystal_sym::Bool
    time_rev::Bool
end


function KPoints(system::String, SpinPol::String, kmesh::Tuple{Signed,Signed,Signed}, time_rev::Bool; KP_flag="Gcenter")
    
    if system == "Crystal"
        if SpinPol ≠ "nc"
            # only using time reversal symmetry
            kpoints = KPoints(kmesh, time_rev; KP_flag)
        else
            # NonCollinear case is no time reversal symmetry
            kpoints = KPoints(kmesh, false; KP_flag)
        end
    elseif system ∈ ("Atom", "Cluster")
		kpoints = nothing
	else
		error("now system = $system. please check system")
	end

    
    return kpoints
end



"""
    MPI_KPoints(...)

Generate mesh in BZ.

Mandatory arguments:

- `LatVecs`: Lattice vectrors
- `kmesh`: mesh number in BZ

The following is the most commonly used optional arguments:
- `time_rev` : is time reversing
"""
function KPoints(kmesh::Tuple{Signed,Signed,Signed}, time_rev::Bool; KP_flag="Gcenter")

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    AllNkpt = prod(kmesh)
    Nkpt, kpts, kweight = Gen_KPoints(kmesh, time_rev; KP_flag)

    MPI_krange = split_evenly(1:Nkpt, nprocs)

    MPI_Nkpt = length(MPI_krange[myrank+1])
    MPI_kpts = kpts[MPI_krange[myrank+1]]
    MPI_kweight = kweight[MPI_krange[myrank+1]]

    MPI_Nkptsize = zeros(Int32, nprocs)
    MPkpts = zeros(Int32, nprocs)
    MPI_Nkptsize[myrank+1] = MPI_Nkpt
    
    MPI.Allreduce!(MPI_Nkptsize, MPI.SUM, comm)

    MPkpts = zeros(Int32, nprocs)
    Sum = 0
    for id = 1:nprocs
        MPkpts[id] = Sum
        Sum += MPI_Nkptsize[id]
    end
   

    return KPoints(AllNkpt, Nkpt, MPI_Nkpt, 
                   kmesh, MPI_kpts, 
                   MPI_kweight, kweight, 
                   MPI_krange, MPkpts,
                   KP_flag, false, time_rev)
end


function Gen_KPoints(kmesh::Tuple{Signed,Signed,Signed}, time_rev::Bool; Shift_K_Point=1.0e-12, KP_flag="Gcenter")

    if lowercase(KP_flag) == "gcenter"
        if time_rev
            return _Gen_KPoints_Gcenter_TRS(kmesh, Shift_K_Point)
        else
            return _Gen_KPoints_Gcenter_noTRS(kmesh, Shift_K_Point)
        end
    elseif lowercase(KP_flag) == "mp"
        if time_rev
            return _Gen_KPoints_MP_TRS(kmesh, Shift_K_Point)
        else
            return _Gen_KPoints_MP_noTRS(kmesh, Shift_K_Point)
        end
    else
        error("please check KP_flag")
    end
end


function _Gen_KPoints_Gcenter_TRS(kmesh::Tuple{Signed,Signed,Signed}, Shift_K_Point)

    kmesh1, kmesh2, kmesh3 = kmesh
    kop = zeros(Int64, kmesh1, kmesh2, kmesh3)
    fill!(kop, -999)

    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        if kop[i+1,j+1,k+1] == -999
            ii = ifelse(i==0 || 2*i==kmesh1, i, kmesh1-i)
            ij = ifelse(j==0 || 2*j==kmesh2, j, kmesh2-j)
            ik = ifelse(k==0 || 2*k==kmesh3, k, kmesh3-k)
            
            if (i==0 || 2*i==kmesh1) && (j==0 || 2*j==kmesh2) && (k==0 || 2*k==kmesh3)
                kop[i+1,j+1,k+1] = 1
            else
                kop[i+1,j+1,k+1] = 2
                kop[ii+1,ij+1,ik+1] = 0
            end
        end
    end

    Nkpt = 0
    for i = 1:kmesh1, j = 1:kmesh2, k = 1:kmesh3
        if kop[i,j,k] > 0
            Nkpt += 1
        end
    end

    kweight = zeros(Int64, Nkpt)
    kpts = Vector{Vector{Float64}}(undef, Nkpt)
    for ik = 1:Nkpt
        kpts[ik] = zeros(Float64, 3)
    end


    ik = 0
    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        if kop[i+1,j+1,k+1] > 0
            ik += 1
            k1 = ifelse(kmesh1==1, 0.0, i/kmesh1+Shift_K_Point)
            k2 = ifelse(kmesh2==1, 0.0, j/kmesh2-Shift_K_Point)
            k3 = ifelse(kmesh3==1, 0.0, k/kmesh3+2*Shift_K_Point)

            kpts[ik][1] = k1
            kpts[ik][2] = k2
            kpts[ik][3] = k3

            kweight[ik] = kop[i+1,j+1,k+1]
        end
    end


    return Nkpt, kpts, kweight
end


function _Gen_KPoints_Gcenter_noTRS(kmesh::Tuple{Signed,Signed,Signed}, Shift_K_Point)

    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    kweight = ones(Int64, Nkpt)
    kpts = Vector{Vector{Float64}}(undef, Nkpt)
    for ik = 1:Nkpt
        kpts[ik] = zeros(Float64, 3)
    end


    ik = 0
    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        ik += 1
        k1 = ifelse(kmesh1==1, 0.0, i/kmesh1+Shift_K_Point)
        k2 = ifelse(kmesh2==1, 0.0, j/kmesh2-Shift_K_Point)
        k3 = ifelse(kmesh3==1, 0.0, k/kmesh3+2*Shift_K_Point)

        kpts[ik][1] = k1
        kpts[ik][2] = k2
        kpts[ik][3] = k3
    end


    return Nkpt, kpts, kweight
end


function _Gen_KPoints_MP_TRS(kmesh::Tuple{Signed,Signed,Signed}, Shift_K_Point)

    kmesh1, kmesh2, kmesh3 = kmesh
    kop = zeros(Int64, kmesh1, kmesh2, kmesh3)
    fill!(kop, -999)

    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        if kop[i+1,j+1,k+1] == -999
            ii = kmesh1 - i - 1
            ij = kmesh2 - j - 1
            ik = kmesh3 - k - 1
            
            if isequal(i,ii) && isequal(j,ij) && isequal(k,ik)
                kop[i+1,j+1,k+1] = 1
            else
                kop[i+1,j+1,k+1] = 2
                kop[ii+1,ij+1,ik+1] = 0
            end
        end
    end

    Nkpt = 0
    for i = 1:kmesh1, j = 1:kmesh2, k = 1:kmesh3
        if kop[i,j,k] > 0
            Nkpt += 1
        end
    end

    kweight = zeros(Int64, Nkpt)
    kpts = Vector{Vector{Float64}}(undef, Nkpt)
    for ik = 1:Nkpt
        kpts[ik] = zeros(Float64, 3)
    end


    ik = 0
    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        if kop[i+1,j+1,k+1] > 0
            ik += 1
            k1 = ifelse(kmesh1==1, 0.0, -0.5 + (2*i+1)/kmesh1/2 + Shift_K_Point)
            k2 = ifelse(kmesh2==1, 0.0, -0.5 + (2*j+1)/kmesh2/2 - Shift_K_Point)
            k3 = ifelse(kmesh3==1, 0.0, -0.5 + (2*k+1)/kmesh3/2 + 2*Shift_K_Point)

            kpts[ik][1] = k1
            kpts[ik][2] = k2
            kpts[ik][3] = k3

            kweight[ik] = kop[i+1,j+1,k+1]
        end
    end

    return Nkpt, kpts, kweight
end


function _Gen_KPoints_MP_noTRS(kmesh::Tuple{Signed,Signed,Signed}, Shift_K_Point)

    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)

    kweight = ones(Int64, Nkpt)
    kpts = Vector{Vector{Float64}}(undef, Nkpt)
    for ik = 1:Nkpt
        kpts[ik] = zeros(Float64, 3)
    end


    ik = 0
    for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
        ik += 1
        k1 = ifelse(kmesh1==1, 0.0, -0.5 + (2*i+1)/kmesh1/2 + Shift_K_Point)
        k2 = ifelse(kmesh2==1, 0.0, -0.5 + (2*j+1)/kmesh2/2 - Shift_K_Point)
        k3 = ifelse(kmesh3==1, 0.0, -0.5 + (2*k+1)/kmesh3/2 + 2*Shift_K_Point)

        kpts[ik][1] = k1
        kpts[ik][2] = k2
        kpts[ik][3] = k3
    end
    

    return Nkpt, kpts, kweight
end
