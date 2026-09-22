struct w90_Chk
    header::String
    BANDNUM::Int32
    n_exclude_bands::Int32
    exclude_bands::Vector{Int64}
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    Nkpt::Int32
    kmesh::Vector{Int32}
    kpts::Vector{Vector{Float64}}
    n_bvecs::Int32
    WANNUM::Int32
    checkpoint::String
    have_disentangled::Bool
    ΩI::Float64
    dis_bands::Vector{BitVector}
    n_dis::Vector{Int32}
    Udis::Vector{Matrix{ComplexF64}}
    Uml::Vector{Matrix{ComplexF64}}
    Mmnkb::Vector{Vector{Matrix{ComplexF64}}}
    WannierCenter::Matrix{Float64}
    Omega_tot::Vector{Float64}
end


function Print_chk(chk::w90_Chk)

    header = chk.header
    WANNUM = chk.WANNUM
    have_disentangled = chk.have_disentangled
    ΩI = chk.ΩI
    WannierCenter = chk.WannierCenter
    Omega_tot = chk.Omega_tot
end


function parse_bool(i::Integer)
    return i != 0
end

