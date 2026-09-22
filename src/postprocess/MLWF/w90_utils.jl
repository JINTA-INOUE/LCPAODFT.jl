function gauge_matrices_dis(chk::w90_Chk)
    Nkpt = chk.Nkpt
    BANDNUM = chk.BANDNUM
    WANNUM = chk.WANNUM

    T = eltype(chk.Uml[1])
    if !chk.have_disentangled
        return [diagm(0 => ones(T, WANNUM)) for _ in 1:Nkpt]
    end

    Iᵏ = Matrix{T}(I, BANDNUM, BANDNUM)

    return map(1:Nkpt) do ik
        p = sortperm(chk.dis_bands[ik]; order = Base.Order.Reverse)
        Iᵏ[:, p] * chk.Udis[ik]
    end
end


function gauge_matrices(chk::w90_Chk)

    Nkpt = chk.Nkpt
    have_disentangled = chk.have_disentangled
    Uml = deepcopy(chk.Uml)
    Udis = gauge_matrices_dis(chk)

    if !have_disentangled
        return Uml
    end


    Umnk = deepcopy(Udis)
    for ik = 1:Nkpt
        _Udis = Udis[ik]
        _Uml = Uml[ik]
        mul!(Umnk[ik], _Udis, _Uml)
    end

    return Umnk
end