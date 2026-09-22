function Set_Nstate(material::LCPAO_model)
    SpinPol = material.SpinPol
    fsize = sum(material.Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    return Nfsize
end


function Set_Nstate(material::CWF_model)
    Ngsize = material.Ngsize
    return Ngsize
end


"""
Sort tetrahedron vertex energies and carry their associated observable values.

This is the Julia equivalent of OpenMX's `OrderE` routine in
`Tetrahedron_Blochl.c`.
"""
function OrderE!(energies, values, n::Integer)
    @inbounds for i = 1:(Int(n) - 1)
        for j = i:Int(n)
            if energies[j] < energies[i]
                energies[i], energies[j] = energies[j], energies[i]
                values[i], values[j] = values[j], values[i]
            end
        end
    end
    return nothing
end


"""Analytical tetrahedron density of states from Bloechl et al."""
@inline function ATM_Dos(et, energy)
    e1, e2, e3, e4 = et
    energy <= e1 && return 0.0
    energy >= e4 && return 0.0

    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2
    e41 = e4 - e1
    e42 = e4 - e2
    e43 = e4 - e3

    if energy < e2
        delta = energy - e1
        return 3.0 * delta * delta / (e21 * e31 * e41)
    elseif energy < e3
        delta = energy - e2
        return (3.0 * e21 + 6.0 * delta -
                3.0 * (e31 + e42) * delta * delta / (e32 * e42)) /
               (e31 * e41)
    else
        delta = e4 - energy
        return 3.0 * delta * delta / (e41 * e42 * e43)
    end
end


"""
Improved analytical tetrahedron spectrum.

The equations follow Appendix B of P. E. Bloechl, O. Jepsen and
O. K. Andersen, Phys. Rev. B 49, 16223 (1994), as implemented by OpenMX in
`Tetrahedron_Blochl.c`.
"""
@inline function ATM_Spectrum(et, at, energy)
    dos = ATM_Dos(et, energy)
    iszero(dos) && return 0.0

    e1, e2, e3, e4 = et
    a1, a2, a3, a4 = at
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2
    e41 = e4 - e1
    e42 = e4 - e2
    e43 = e4 - e3
    a21 = a2 - a1
    a31 = a3 - a1
    a41 = a4 - a1
    a42 = a4 - a2
    a43 = a4 - a3
    third = 1.0 / 3.0

    if energy < e2
        average = a1 + (energy - e1) * third *
            (a21 / e21 + a31 / e31 + a41 / e41)
    elseif energy < e3
        average = a1 +
            (a21 + e21 * a31 / e31 + e21 * a41 / e41) * third *
            (e3 - energy) / e32 +
            (a4 - a1 -
             (e43 * a41 / e41 + e43 * a42 / e42 + a43) * third) *
            (energy - e2) / e32
    else
        average = a4 + (energy - e4) * third *
            (a41 / e41 + a42 / e42 + a43 / e43)
    end
    return dos * average
end


"""
Return the four linear interpolation weights of `ATM_Spectrum`.

For sorted tetrahedron energies `et`, the spectrum for any vertex observable
`at` is `sum(ATM_Spectrum_weights(et, energy) .* at)`.  Computing these
weights once lets all tensor components and orbital projections reuse the
same tetrahedron algebra.
"""
@inline function ATM_Spectrum_weights(et, energy)
    dos = ATM_Dos(et, energy)
    iszero(dos) && return (0.0, 0.0, 0.0, 0.0)

    e1, e2, e3, e4 = et
    e21 = e2 - e1
    e31 = e3 - e1
    e32 = e3 - e2
    e41 = e4 - e1
    e42 = e4 - e2
    e43 = e4 - e3
    third = 1.0 / 3.0

    if energy < e2
        scale = (energy - e1) * third
        w2 = scale / e21
        w3 = scale / e31
        w4 = scale / e41
        w1 = 1.0 - w2 - w3 - w4
    elseif energy < e3
        lower = third * (e3 - energy) / e32
        upper = (energy - e2) / e32
        r31 = e21 / e31
        r41 = e21 / e41
        r43_41 = e43 / e41
        r43_42 = e43 / e42

        w1 = 1.0 + lower * (-1.0 - r31 - r41) +
             upper * (-1.0 + third * r43_41)
        w2 = lower + upper * third * r43_42
        w3 = lower * r31 + upper * third
        w4 = lower * r41 +
             upper * (1.0 - third * (r43_41 + r43_42 + 1.0))
    else
        scale = (energy - e4) * third
        w1 = -scale / e41
        w2 = -scale / e42
        w3 = -scale / e43
        w4 = 1.0 + scale * (1.0 / e41 + 1.0 / e42 + 1.0 / e43)
    end
    return (dos * w1, dos * w2, dos * w3, dos * w4)
end


function Set_dfdE!(fermi_T_dE, mu, TDFE, TDF_N, kBT; MaxExp = 36.0)
    
    for ie = 1:TDF_N
        ene = TDFE[ie]
        betaE = (ene-mu)/kBT
        if abs(betaE) > MaxExp
            fermi_T_dE[ie] = 0.0
        else
            fermi_T_dE[ie] = 1/kBT*exp(betaE)/((exp(betaE) + 1.0)^2)
        end
    end
end


function Seebeck_inv2!(A::Matrix{Float64}, B::Matrix{Float64})
    
    det = A[1,1]*A[2,2] - A[1,2]*A[2,1]

    B[1,1] = A[2,2]
    B[1,2] = -A[1,2]
    B[2,1] = -A[2,1]
    B[2,2] = A[1,1]

    return det
end


function Seebeck_inv3!(A::Matrix{Float64}, B::Matrix{Float64})
    
    B[1,1] = A[2,2]*A[3,3] - A[3,2]*A[2,3]
    B[1,2] = A[2,3]*A[3,1] - A[3,3]*A[2,1]
    B[1,3] = A[2,1]*A[3,2] - A[3,1]*A[2,2]
    B[2,1] = A[3,2]*A[1,3] - A[1,2]*A[3,3]
    B[2,2] = A[3,3]*A[1,1] - A[1,3]*A[3,1]
    B[2,3] = A[3,1]*A[1,2] - A[1,1]*A[3,2]
    B[3,1] = A[1,2]*A[2,3] - A[2,2]*A[1,3]
    B[3,2] = A[1,3]*A[2,1] - A[2,3]*A[1,1]
    B[3,3] = A[1,1]*A[2,2] - A[2,1]*A[1,2]

    det = A[1,1]*B[1,1] + A[1,2]*B[1,2] + A[1,3]*B[1,3]

    return det
end
