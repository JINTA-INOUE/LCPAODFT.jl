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