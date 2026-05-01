function LSDA_CA(denup::Float64, dendown::Float64, switch::Integer)

    tden = denup + dendown
    min_den = 1e-15
    if tden <= min_den
        return zeros(Float64, 2)
    else
        zeta = (denup-dendown)/tden
        zeta = ifelse(zeta>1.0, 1-min_den, zeta)
        zeta = ifelse(zeta<-1.0, -1+min_den, zeta)
        
        coe = 0.6203504908994

        rs = ifelse(tden<min_den, 62035.04908994, coe*tden^(-1/3))

        tmp0 = 0.458165293632163/rs
        ExP = -tmp0
        dExP = tmp0/rs

        ExF = 1.25992104989487*ExP
        dExF = 1.25992104989487*dExP

        if rs >= 1.0
            tmp0 = sqrt(rs)
            dum = (1.0 + 1.0529*tmp0 + 0.3334*rs)
            tmp1 = 0.1423/dum
            EcP = -tmp1
            dEcP = tmp1/dum*(0.52645/tmp0 + 0.3334)

            dum = (1.0 + 1.3981*tmp0 + 0.2611*rs)
            tmp1 = 0.0843/dum
            EcF = -tmp1
            dEcF = tmp1/dum*(0.69905/tmp0 + 0.2611)
        else
            tmp0 = log(rs)
            EcP = -0.0480 + 0.0311*tmp0 + rs*(0.0020*tmp0 - 0.0116)
            dEcP = 0.0311/rs + 0.0020*tmp0 - 0.0096

            EcF = -0.0269 + 0.01555*tmp0 + rs*(0.0007*tmp0 - 0.0048)
            dEcF = 0.01555/rs + 0.0007*tmp0 - 0.0041
        end

        if switch == 0

            z0 = (1 + zeta)^(4/3)
            z1 = (1 - zeta)^(4/3)
            fzeta = 1.92366105093154*(z0 + z1 - 2.0)

            Exc = ExP + EcP + (ExF + EcF - ExP - EcP)*fzeta

            return Exc, Exc

        elseif switch == 1

            z0 = (1 + zeta)^(1/3)
            z1 = (1 - zeta)^(1/3)

            z02 = z0*z0
            z04 = z02*z02
            z12 = z1*z1
            z14 = z12*z12
            fzeta  = 1.92366105093154*(z04 + z14 - 2.0)
            dfzeta = 2.56488140124205*(z0 - z1)

            tmp0 = ExF + EcF - ExP - EcP
            Exc = ExP + EcP + tmp0*fzeta
            dExc = dExP + dEcP + (dExF + dEcF - dExP - dEcP)*fzeta
            Vxc = Exc - 1/3*rs*dExc

            
            return Vxc + tmp0*(1-zeta)*dfzeta, Vxc + tmp0*(-1-zeta)*dfzeta

        elseif switch == 2

            z0 = (1 + zeta)^(1/3)
            z1 = (1 - zeta)^(1/3)

            z02 = z0*z0
            z04 = z02*z02
            z12 = z1*z1
            z14 = z12*z12
            fzeta  = 1.92366105093154*(z04 + z14 - 2.0)
            dfzeta = 2.56488140124205*(z0 - z1)

            tmp0 = ExF + EcF - ExP - EcP
            Exc = ExP + EcP + tmp0*fzeta
            dExc = dExP + dEcP + (dExF + dEcF - dExP - dEcP)*fzeta
            tmp1 = 1/3*rs*dExc
            tmp2 = tmp0*dfzeta
            Vxc = Exc - tmp1

            return tmp1 - tmp2*(1-zeta), tmp1 - tmp2*(-1-zeta)
        else
            error("please check switch")
        end
    end
end
