"""
    calulation VXC using Ceperly Alder
"""
function LDA_CA(den::Float64, switch::Integer)
    
    if den <= 1e-15
        return 0.0
    else
        # rs = ((3/4/pi))^(1/3)*den^(-1/3)
        coe = 0.6203504908994
        rs = coe*den^(-1/3)

        tmp0 = 0.458165293632163/rs
        Ex = -tmp0
        dEx = tmp0/rs

        if rs >= 1.0
            tmp0 = sqrt(rs)
            dum = (1.0 + 1.0529*tmp0 + 0.3334*rs)
            tmp1 = 0.1423/dum
            Ec = -tmp1
            dEc = tmp1/dum*(0.52645/tmp0 + 0.3334)
        else
            tmp0 = log(rs)
            Ec = -0.0480 + 0.0311*tmp0 + rs*(0.0020*tmp0 - 0.0116)
            dEc = 0.0311/rs + 0.0020*tmp0 - 0.0096
        end

        if switch == 0
           return Ex + Ec
        elseif switch == 1
            return Ex + Ec - 1/3*rs*(dEx + dEc)
        elseif switch == 2
            return 1/3*rs*(dEx + dEc)
        elseif switch == 3
            return -1/3/coe/coe/coe*rs*rs*rs*rs*(dEx + dEc)
        else
            error("please check switch")
        end
    end
end
