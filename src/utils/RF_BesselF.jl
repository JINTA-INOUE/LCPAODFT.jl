function RF_BesselF(k, NormK, Spe_RF_BesselF)

    m = 0
    mp_min = 1
    mp_max = NkGrid

    if k < NormK[begin]
        m = 1
    elseif NormK[end] < k
        return 0.0
    else
        while (mp_max - mp_min) ≠ 1
            
            m = div(mp_min + mp_max, 2)
            if NormK[m] < k
                mp_min = m
            else
                mp_max = m
            end
        end
        m = mp_max

        if m < 3
            m = 3
        elseif m >= NkGrid+1
            m = NkGrid-1
        end
    end
        
        
    if m == 2
        h2 = NormK[m] - NormK[m-1]
        h3 = NormK[m+1] - NormK[m]
        f2 = Spe_RF_BesselF[m-1]
        f3 = Spe_RF_BesselF[m]
        f4 = Spe_RF_BesselF[m+1]
        h1 = -(h2 + h3)
        f1 = f4
    elseif m == NkGrid
        h1 = NormK[m-1] - NormK[m-2]
        h2 = NormK[m] - NormK[m-1]
        f1 = Spe_RF_BesselF[m-2]
        f2 = Spe_RF_BesselF[m-1]
        f3 = Spe_RF_BesselF[m]
        h3 = -(h1 + h2)
        f4 = f1
    else
        h1 = NormK[m-1] - NormK[m-2]
        h2 = NormK[m] - NormK[m-1]
        h3 = NormK[m+1] - NormK[m]
        f1 = Spe_RF_BesselF[m-2]
        f2 = Spe_RF_BesselF[m-1]
        f3 = Spe_RF_BesselF[m]
        f4 = Spe_RF_BesselF[m+1]
    end

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)
    
    x1 = k - NormK[m-1]
    x2 = k - NormK[m]
    y1 = x1/h2
    y2 = x2/h2
    
    return y2*y2*(3*f2 + h2*g1 + (2*f2 + h2*g1)*y2) + y1*y1*(3*f3 - h2*g2 - (2*f3 - h2*g2)*y1)
end