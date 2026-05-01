function Nonlocal_RadialF(GL, Num_Mesh_VPS, r, VPS_RV, Spe_VNL)

    mp_min = 1
    mp_max = Num_Mesh_VPS

    if r < VPS_RV[1]
        if GL == 0
            m = 2
        else
            return 0.0
        end
    elseif VPS_RV[end] < r
        return 0.0
    else
        while (mp_max - mp_min) ≠ 1
            
            m = div(mp_min + mp_max, 2)
            if VPS_RV[m] < r
                mp_min = m
            else
                mp_max = m
            end
        end
        m = mp_max
    end

    if m == 2
        h2 = VPS_RV[m] - VPS_RV[m-1]
        h3 = VPS_RV[m+1] - VPS_RV[m]
        f2 = Spe_VNL[m-1]
        f3 = Spe_VNL[m]
        f4 = Spe_VNL[m+1]
        h1 = -(h2 + h3)
        f1 = f4
    elseif m == Num_Mesh_VPS
        h1 = VPS_RV[m-1] - VPS_RV[m-2]
        h2 = VPS_RV[m] - VPS_RV[m-1]
        f1 = Spe_VNL[m-2]
        f2 = Spe_VNL[m-1]
        f3 = Spe_VNL[m]
        h3 = -(h1 + h2)
        f4 = f1
    else
        h1 = VPS_RV[m-1] - VPS_RV[m-2]
        h2 = VPS_RV[m] - VPS_RV[m-1]
        h3 = VPS_RV[m+1] - VPS_RV[m]
        f1 = Spe_VNL[m-2]
        f2 = Spe_VNL[m-1]
        f3 = Spe_VNL[m]
        f4 = Spe_VNL[m+1]
    end

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)
    
    x1 = r - VPS_RV[m-1]
    x2 = r - VPS_RV[m]
    y1 = x1/h2
    y2 = x2/h2
    
    return y2*y2*(3*f2 + h2*g1 + (2*f2 + h2*g1)*y2) + y1*y1*(3*f3 - h2*g2 - (2*f3 - h2*g2)*y1)
end