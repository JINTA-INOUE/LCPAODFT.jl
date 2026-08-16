function VNAF(Num_Mesh_VPS, r, Spe_Atom_Cut1, VPS_RV, Vna)

    mp_min = 1
    mp_max = Num_Mesh_VPS

    if Spe_Atom_Cut1 < r
        return 0.0
    elseif r < VPS_RV[1]

        m = 5
        rm = VPS_RV[m]

        h1 = VPS_RV[m-1] - VPS_RV[m-2]
        h2 = VPS_RV[m]   - VPS_RV[m-1]
        h3 = VPS_RV[m+1] - VPS_RV[m]

        f1 = Vna[m-2]
        f2 = Vna[m-1]
        f3 = Vna[m]
        f4 = Vna[m+1]

        g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
        g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)

        x1 = rm - VPS_RV[m-1]
        x2 = rm - VPS_RV[m]
        y1 = x1/h2
        y2 = x2/h2
        y12 = y1*y1
        y22 = y2*y2

        A = 3*f2 + h2*g1
        B = 2*f2 + h2*g1
        C = 3*f3 - h2*g2
        D = 2*f3 - h2*g2
        E = A + B*y2
        F = C - D*y1

        f = y22*E + y12*F
        df = (2*y2*E + y22*B + 2*y1*F - y12*D)/h2

        a = 0.5*df/rm
        b = f - a*rm^2

        return a*r*r + b
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

        if m < 3
            m = 3
        elseif m >= Num_Mesh_VPS
            m = Num_Mesh_VPS-1
        end

        h1 = VPS_RV[m-1] - VPS_RV[m-2]
        h2 = VPS_RV[m]   - VPS_RV[m-1]
        h3 = VPS_RV[m+1] - VPS_RV[m]

        f1 = Vna[m-2]
        f2 = Vna[m-1]
        f3 = Vna[m]
        f4 = Vna[m+1]

        if m == 2
            h1 = -(h2+h3)
            f1 = f4
        end

        if m == Num_Mesh_VPS
            h3 = -(h1+h2)
            f4 = f1
        end

        g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
        g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)

        x1 = r - VPS_RV[m-1]
        x2 = r - VPS_RV[m]
        y1 = x1/h2
        y2 = x2/h2
        A = h2*g1
        B = h2*g2
        result = y2*y2*(3*f2 + A + (2*f2 + A)*y2) + y1*y1*(3*f3 - B - (2*f3 - B)*y1)

        return result
    end
end
