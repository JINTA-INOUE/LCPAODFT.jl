function RadialF(l, Num_Mesh_PAO, r, PAO_RV, PAO_RWF)

    mp_min = 1
    mp_max = Num_Mesh_PAO

    if r < PAO_RV[1]
        if l == 0
            return PAO_RWF[1]
        else
            return 0.0
        end

    elseif PAO_RV[end] < r
        return 0.0
    else
        while (mp_max - mp_min) ≠ 1
            
            m = div(mp_min + mp_max, 2)
            if PAO_RV[m] < r
                mp_min = m
            else
                mp_max = m
            end
        end
        m = mp_max

        if m < 3
            m = 3
        elseif m >= Num_Mesh_PAO
            m = Num_Mesh_PAO-1
        end
        
        if m == 2
            h2 = PAO_RV[m] - PAO_RV[m-1]
            h3 = PAO_RV[m+1] - PAO_RV[m]
            f2 = PAO_RWF[m-1]
            f3 = PAO_RWF[m]
            f4 = PAO_RWF[m+1]
            h1 = -(h2 + h3)
            f1 = f4
        elseif m == Num_Mesh_PAO
            h1 = PAO_RV[m-1] - PAO_RV[m-2]
            h2 = PAO_RV[m] - PAO_RV[m-1]
            f1 = PAO_RWF[m-2]
            f2 = PAO_RWF[m-1]
            f3 = PAO_RWF[m]
            h3 = -(h1 + h2)
            f4 = f1
        else
            h1 = PAO_RV[m-1] - PAO_RV[m-2]
            h2 = PAO_RV[m] - PAO_RV[m-1]
            h3 = PAO_RV[m+1] - PAO_RV[m]
            f1 = PAO_RWF[m-2]
            f2 = PAO_RWF[m-1]
            f3 = PAO_RWF[m]
            f4 = PAO_RWF[m+1]
        end

        g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
        g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)
    
        x1 = r - PAO_RV[m-1]
        x2 = r - PAO_RV[m]
        y1 = x1/h2
        y2 = x2/h2
    
        return y2*y2*(3*f2 + h2*g1 + (2*f2 + h2*g1)*y2) + y1*y1*(3*f3 - h2*g2 - (2*f3 - h2*g2)*y1)
    end
end