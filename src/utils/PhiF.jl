function PhiF(N, r, rv, yv)
    
    mp_min = 1
    mp_max = N

    if r > rv[end]
        return 0.0
    elseif r < rv[begin]
        
        m = 5
        rm = rv[m]

        h1 = rv[m-1] - rv[m-2]
        h2 = rv[m]   - rv[m-1]
        h3 = rv[m+1] - rv[m]

        f1 = yv[m-2]
        f2 = yv[m-1]
        f3 = yv[m]
        f4 = yv[m+1]

        g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
        g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)

        x1 = rm - rv[m-1]
        x2 = rm - rv[m]
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
            if rv[m] < r
                mp_min = m
            else
                mp_max = m
            end
        end
        m = mp_max

        if m < 3
            m = 3
        end
        
        if m >= N
            m = N-1
        end

        h1 = rv[m-1] - rv[m-2]
        h2 = rv[m]   - rv[m-1]
        h3 = rv[m+1] - rv[m]

        f1 = yv[m-2]
        f2 = yv[m-1]
        f3 = yv[m]
        f4 = yv[m+1]

        if m == 2
            h1 = -(h2+h3)
            f1 = f4
        end

        if m == N
            h3 = -(h1+h2)
            f4 = f1
        end

        g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2)
        g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3)

        x1 = r - rv[m-1]
        x2 = r - rv[m]
        y1 = x1/h2
        y2 = x2/h2

        return y2*y2*(3*f2 + h2*g1 + (2*f2 + h2*g1)*y2) + y1*y1*(3*f3 - h2*g2 - (2*f3 - h2*g2)*y1)
    end 
end
