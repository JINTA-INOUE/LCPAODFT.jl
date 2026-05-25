@inline function Dr_VH_AtomF(Spe_Core_Charge, N, x, r, xv, rv, yv)
    
    xmin = xv[1]
    xmax = xv[N]

    if xmax <= x
        return -Spe_Core_Charge/r^2
    elseif r < rv[1]
        
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

        df = (2*y2*E + y22*B + 2*y1*F - y12*D)/h2

        a = 0.5*df/rm

        return 2*a*r
    else
        x = max(x, xmin)
        tmp1 = (N - 1)/(xmax - xmin)
        tmp2 = (x - xmin)*tmp1
        i = floor(Int, tmp2)
        dt = tmp2 - i
        i = i + 1

        d0 = yv[i+1] - yv[i]
        d1 = yv[i+2] - yv[i+1]
        d2 = yv[i+3] - yv[i+2]

        a = d0 - 2*d1 + d2
        b = -2*d0 + 3*d1 - d2
        c = d0 + d1

        return 0.5*((3*a*dt + 2*b)*dt + c)*tmp1/r
    end 
end
