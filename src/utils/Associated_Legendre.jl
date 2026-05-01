function Associated_Legendre(l, m, x; cut0=1e-24)
    
    if m < 0 || m > l || abs(x) > 1.0
        error("Invalid l, m, x")
    elseif (1.0-cut0)<abs(x)
        x = sign(x)*(1.0 - cut0)
    end


    Pm = 1.0
    if m > 0
        f = 1.0
        tmp0 = sqrt((1-x)*(1+x))
        for _ = 1:m
          Pm = -Pm*f*tmp0
          f += 2.0
        end
    end
        
    if isequal(l, m)
        p0 = Pm
        return p0
    else

        # calculate Pm1
        Pm1 = x*(2*m + 1)*Pm
        if l == (m+1)
            p0 = Pm1
            return p0
        else
            for ll = m+2:l
                tmp0 = (x*(2*ll-1)*Pm1 - (ll+m-1)*Pm)/(ll-m)
                Pm  = Pm1
                Pm1 = tmp0
            end
            p0 = Pm1
            return p0
        end
    end
end


function Associated_Legendre2(l, m, x; cut0=1e-24, cut1=1e-12)
    
    if m < 0 || m > l || abs(x) > 1.0
        error("Invalid l, m, x")
    elseif (1.0-cut0)<abs(x)
        x = sign(x)*(1.0 - cut0)
    end


    Pm = 1.0
    if m > 0
        f = 1.0
        tmp0 = sqrt((1.0-x)*(1.0+x))
        for _ = 1:m
            Pm = -Pm*f*tmp0
            f += 2.0
        end
    end
        

    if isequal(l, m)
        p0 = Pm
        p1 = 0.0
        tmp0 = sqrt(1.0-x^2)
        dP = ifelse(tmp0 > cut1, (l*x*p0 - (l+m)*p1)/tmp0, 0.0)

        return p0, dP
    else

        # calculate Pm1
        Pm1 = x*(2*m + 1)*Pm
        if l == (m+1)
            p0 = Pm1
            p1 = Pm
            tmp0 = sqrt(1.0-x^2)
            dP = ifelse(tmp0 > cut1, (l*x*p0 - (l+m)*p1)/tmp0, 0.0)

            return p0, dP
        else
            for ll = m+2:l
                tmp0 = (x*(2*ll-1)*Pm1 - (ll+m-1)*Pm)/(ll-m)
                Pm  = Pm1
                Pm1 = tmp0
            end
            p0 = Pm1
            p1 = Pm

            tmp0 = sqrt(1.0-x^2)
            dP = ifelse(tmp0 > cut1, (l*x*p0 - (l+m)*p1)/tmp0, 0.0)

            return p0, dP
        end
    end
end