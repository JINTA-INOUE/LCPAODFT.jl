function Gauss_Legendre(N)
    
    EPS = 1e-13
    x = zeros(Float64, N)
    weight = zeros(Float64, N)

    m = div(N+1, 2)
    xm = 0.0
    xl = 1.0

    for i = 1:m

        z = cos(pi*(i-0.25)/(N+0.5))
        z1 = 0.0
        pp = 1.0

        while abs(z-z1) > EPS

            p1 = 1.0
            p2 = 0.0

            for j = 1:N
                p3 = p2
                p2 = p1
                p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
            end

            pp = N*(z*p1-p2)/(z^2 - 1)
            z1 = z
            z = z1 - p1/pp
        end

        p1 = 1.0
        p2 = 0.0

        for j = 1:N
            p3 = p2
            p2 = p1
            p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
        end

        pp = N*(z*p1-p2)/(z^2 - 1)
        z1 = z
        z = z1 - p1/pp

        x[i] = xm - xl*z
        x[N-i+1] = xm + xl*z
        weight[i] = 2*xl/((1-z^2)*pp^2)
        weight[N-i+1] = weight[i]    
    end


    return x, weight
end


function Gauss_Legendre_x(N)
    
    EPS = 1e-13
    x = zeros(Float64, N)
    m = div(N+1, 2)
    xm = 0.0
    xl = 1.0

    for i = 1:m

        z = cos(pi*(i-0.25)/(N+0.5))
        z1 = 0.0
        pp = 1.0

        while abs(z-z1) > EPS

            p1 = 1.0
            p2 = 0.0

            for j = 1:N
                p3 = p2
                p2 = p1
                p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
            end

            pp = N*(z*p1-p2)/(z^2 - 1)
            z1 = z
            z = z1 - p1/pp
        end

        p1 = 1.0
        p2 = 0.0

        for j = 1:N
            p3 = p2
            p2 = p1
            p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
        end

        pp = N*(z*p1-p2)/(z^2 - 1)
        z1 = z
        z = z1 - p1/pp

        x[i] = xm - xl*z
        x[N-i+1] = xm + xl*z 
    end


    return x
end
