function dampingF(rcut::Float64, r::Float64)

    c3 = -10.0
    c4 = -15.0
    c5 = -6.0

    r0 = r - rcut
    r3 = r0^3
    r4 = r0^4
    r5 = r0^5

    if rcut < r
        return 0.0
    elseif r < rcut-1.0
        return 1.0
    else
        return c3*r3 + c4*r4 + c5*r5
    end
end


function deri_dampingF(rcut, r)

    buf = 1.0

    r01 = -buf
    r02 = r01*r01
    r03 = r02*r01
    r04 = r03*r01
    r05 = r04*r01

    A = zeros(Float64,3,3)
    A[1,1] = r03
    A[1,2] = r04
    A[1,3] = r05
    A[2,1] = 3*r02
    A[2,2] = 4*r03
    A[2,3] = 5*r04
    A[3,1] = 6*r01
    A[3,2] = 12*r02
    A[3,3] = 20*r03

    Ainv = inv(A)

    c3 = Ainv[1,1]
    c4 = Ainv[2,1]
    c5 = Ainv[3,1]

    r0 = r - rcut
    r2 = r0^2
    r3 = r0^3
    r4 = r0^4
    r5 = r0^5

    if rcut < r
        return 0.0
    elseif r < rcut-buf
        return 0.0
    else
        return 3*c3*r2 + 4*c4*r3 + 5*c5*r4
    end
end