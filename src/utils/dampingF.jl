function dampingF(rcut::Float64, r::Float64)

    buf = 1.0

    r01 = -buf
    r02 = r01*r01
    r03 = r02*r01
    r04 = r03*r01
    r05 = r04*r01

    A11 = r03
    A12 = r04
    A13 = r05
    A21 = 3*r02
    A22 = 4*r03
    A23 = 5*r04
    A31 = 6*r01
    A32 = 12*r02
    A33 = 20*r03

    detA = A11*A22*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31 - A12*A21*A33 - A11*A23*A32
  
    c3 =  (A22*A33-A23*A32)/detA
    c4 = -(A21*A33-A23*A31)/detA
    c5 =  (A21*A32-A22*A31)/detA

    r0 = r - rcut
    r3 = r0^3
    r4 = r0^4
    r5 = r0^5

    if rcut < r
        return 0.0
    elseif r < rcut-buf
        return 1.0
    else
        return c3*r3 + c4*r4 + c5*r5
    end
end


function deri_dampingF(rcut::Float64, r::Float64)

    buf = 1.0

    r01 = -buf
    r02 = r01*r01
    r03 = r02*r01
    r04 = r03*r01
    r05 = r04*r01

    A11 = r03
    A12 = r04
    A13 = r05
    A21 = 3*r02
    A22 = 4*r03
    A23 = 5*r04
    A31 = 6*r01
    A32 = 12*r02
    A33 = 20*r03

    detA = A11*A22*A33 + A12*A23*A31 + A13*A21*A32 - A13*A22*A31 - A12*A21*A33 - A11*A23*A32
  
    c3 =  (A22*A33-A23*A32)/detA
    c4 = -(A21*A33-A23*A31)/detA
    c5 =  (A21*A32-A22*A31)/detA

    r0 = r - rcut
    r2 = r0^2
    r3 = r0^3
    r4 = r0^4

    if rcut < r
        return 0.0
    elseif r < rcut-buf
        return 0.0
    else
        return 3*c3*r2 + 4*c4*r3 + 5*c5*r4
    end
end