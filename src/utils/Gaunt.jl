function _Set_f_for_Gaunt!(f::Vector{Float64})
    mult = 1.0
    f[1] = 1.0
    for k = 2:S3J_MAX_FACT
        f[k] = f[k-1]*mult
        mult = mult + 1.0
    end
end


@inline function Gaunt(f::Vector{Float64}, l::Int, m::Int, l1::Int, m1::Int, l2::Int, m2::Int)

    cleb1 = Clebsch_Gordan(f, l1, 0, l2, 0, l, 0)
    cleb2 = Clebsch_Gordan(f, l1, m1, l2, m2, l, m)

    temp0 = 2*l1 + 1
    temp1 = 2*l2 + 1
    temp2 = 4*pi*(2*l + 1)
    temp3 = sqrt(temp0*temp1/temp2)

    return temp3*cleb1*cleb2
end


@inline function Clebsch_Gordan(f::Vector{Float64}, j1::Int, m1::Int, j2::Int, m2::Int, j::Int, m::Int)
    
    esp = j1 - j2 + m
    if !(abs(esp - (j1 - j2 + m)) < 1e-10)
        return 0.0
    end

    cgris = ifelse(esp%2 == 0, 1.0, -1.0)
    cgris *= sqrt(2*j+1)*s3j(f, j1, j2, j, m1, m2, -m)

    return cgris
end


@inline function s3j(f::Vector{Float64}, j1::Int, j2::Int, j3::Int, m1::Int, m2::Int, m3::Int)
   
    jpm1 = j1 + m1
    jpm2 = j2 + m2
    jpm3 = j3 + m3
    jmm1 = j1 - m1
    jmm2 = j2 - m2
    jmm3 = j3 - m3
    
    if (jpm1-jmm1+jpm2-jmm2+jpm3-jmm3) ≠ 0
        return 0.0
    end


    j1pj2mj3 = jpm1 + jpm2 - jmm3
    j3mj2pm1 = jmm3 - jpm2
    j3mj1mm2 = jpm3 - jmm1

    
    kmin = max(-j3mj2pm1, -j3mj1mm2, 0)
    kmax = min(j1pj2mj3, jmm1, jpm2)

    if kmin > kmax
        return 0.0
    end

    mult = 1.0
    if iszero(kmin%2)
        mult = 1.0
    else
        mult = -1.0
    end
    ris = 0.0

    for k = kmin:kmax
        ris += mult/(f[k+1]*f[j1pj2mj3-k+1]*f[jmm1-k+1]*f[jpm2-k+1]*f[j3mj2pm1+k+1]*f[j3mj1mm2+k+1])
        mult = -mult
    end

    if ((jpm1-jmm2)%2) ≠ 0
        ris = -ris
    end

    ris *= sqrt(f[j1pj2mj3+1]*f[jpm1-jmm2+jpm3+1]*f[-jmm1+jpm2+jpm3+1]*f[jpm1+1]*f[jpm2+1]*f[jpm3+1]*f[jmm1+1]*f[jmm2+1]*f[jmm3+1]/f[jpm1+jpm2+jpm3+1+1])

    return ris
end