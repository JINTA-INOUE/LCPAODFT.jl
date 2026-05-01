include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test

@testset "Gaunt" begin
    Lmax = 2
    S3J_MAX_FACT = 40
    f = zeros(Float64, S3J_MAX_FACT)
    mult = 1.0
    f[1] = 1.0
    for k = 2:S3J_MAX_FACT
        f[k] = f[k-1]*mult
        mult = mult + 1.0
    end

    check_triangle(l,l1,l2) = (abs(l1 - l2) <= l <= (l1 + l2))
    check_parity(l,l1,l2) = ((l + l1 + l2) % 2 == 0)

    for l = 0:Lmax, m = -l:l, l1 = 0:Lmax, m1 = -l1:l1, l2 = 0:Lmax, m2 = -l2:l2
        if !check_triangle(l,l1,l2) || !check_parity(l,l1,l2)
            g = Gaunt(f,l,m,l1,m1,l2,m2)
            @test isapprox(g, 0.0; atol=1e-12)
        end
    end
end
