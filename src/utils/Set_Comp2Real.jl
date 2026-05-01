"""
    Set_Comp2Real.jl
    ylm :: Real Spherical Harmonics
    Ylm :: Complex Spherical Harmonics

    l = 0
    s = y₀₀ = Y₀₀

    l = 1
    py = y₁₋₁ = im/√2*( Y₁₋₁ + Y₁₁ )
    pz = y₁₀  = Y₁₀
    px = y₁₁  = 1/√2*( Y₁₋₁ - Y₁₁ )
    
    l = 2
    dxy      = y₂₋₂ = im/√2*( Y₂₋₂ - Y₂₂ )
    dyz      = y₂₋₁ = im/√2*( Y₂₋₁ + Y₂₁ )
    dz^2     = y₂₀  = Y₂₀
    dxz      = y₂₁  = 1/√2*( Y₂₋₁ - Y₂₁ )
    dx^2-y^2 = y₂₂  = 1/√2*( Y₂₋₂ + Y₂₂ )
"""
function Set_Comp2Real(l::Int)

    if l == 0
        C = [1.0+im*0.0;;]
        return C
    elseif l == 1
        # [px py pz] = [y₁₁ y₁₋₁ y₁₀] = C*[Y₁₋₁ Y₁₀ Y₁₁]
        C = [1/√2  0.0    -1/√2;
            im/√2  0.0    im/√2;
             0.0   1.0      0.0]
        return C
    elseif l == 2
        # [dz^2 dx^2-y^2 dxy dxz dyz] = [y₂₀ y₂₂ y₂₋₂ y₂₁ y₂₋₁] = C*[Y₂₋₂ Y₂₋₁ Y₂₀ Y₂₁ Y₂₂]
        C = [0.0     0.0     1.0     0.0       0.0;
             1/√2    0.0     0.0     0.0      1/√2;
             im/√2   0.0     0.0     0.0    -im/√2;
             0.0     1/√2    0.0    -1/√2      0.0;
             0.0    im/√2    0.0    im/√2      0.0]
        return C
    elseif l == 3
        C = [0.0     0.0     0.0     1.0     0.0       0.0       0.0;
             0.0     0.0     1/√2    0.0    -1/√2      0.0       0.0;
             0.0     0.0     im/√2   0.0    im/√2      0.0       0.0;
             0.0     1/√2    0.0     0.0     0.0       1/√2      0.0;
             0.0    im/√2    0.0     0.0     0.0      -im/√2     0.0;
             1/√2    0.0     0.0     0.0     0.0       0.0     -1/√2; 
             im/√2   0.0     0.0     0.0     0.0       0.0      im/√2]
        return C
    else
        return Calc_Comp2Real(l)
    end
end


function Set_Comp2Real!(C::Matrix{ComplexF64}, l::Int)
    if l == 0
        @. C = [1.0 + im*0.0;;]
    elseif l == 1
        # [px py pz] = [y₁₁ y₁₋₁ y₁₀] = C*[Y₁₋₁ Y₁₀ Y₁₁]
        @. C = [1/√2   0.0   -1/√2;
                im/√2  0.0   im/√2;
                0.0    1.0     0.0]
    elseif l == 2
        # [dz^2 dx^2-y^2 dxy dxz dyz] = [y₂₀ y₂₂ y₂₋₂ y₂₁ y₂₋₁] = C*[Y₂₋₂ Y₂₋₁ Y₂₀ Y₂₁ Y₂₂]
        @. C = [0.0     0.0     1.0     0.0       0.0;
                1/√2    0.0     0.0     0.0      1/√2;
                im/√2   0.0     0.0     0.0    -im/√2;
                0.0     1/√2    0.0    -1/√2      0.0;
                0.0    im/√2    0.0    im/√2      0.0]
    elseif l == 3
        @. C = [0.0     0.0     0.0     1.0     0.0       0.0       0.0;
             0.0     0.0     1/√2    0.0    -1/√2      0.0       0.0;
             0.0     0.0     im/√2   0.0    im/√2      0.0       0.0;
             0.0     1/√2    0.0     0.0     0.0       1/√2      0.0;
             0.0    im/√2    0.0     0.0     0.0      -im/√2     0.0;
             1/√2    0.0     0.0     0.0     0.0       0.0     -1/√2; 
             im/√2   0.0     0.0     0.0     0.0       0.0      im/√2]
    else
        Calc_Comp2Real!(C, l)
    end
end


function Calc_Comp2Real!(C::Matrix{ComplexF64}, L::Int)

    C[1,L+1] = 1.0 + 0.0*im
    inv_sqrt2 = 1/sqrt(2)

    j = -1
    for i = 1:4:2*L
        j += 1
        C[i+1,L-2*j]   = ComplexF64( inv_sqrt2, 0.0)
        C[i+1,L+2*j+2] = ComplexF64(-inv_sqrt2, 0.0)
    end

    j = 0
    for i = 3:4:2*L
        j += 1
        C[i+1,L-2*j+1] = ComplexF64( inv_sqrt2, 0.0)
        C[i+1,L+2*j+1] = ComplexF64( inv_sqrt2, 0.0)
    end

    j = -1
    for i = 2:4:2*L
        j += 1
        C[i+1,L-2*j] =   ComplexF64( 0.0, inv_sqrt2)
        C[i+1,L+2*j+2] = ComplexF64( 0.0, inv_sqrt2)
    end

    j = 0
    for i = 4:4:2*L
        j += 1
        C[i+1,L-2*j+1] = ComplexF64( 0.0,  inv_sqrt2)
        C[i+1,L+2*j+1] = ComplexF64( 0.0, -inv_sqrt2)
    end    
end


function Calc_Comp2Real(L::Int)
    C = zeros(ComplexF64, 2*L+1, 2*L+1)
    Calc_Comp2Real!(C, L)
    return C
end


function Set_Comp2Real!(C::Matrix{ComplexF64}, Spe_MaxL_Basis, Spe_Num_Basis)
    Sum = 1
    for l = 0:Spe_MaxL_Basis
        Cl = Set_Comp2Real(l)
        for _ = 1:Spe_Num_Basis[l+1]
            @. C[Sum:Sum+2*l, Sum:Sum+2*l] = Cl
            Sum += 2*l+1
        end
    end
end


function Set_NLComp2Real!(C::Matrix{ComplexF64}, VPS_List)
    
    MaxL = maximum(VPS_List)
    Nl = zeros(Int64, MaxL+1)
    for l = 0:MaxL
        Nl[l+1] = length(filter(x->x==l, VPS_List))
    end

    tot = 1
    for l = 0:MaxL, _ = 1:Nl[l+1]
        C[tot:tot+2*l, tot:tot+2*l] = Set_Comp2Real(l)
        tot += 2*l+1
    end
end


# for Set_ProExpn_VNA.jl
function Set_VNAComp2Real(maxL)
    
    C = Vector{Matrix{ComplexF64}}(undef, maxL+1)
    for l = 0:maxL
        C[l+1] = zeros(ComplexF64, 2*l+1, 2*l+1)
        Set_Comp2Real!(C[l+1], l)
    end

    return C
end

