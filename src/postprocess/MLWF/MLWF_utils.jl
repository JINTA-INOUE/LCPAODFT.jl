function Set_MLWF_Projector(MLWF_Orbs)

    MLWF_Num_Kinds_Projectors = length(MLWF_Orbs)

    MLWF_PAO_List1 = ["s","px","py","pz","dz2","dx2-y2","dxy","dxz","dyz","fz3","fxz2","fyz2","fzx2","fxyz","fx3-3xy2","f3yx2-y3"]
    MLWF_PAO_List2 = ["sp"]
    MLWF_PAO_List3 = ["p","sp2"]
    MLWF_PAO_List4 = ["sp3"]
    MLWF_PAO_List5 = ["sp3dz2", "d"]
    MLWF_PAO_List6 = ["sp3deg","sp2"]
    MLWF_PAO_List7 = ["f"]

    MLWF_Pro_lName = Vector{String}(undef, MLWF_Num_Kinds_Projectors)
    MLWF_ProNumBasis = zeros(Int32, MLWF_Num_Kinds_Projectors)
    MLWF_Symbol = Vector{String}(undef, MLWF_Num_Kinds_Projectors)
    MLWF_EachNum = zeros(Int32, MLWF_Num_Kinds_Projectors)
    MLWF_Func_Num = 0
    num = 0
    for orb in MLWF_Orbs

        num += 1
        mlwf_symbol, mlwf_pao = split(orb, "-")
        if length(mlwf_pao) <= 1
            @show MLWF_Orbs
            error("please check MLWF_Orbs")
        end


        mlwf_l = string(mlwf_pao[2:end])
        MLWF_Pro_lName[num] = mlwf_l
        MLWF_ProNumBasis[num] = parse(Int32, mlwf_pao[1])

        if mlwf_symbol[end] in ("S", "H")
            MLWF_Symbol[num] = mlwf_symbol[1:end-1]
            # mlwf_extrasymbol = mlwf_symbol[end:end]
        else
            MLWF_Symbol[num] = mlwf_symbol
        end
        @show mlwf_l

        if mlwf_l ∈ MLWF_PAO_List1
            MLWF_Func_Num += 1
            MLWF_EachNum[num] = 1
        elseif mlwf_l ∈ MLWF_PAO_List2
            MLWF_Func_Num += 2
            MLWF_EachNum[num] = 2
        elseif mlwf_l ∈ MLWF_PAO_List3
            MLWF_Func_Num += 3
            MLWF_EachNum[num] = 3
        elseif mlwf_l ∈ MLWF_PAO_List4
            MLWF_Func_Num += 4
            MLWF_EachNum[num] = 4
        elseif mlwf_l ∈ MLWF_PAO_List5
            MLWF_Func_Num += 5
            MLWF_EachNum[num] = 5
        elseif mlwf_l ∈ MLWF_PAO_List6
            MLWF_Func_Num += 6
            MLWF_EachNum[num] = 6
        elseif mlwf_l ∈ MLWF_PAO_List7
            MLWF_Func_Num += 7
            MLWF_EachNum[num] = 7
        else
            error("please check MLWF_Orbs")
        end
    end


    return MLWF_Pro_lName, MLWF_ProNumBasis, MLWF_Symbol, MLWF_EachNum, MLWF_Func_Num
end


function Set_MLWF_ProjName(MLWF_Num_Kinds_Projectors, MLWF_Pro_lName)
    
    MLWF_PAO_List1 = ["s","px","py","pz","dz2","dx2-y2","dxy","dxz","dyz","fz3","fxz2","fyz2","fzx2","fxyz","fx3-3xy2","f3yx2-y3"]
    MLWF_PAO_List2 = ["sp"]
    MLWF_PAO_List3 = ["p","sp2"]
    MLWF_PAO_List4 = ["sp3"]
    MLWF_PAO_List5 = ["sp3dz2", "d"]
    MLWF_PAO_List6 = ["sp3deg","sp2"]
    MLWF_PAO_List7 = ["f"]

    MLWF_ProjName = String[]
    for proj = 1:MLWF_Num_Kinds_Projectors

        mlwf_l = MLWF_Pro_lName[proj]

        if mlwf_l ∈ MLWF_PAO_List1
            push!(MLWF_ProjName, mlwf_l)
        elseif mlwf_l ∈ MLWF_PAO_List2
            push!(MLWF_ProjName, "s+px")
            push!(MLWF_ProjName, "s-px")
        elseif mlwf_l ∈ MLWF_PAO_List3
            if mlwf_l == "p"
                push!(MLWF_ProjName, "px")
                push!(MLWF_ProjName, "py")
                push!(MLWF_ProjName, "pz")
            elseif mlwf_l == "sp"
                push!(MLWF_ProjName, "s-px+py")
                push!(MLWF_ProjName, "s-px-py")
                push!(MLWF_ProjName, "s+px")
            end
        elseif mlwf_l ∈ MLWF_PAO_List4
            push!(MLWF_ProjName, "1/√2(s+px+py+pz)")
            push!(MLWF_ProjName, "1/√2(s+px-py-pz)")
            push!(MLWF_ProjName, "1/√2(s-px+py-pz)")
            push!(MLWF_ProjName, "1/√2(s-px-py+pz)")
        elseif mlwf_l ∈ MLWF_PAO_List5
            if mlwf_l == "d"
                push!(MLWF_ProjName, "dz2")
                push!(MLWF_ProjName, "dx2-y2")
                push!(MLWF_ProjName, "dxy")
                push!(MLWF_ProjName, "dxz")
                push!(MLWF_ProjName, "dyz")
            elseif mlwf_l == "sp3dz2"
                error("please add code.")
            end
        elseif mlwf_l ∈ MLWF_PAO_List6
            error("please add code.")
        elseif mlwf_l ∈ MLWF_PAO_List7
            push!(MLWF_ProjName, "fz3")
            push!(MLWF_ProjName, "fxz2")
            push!(MLWF_ProjName, "fyz2")
            push!(MLWF_ProjName, "fzx2")
            push!(MLWF_ProjName, "fxyz")
            push!(MLWF_ProjName, "fx3-3xy2")
            push!(MLWF_ProjName, "f3yx2-y3")
        else
            error("please check MLWF_Orbs")
        end
    end


    return MLWF_ProjName
end


function Set_MLWF_NumP_Pro(MLWF_Num_Kinds_Projectors, MLWF_ProNumBasis, MLWF_Pro_lName)

    MLWF_NumP_Pro = Vector{Vector{Int32}}(undef, MLWF_Num_Kinds_Projectors)

    for p = 1:MLWF_Num_Kinds_Projectors

        pnum = MLWF_ProNumBasis[p]
        SpeName = MLWF_Pro_lName[p]

        if SpeName == "s"
            MLWF_NumP_Pro[p] = [pnum, 1, 1, 1]
        elseif SpeName ∈ ("p", "px", "py", "pz")
            MLWF_NumP_Pro[p] = [1, pnum, 1, 1]
        elseif SpeName ∈ ("d", "dz2", "dx2-y2", "dxy", "dxz", "dyz")
            MLWF_NumP_Pro[p] = [1, 1, pnum, 1]
        elseif SpeName ∈ ("f", "fz3", "fxz2", "fyz2", "fzx2", "fxyz", "fx3-3xy2", "f3yx2-y3")
            MLWF_NumP_Pro[p] = [1, 1, 1, pnum]
        elseif SpeName ∈ ("sp", "sp2", "sp3")
            MLWF_NumP_Pro[p] = [pnum, pnum, 1, 1]
        elseif SpeName ∈ ("sp3dz2", "sp3deg")
            MLWF_NumP_Pro[p] = [pnum, pnum, pnum, 1]
        else
            @show SpeName
            error("not find SpeName in Set_MLWF_NumP_Pro. pleace check input files")
        end
    end


    return MLWF_NumP_Pro
end


function Analyze_Wannier_Projectors(MLWF_Num_Kinds_Projectors, MLWF_Pro_lName)

    MLWF_Num_Pro = zeros(Int32, MLWF_Num_Kinds_Projectors)
    MLWF_NumL_Pro = Vector{Vector{Int32}}(undef, MLWF_Num_Kinds_Projectors)
    MLWF_Select_Matrix = Vector{Vector{Int32}}(undef, MLWF_Num_Kinds_Projectors)
    MLWF_Projector_Hybridize_Matrix = Vector{Matrix{Float64}}(undef, MLWF_Num_Kinds_Projectors)

    for p = 1:MLWF_Num_Kinds_Projectors

        SpeName = MLWF_Pro_lName[p]

        if SpeName == "s"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [1, 0, 0, 0]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "px"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "py"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "pz"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "dz2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "dx2-y2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [1]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "dxy"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [2]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "dxz"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [3]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "dyz"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [4]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fz3"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [0]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fxz2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [1]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fyz2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [2]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fzx2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [3]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fxyz"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [4]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "fx3-3xy2"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [5]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "f3yx2-y3"
            MLWF_Num_Pro[p] = 1
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [6]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0;;]
        elseif SpeName == "sp"
            MLWF_Num_Pro[p] = 2
            MLWF_NumL_Pro[p] = [1, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0, 1]
            MLWF_Projector_Hybridize_Matrix[p] = [1 1; 1 -1]/sqrt(2)
        elseif SpeName == "sp2"
            MLWF_Num_Pro[p] = 3
            MLWF_NumL_Pro[p] = [1, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2]
            MLWF_Projector_Hybridize_Matrix[p] = [1/sqrt(3) -1/sqrt(6) 1/sqrt(2); 1/sqrt(3) -1/sqrt(6) -1/sqrt(2); 1/sqrt(3) 2/sqrt(6) 0.0]
        elseif SpeName == "sp3"
            MLWF_Num_Pro[p] = 4
            MLWF_NumL_Pro[p] = [1, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2, 3]
            MLWF_Projector_Hybridize_Matrix[p] = [1 1 1 1; 1 1 -1 -1; 1 -1 1 -1; 1 -1 -1 1]/2
        elseif SpeName == "sp3dz2"
            MLWF_Num_Pro[p] = 5
            MLWF_NumL_Pro[p] = [1, 1, 1, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2, 3, 4]
            MLWF_Projector_Hybridize_Matrix[p] = [1/sqrt(3) -1/sqrt(6) 1/sqrt(2) 0 0; 
                                                  1/sqrt(3) -1/sqrt(6) -1/sqrt(2) 0 0; 
                                                  1/sqrt(3) 2/sqrt(6) 0 0 0; 
                                                  0 0 0 1/sqrt(2) 1/sqrt(2); 
                                                  0 0 0 -1/sqrt(2) 1/sqrt(2)]
        elseif SpeName == "sp3deg"
            MLWF_Num_Pro[p] = 5
            MLWF_NumL_Pro[p] = [1, 1, 1, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2, 3, 4, 5]
            MLWF_Projector_Hybridize_Matrix[p] = [1/sqrt(6) -1/sqrt(2) 0.0 0.0 -1/sqrt(12) 0.5;
                                                1/sqrt(6)  1/sqrt(2) 0.0 0.0 -1/sqrt(12) 0.5;
                                                1/sqrt(6)  0.0  -1/sqrt(2) 0.0 -1/sqrt(12) -0.5;
                                                1/sqrt(6)  0.0   1/sqrt(2) 0.0 -1/sqrt(12) -0.5;
                                                1/sqrt(6)  0.0   0.0 -1/sqrt(2) 1/sqrt(3) 0.0;
                                                1/sqrt(6)  0.0   0.0  1/sqrt(2) 1/sqrt(3) 0.0;]
        elseif SpeName == "p"
            MLWF_Num_Pro[p] = 3
            MLWF_NumL_Pro[p] = [0, 1, 0, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0 0.0 0.0; 
                                                0.0 1.0 0.0; 
                                                0.0 0.0 1.0]
        elseif SpeName == "d"
            MLWF_Num_Pro[p] = 5
            MLWF_NumL_Pro[p] = [0, 0, 1, 0]
            MLWF_Select_Matrix[p] = [0, 1, 2, 3, 4]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0 0.0 0.0 0.0 0.0; 
                                                0.0 1.0 0.0 0.0 0.0; 
                                                0.0 0.0 1.0 0.0 0.0; 
                                                0.0 0.0 0.0 1.0 0.0; 
                                                0.0 0.0 0.0 0.0 1.0]
        elseif SpeName == "f"
            MLWF_Num_Pro[p] = 7
            MLWF_NumL_Pro[p] = [0, 0, 0, 1]
            MLWF_Select_Matrix[p] = [0, 1, 2, 3, 4, 5, 6]
            MLWF_Projector_Hybridize_Matrix[p] = [1.0 0.0 0.0 0.0 0.0 0.0 0.0; 
                                                0.0 1.0 0.0 0.0 0.0 0.0 0.0; 
                                                0.0 0.0 1.0 0.0 0.0 0.0 0.0; 
                                                0.0 0.0 0.0 1.0 0.0 0.0 0.0; 
                                                0.0 0.0 0.0 0.0 1.0 0.0 0.0;
                                                0.0 0.0 0.0 0.0 0.0 1.0 0.0;
                                                0.0 0.0 0.0 0.0 0.0 0.0 1.0]
        else
            error("pleace check input files")
        end
    end

    return MLWF_Num_Pro, MLWF_NumL_Pro, MLWF_Select_Matrix, MLWF_Projector_Hybridize_Matrix
end


function Get_Wannier_Euler_Rotation_Angle(MLWF_Num_Kinds_Projectors, MLWF_Z_Direction, MLWF_X_Direction)

    MLWF_Euler_Rotation_Angle = zeros(Float64, MLWF_Num_Kinds_Projectors, 3)
    for p = 1:MLWF_Num_Kinds_Projectors

        Norm = norm(MLWF_Z_Direction[p,:])
        zx = MLWF_Z_Direction[p,1]/Norm
        zy = MLWF_Z_Direction[p,2]/Norm
        zz = MLWF_Z_Direction[p,3]/Norm

        Norm = norm(MLWF_X_Direction[p,:])
        xx = MLWF_X_Direction[p,1]/Norm
        xy = MLWF_X_Direction[p,2]/Norm
        xz = MLWF_X_Direction[p,3]/Norm

        coszx = zx*xx+zy*xy+zz*xz
        if abs(coszx) > 1e-6
            error("please check")
        end

        yx=zy*xz-zz*xy
        yy=zz*xx-zx*xz
        yz=zx*xy-zy*xx

        Norm=sqrt(yx*yx+yy*yy+yz*yz)
        yx=yx/Norm
        yy=yy/Norm
        yz=yz/Norm

        if abs(abs(zz) - 1.0) < 1e-6
            beta = ifelse(zz > 0.0, 0.0, pi)
            gamma = 0.0
            alpha = asin(xy/cos(beta))
            if xx/zz < 0.0
                alpha = pi - alpha
            end
            if alpha < 0.0
                alpha = 2*pi + alpha
            end
        else
            beta = acos(zz)
            alpha = asin(zy/sqrt(zx*zx + zy*zy))
            if zx < 0.0
                alpha = pi - alpha
            end
            if alpha < 0.0
                alpha = 2*pi + alpha
            end

            if abs(abs(-xz/sin(beta)) - 1.0) < 1e-5
                gamma = ifelse(-xz/sin(beta)<0.0, pi, 0.0)
            else
                gamma = acos(-xz/sin(beta))
            end

            tmp1 = -(xx-(-xz/sin(beta))*zz*cos(alpha))/sin(alpha)
            if tmp1 < 0.0
                gamma = 2*pi - gamma
            end
        end

        MLWF_Euler_Rotation_Angle[p,1] = alpha
        MLWF_Euler_Rotation_Angle[p,2] = beta
        MLWF_Euler_Rotation_Angle[p,3] = gamma

        println("x-axis, $xx, $xy, $xz")
        println("y-axis, $yx, $yy, $yz")
        println("z-axis $zx, $zy, $zz")
        println("Euler Angles are $(alpha/pi*180.0), $(beta/pi*180.0), $(gamma/pi*180.0) .(in degree)")
        println("Euler Angles are $alpha, $beta, $gamma.(in rad)\n")
    end


    return MLWF_Euler_Rotation_Angle
end


function Get_Rotational_Matrix!(L, Euler_Rotation_Angle, RotMat_for_Real_Func)
    
    j = L
    alpha = Euler_Rotation_Angle[1]
    beta = Euler_Rotation_Angle[2]
    gamma = Euler_Rotation_Angle[3]

    dj = zeros(Float64, 2*L+1, 2*L+1)
    Dlm = zeros(ComplexF64, 2*L+1, 2*L+1)
    for i = 1:2*L+1
        Dlm[i,i] = 1.0 + 0.0*im
        RotMat_for_Real_Func[i,i] = 1.0+0.0*im
    end


    for jmp = 0:2*j, jm = jmp:2*j
        mp = j - jmp
        m = j - jm

        tmp1 = factorial(j+m)*factorial(j-m)*factorial(j+mp)*factorial(j-mp)
        fac1 = sqrt(tmp1)
        maxk = min(j+m, j-mp)

        Sumk = 0.0
        for k = 0:maxk
            tmp1 = factorial(j-mp-k)*factorial(j+m-k)*factorial(k+mp-m)*factorial(k)
            tmp2 = ifelse(abs(cos(beta/2))<1e-8 && iszero(2*j+m-mp-2*k), 1.0, exp(log(abs(cos(beta/2)))*(2*j+m-mp-2*k)))
            tmp3 = ifelse(abs(sin(beta/2))<1e-8 && iszero(2*k+mp-m), 1.0, exp(log(abs(sin(beta/2)))*(2*k+mp-m)))
            fac2 = tmp2*tmp3/tmp1
            Sumk = Sumk + ifelse(iszero(mod(k+mp-m,2)), fac2, -fac2)
        end

        dj[jmp+1,jm+1] = Sumk*fac1
        Dlm[jmp+1,jm+1] = Sumk*fac1*Complex(cos(mp*alpha+m*gamma), -sin(m*gamma+mp*alpha))
    end


    for jmp = 0:2*j, jm = 0:jmp-1
        mp = j-jmp
        m = j-jm
        if iszero(mod(mp-m,2))
            dj[jmp+1,jm+1] = dj[jm+1,jmp+1]
        else
            dj[jmp+1,jm+1] = -dj[jm+1,jmp+1]
        end

        tmp1 =  cos(mp*alpha)*cos(m*gamma) - sin(mp*alpha)*sin(m*gamma)
        tmp2 = -cos(mp*alpha)*sin(m*gamma) - sin(mp*alpha)*cos(m*gamma)

        Dlm[jmp+1,jm+1] = dj[jmp+1,jm+1] * Complex(tmp1, tmp2)
    end

    

    #=
    The rotation matrix connecting real function is defined as M=U^(-1)*Dlm^(T)*U, U is the transfer matrix
    from real to imaginary function for orbitals
     p1         px
     p0   =  U* py
     p-1        pz

     d2           dz2
     d1           dx2-y2
     d0    =  U*  dxy
     d-1          dxz
     d-2          dyz

     f3           fz3
     f2           fz3
     f1           fz3
     f0    =   U* fz3
     f-1          fz3
     f-2          fz3
     f-3          fz3

     px'         px
     py'  =  M * py
     pz'         pz

     dz2'          dz2
     dx2-y2'       dx2-y2
     dxy'    = M * dxy
     dxz'          dxz
     dyz'          dyz
    Here ' means those in the rotated coordinate, without ' means those in original coordinate.
    =#
    if j == 1
        Umat = [-1/sqrt(2) -im/sqrt(2) 0.0;
                  0.0         0.0      1.0;
                 1/sqrt(2) -im/sqrt(2) 0.0]
    elseif j == 2
        Umat = [  0.0   1/sqrt(2)  im/sqrt(2)  0.0           0.0;
                  0.0      0.0       0.0    -1/sqrt(2) -im/sqrt(2);
                  1.0      0.0       0.0       0.0           0.0;
                  0.0      0.0       0.0     1/sqrt(2) -im/sqrt(2);
                  0.0   1/sqrt(2)  -im/sqrt(2) 0.0          0.0]
    elseif j == 3
        Umat = [  0.0   1/sqrt(2)  im/sqrt(2)  0.0           0.0;
                  0.0      0.0       0.0    -1/sqrt(2) -im/sqrt(2);
                  1.0      0.0       0.0       0.0           0.0;
                  0.0      0.0       0.0     1/sqrt(2) -im/sqrt(2);
                  0.0   1/sqrt(2)  -im/sqrt(2) 0.0          0.0]
    end


    Utmp = zeros(ComplexF64, 2*j+1, 2*j+1)
    Umat_inv = zeros(ComplexF64, 2*j+1, 2*j+1)
    for jmp = 1:2*j+1, jm = 1:2*j+1
        Umat_inv[jmp,jm] = conj(Umat[jm,jmp])
    end

    for jmp = 1:2*j+1, jm = 1:2*j+1
        Sum = 0.0 + im*0.0
        for i = 1:2*j+1
            Sum += Umat_inv[jmp,i]*Dlm[jm,i]
        end
        Utmp[jmp,jm] = Sum
    end

    for jmp = 1:2*j+1, jm = 1:2*j+1
        Sum = 0.0 + im*0.0
        for i = 1:2*j+1
            Sum += Utmp[jmp,i]*Umat[i,jm]
        end
        Umat_inv[jmp,jm] = Sum
    end


    for jmp = 1:2*j+1, jm = 1:2*j+1
        if abs(imag(Umat_inv[jmp,jm])) > 1e-8
            error("please check rotational matrix")
        end
        @. RotMat_for_Real_Func = real(Umat_inv)
    end
end


function Find_NN_Projectors_Basis(MLWF_Atom_Cut1, MLWF_Pos, TCpyCell, Natom, Gxyz, atv, Atoms_Cut1)

    MLWF_Num_Kinds_Projectors = length(MLWF_Atom_Cut1)

    FNAN_WP = zeros(Int32, MLWF_Num_Kinds_Projectors)
    for atom = 1:MLWF_Num_Kinds_Projectors

        rcutA = MLWF_Atom_Cut1[atom]
        FNAN_WP[atom] = 0

        for Rn = 0:TCpyCell
            for jatom = 1:Natom

                rcutB = Atoms_Cut1[jatom]

                dx = MLWF_Pos[atom][1] - Gxyz[jatom][1] - atv[Rn+1][1]
                dy = MLWF_Pos[atom][2] - Gxyz[jatom][2] - atv[Rn+1][2]
                dz = MLWF_Pos[atom][3] - Gxyz[jatom][3] - atv[Rn+1][3]
                
                r = sqrt(abs(dx^2 + dy^2 + dz^2))
                if r <= rcutA + rcutB
                    FNAN_WP[atom] = FNAN_WP[atom] + 1
                end
            end
        end
    end


    natn_WP = Vector{Vector{Int32}}(undef, MLWF_Num_Kinds_Projectors)
    ncn_WP = Vector{Vector{Int32}}(undef, MLWF_Num_Kinds_Projectors)
    for atom = 1:MLWF_Num_Kinds_Projectors
        natn_WP[atom] = zeros(Int32, FNAN_WP[atom])
        ncn_WP[atom] = zeros(Int32, FNAN_WP[atom])
    end


    for atom = 1:MLWF_Num_Kinds_Projectors

        rcutA = MLWF_Atom_Cut1[atom]
        FNAN_WP[atom] = 0

        for Rn = 0:TCpyCell
            for jatom = 1:Natom

                rcutB = Atoms_Cut1[jatom]

                dx = MLWF_Pos[atom][1] - Gxyz[jatom][1] - atv[Rn+1][1]
                dy = MLWF_Pos[atom][2] - Gxyz[jatom][2] - atv[Rn+1][2]
                dz = MLWF_Pos[atom][3] - Gxyz[jatom][3] - atv[Rn+1][3]
                
                r = sqrt(abs(dx^2 + dy^2 + dz^2))
                if r <= rcutA + rcutB
                    FNAN_WP[atom] = FNAN_WP[atom]+1
                    natn_WP[atom][FNAN_WP[atom]] = jatom
                    ncn_WP[atom][FNAN_WP[atom]] = Rn
                end
            end
        end
    end


    return FNAN_WP, natn_WP, ncn_WP
end
