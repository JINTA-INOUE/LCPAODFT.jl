function Set_MLWF_kgrid(Latvecs::Matrix{Float64}, kmesh, MAXSHELL)

    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    Recvecs = 2*pi*inv(Latvecs')
    Recvecs_len = zeros(Float16, 3)
    Recvecs_len1 = sqrt(Recvecs[1,1]^2 + Recvecs[1,2]^2 + Recvecs[1,3]^2)
    Recvecs_len2 = sqrt(Recvecs[2,1]^2 + Recvecs[2,2]^2 + Recvecs[2,3]^2)
    Recvecs_len3 = sqrt(Recvecs[3,1]^2 + Recvecs[3,2]^2 + Recvecs[3,3]^2)


    metric = zeros(Float64, 3, 3)
    for j = 1:3, i = 1:j
        Sumr = 0.0
        for l = 1:3
            Sumr += Latvecs[i,l]*Latvecs[j,l]
        end
        metric[i,j] = Sumr
        if i < j
            metric[j,i] = metric[i,j]
        end
    end


    rvect = Vector{Vector{Int32}}(undef, 1)
    for cell = 1:1
        rvect[cell] = zeros(Int32, 3)
    end
    ndegen = zeros(Int32, 1)
    NCell = Wigner_Seitz_Vectors!(metric, kmesh, -1, rvect, ndegen)
    rvect = Vector{Vector{Int32}}(undef, NCell)
    for cell = 1:NCell
        rvect[cell] = zeros(Int32, 3)
    end
    ndegen = zeros(Int32, NCell)
    NCell = Wigner_Seitz_Vectors!(metric, kmesh, 0, rvect, ndegen)


    kg = Vector{Vector{Float64}}(undef, Nkpt)
    ik = 1
    for i = kmesh1:-1:1, j = kmesh2:-1:1, k = kmesh3:-1:1
        kg[ik] = zeros(Float64, 3)
        kg[ik][1] = (kmesh[1]-i)/kmesh1 + 1e-12
        kg[ik][2] = (kmesh[2]-j)/kmesh2 + 1e-12
        kg[ik][3] = (kmesh[3]-k)/kmesh3 + 1e-12
        ik += 1
    end


    klatt = Vector{Vector{Float64}}(undef, 3)
    for i = 1:3
        klatt[i] = zeros(Float64, 3)
    end

    for i = 1:3
        klatt[1][i] = Recvecs[1,i]/kmesh1
        klatt[2][i] = Recvecs[2,i]/kmesh2
        klatt[3][i] = Recvecs[3,i]/kmesh3
    end

    

    tmp_M_s = zeros(Int32, MAXSHELL)
    tmp_bvector = Vector{Vector{Int32}}(undef, 8*MAXSHELL^3)
    for i = 1:8*MAXSHELL^3
        tmp_bvector[i] = zeros(Int32, 3)
    end

    shell_num = Shell_Structure!(klatt, tmp_M_s, tmp_bvector, MAXSHELL)
    if shell_num == 0
        printf("******************************Error********************************")
        println("*    Can not find proper b vectors, please increase parameter     *")
        println("*    MAXSHELL OR change Wannier.Kgrids.                           *")
        println("******************************Error********************************")
        println("***********************************INFO**************************************")
        @printf("Reciprocal Lattices lengths are:%10.6f %10.6f %10.6f\n", Recvecs_len1, Recvecs_len2, Recvecs_len3)
        @printf("The ratio among them are: b1:b2=%10.6f b1:b3=%10.6f b2:b3=%10.6f\n", Recvecs_len1/Recvecs_len2, Recvecs_len1/Recvecs_len3, Recvecs_len2/Recvecs_len3)
        println("Message: Please try to set Wannier.Kgrid has the similar ratio as above.")
        println("************************************INFO*************************************")
    end



    Reject_Shell = zeros(Int32, shell_num)
    tmp_wb = zeros(Float64, shell_num)
    find_w, shell_num, searched_shell = Cal_Weight_of_Shell!(klatt, tmp_M_s, tmp_bvector, shell_num, tmp_wb, Reject_Shell)
    if find_w == 0
        println("*************************** Error ****************************")
        println("*    Weights for b vectors (totally $shell_num) are not found.      *")
        println("*    Please increase MAXSHELL (presently it is $MAXSHELL) OR       *")
        println("*    change Wannier.Kgrids and try again                     *")
        println("*************************** Error ****************************")
        println("***********************************INFO**************************************")
        @printf("Reciprocal Lattices lengths are:%10.6f %10.6f %10.6f\n", Recvecs_len1, Recvecs_len2, Recvecs_len3)
        @printf("The ratio among them are: b1:b2=%10.6f b1:b3=%10.6f b2:b3=%10.6f\n", Recvecs_len1/Recvecs_len2, Recvecs_len1/Recvecs_len3, Recvecs_len2/Recvecs_len3)
        println("Message: Please try to set Wannier.Kgrid has the similar ratio as above.")
        println("************************************INFO*************************************")
        error("stop")
    else
        tot_bvector, bvector, frac_bv_int, frac_bv, wb = Set_bvectors(Recvecs, kmesh, shell_num, searched_shell, Reject_Shell, tmp_bvector, tmp_M_s, tmp_wb)
    end

    tot_bvector -= 1
    kplusb = Set_kplusb(Nkpt, tot_bvector, kg, frac_bv)


    return kg, tot_bvector, bvector, frac_bv_int, frac_bv, kplusb, wb
end


function Wigner_Seitz_Vectors!(metric, kmesh, r_num, rvect, ndegen)

    nd = zeros(Int32, 125, 3)
    dist = zeros(Float64, 125)

    rnum = 0
    for n1 = -kmesh[1]:kmesh[1], n2 = -kmesh[2]:kmesh[2], n3 = -kmesh[3]:kmesh[3]
        icnt = 1
        for i1 = -2:2, i2 = -2:2, i3 = -2:2
            nd[icnt,1] = n1 - i1*kmesh[1]
            nd[icnt,2] = n2 - i2*kmesh[2]
            nd[icnt,3] = n3 - i3*kmesh[3]
            dist[icnt] = 0.0
            for i = 1:3, j = 1:3
                dist[icnt] += nd[icnt,i]*metric[i,j]*nd[icnt,j]
            end

            icnt += 1
        end

        dist_min = minimum(dist)

        if abs(dist[63] - dist_min) < 0.0000001
            if r_num ≠ -1
                ndegen[rnum+1] = 0
                for i = 1:125
                    if abs(dist[i] - dist_min) < 0.0000001
                        ndegen[rnum+1] += 1
                    end
                end

                rvect[rnum+1][1] = n1
                rvect[rnum+1][2] = n2
                rvect[rnum+1][3] = n3
            end
            rnum += 1
        end
    end

    if r_num == -1
        return rnum
    end


    println("There are $rnum lattice points found in Wigner-Seitz supercell.")

    dist_min = 0.0
    for i = 1:rnum
        dist_min += 1/ndegen[i]
    end

    if abs(dist_min-prod(kmesh))>1e-6
        println("**************************** Error **********************");
        println("*   In Wigner_Seitz_Vectors subroutine, error happens.  *");
        println("*   Please change setting of Wannier.Kgrid.             *");
        println("**************************** Error **********************");
        println("***********************************INFO**************************************");
        # println("Reciprocal Lattices lengths are : $(Recvecs[1,1])  $(Recvecs[1,2])  $(Recvecs[1,3])")
        # printf("The ratio among them are: b1:b2=%10.6f b1:b3=%10.6f b2:b3=%10.6f\n",rtv[0][0]/rtv[0][1],rtv[0][0]/rtv[0][2],rtv[0][1]/rtv[0][2]);
        println("Message: Please try to set Wannier.Kgrid has the similar ratio as above.");
        println("************************************INFO*************************************");
    end


    return rnum
end


function Ascend_Ordering!(xyz_value, ordering, tot_kpt)

    for i = 1:tot_kpt-1, j = i:-1:1
        if xyz_value[j+1] < xyz_value[j]
            tmp_xyz = xyz_value[j+1]
            xyz_value[j+1] = xyz_value[j]
            xyz_value[j] = tmp_xyz
            tmp_order = ordering[j+1]
            ordering[j+1] = ordering[j]
            ordering[j] = tmp_order
        end
    end
end


function Shell_Structure!(klatt, M_s, bvector, MAXSHELL)

    kindx = 0
    combination = zeros(Int32, 8*MAXSHELL^3, 3)
    ordered_com = zeros(Int32, 8*MAXSHELL^3, 3)
    distance = zeros(Float64, 8*MAXSHELL^3)
    ordering = zeros(Int32, 8*MAXSHELL^3)

    for i1 = -MAXSHELL+1:MAXSHELL-1, i2 = -MAXSHELL+1:MAXSHELL-1, i3= -MAXSHELL+1:MAXSHELL-1
        
        dx = i1*klatt[1][1] + i2*klatt[2][1] + i3*klatt[3][1]
        dy = i1*klatt[1][2] + i2*klatt[2][2] + i3*klatt[3][2]
        dz = i1*klatt[1][3] + i2*klatt[2][3] + i3*klatt[3][3]

        distance[kindx+1] = sqrt(abs(dx*dx + dy*dy + dz*dz))
        combination[kindx+1,1] = i1
        combination[kindx+1,2] = i2
        combination[kindx+1,3] = i3
        ordering[kindx+1] = kindx

        kindx += 1
    end

    tot_kpt = kindx
    Ascend_Ordering!(distance, ordering, tot_kpt)

    for kindx = 1:tot_kpt
        ordered_com[kindx,1] = combination[ordering[kindx]+1,1]
        ordered_com[kindx,2] = combination[ordering[kindx]+1,2]
        ordered_com[kindx,3] = combination[ordering[kindx]+1,3]
    end


    i = 0
    dx = distance[1]
    current_shell = 0
    for kindx = 1:tot_kpt
        if abs(dx-distance[kindx]) > 1e-5
            ordering[current_shell+1] = i
            i = 1
            current_shell += 1
            dx = distance[kindx]
        else
            i += 1
        end
    end

    tot_shell = ifelse(current_shell-1 > MAXSHELL, MAXSHELL, 0)
    
    tot_vectors = 0
    for i = 1:tot_shell
        M_s[i] = ordering[i+1]
        tot_vectors += M_s[i]
    end

    for i = 1:tot_vectors, j = 1:3
        bvector[i][j] = ordered_com[i+1,j]
    end

    return tot_shell
end


function Cal_Weight_of_Shell!(klatt, M_s, bvector, num_shell, wb, Reject_Shell)

    shell_num = num_shell

    qvector = zeros(Float64, 6)
    qvector[1] = 1.0
    qvector[4] = 1.0
    qvector[6] = 1.0

    for i = 1:shell_num
        Reject_Shell[i] = 0.0
    end

    find_w = 0
    searced_shell = 0
    current_shell = 1
    shell_size = 1

    while find_w == 0 && current_shell <= shell_num

        Amatrix = zeros(Float64, 6, shell_size)
        copy_Amatrix = zeros(Float64, 6, shell_size)
        Dmatrix = zeros(Float64, 6, shell_size)
        Umatrix = zeros(Float64, 6, 6)
        VTmatrix = zeros(Float64, shell_size, shell_size)

        if shell_size > 6
            dsing = zeros(Float64, 6)
        else
            dsing = zeros(Float64, shell_size)
        end


        for i = 1:shell_size
            wb[i] = 0.0
        end

        realshellindx = -1
        for shellindx = 0:current_shell-1
            if Reject_Shell[shellindx+1] == 1
                println("Shell $(shellindx+1) is rejected.")
                continue
            end

            realshellindx += 1
            for i = 1:6
                Amatrix[i,realshellindx+1] = 0.0
            end

            startbv = 0
            if shellindx == 0
                startbv = 0
            else
                for i = 1:shellindx
                    startbv = startbv + M_s[i]
                end
            end

            for bvindx = startbv:M_s[shellindx+1]+startbv-1
                bx = bvector[bvindx+1][1]*klatt[1][1] + bvector[bvindx+1][2]*klatt[2][1] + bvector[bvindx+1][3]*klatt[3][1]
                by = bvector[bvindx+1][1]*klatt[1][2] + bvector[bvindx+1][2]*klatt[2][2] + bvector[bvindx+1][3]*klatt[3][2]
                bz = bvector[bvindx+1][1]*klatt[1][3] + bvector[bvindx+1][2]*klatt[2][3] + bvector[bvindx+1][3]*klatt[3][3]

                Amatrix[1,realshellindx+1] += bx*bx
                Amatrix[2,realshellindx+1] += by*bx
                Amatrix[3,realshellindx+1] += bz*bx
                Amatrix[4,realshellindx+1] += by*by
                Amatrix[5,realshellindx+1] += bz*by
                Amatrix[6,realshellindx+1] += bz*bz
            end

            copy_Amatrix[1,realshellindx+1] = Amatrix[1,realshellindx+1]
            copy_Amatrix[2,realshellindx+1] = Amatrix[2,realshellindx+1]
            copy_Amatrix[3,realshellindx+1] = Amatrix[3,realshellindx+1]
            copy_Amatrix[4,realshellindx+1] = Amatrix[4,realshellindx+1]
            copy_Amatrix[5,realshellindx+1] = Amatrix[5,realshellindx+1]
            copy_Amatrix[6,realshellindx+1] = Amatrix[6,realshellindx+1]
        end


        # Umat, sigma, Vmat = svd(Amatrix)
        Umat, sigma, Vmat = svd(Amatrix; full = true)

        fill!(Umatrix, 0.0)
        fill!(Dmatrix, 0.0)
        fill!(VTmatrix, 0.0)
    

        if shell_size < 6
            for i = 1:shell_size
                if abs(sigma[i]) < 1.0e-5
                    Reject_Shell[current_shell] = 1
                    break
                end
            end
        else
            for i = 1:6
                if abs(sigma[i]) < 1.0e-5
                    Reject_Shell[current_shell] = 1
                    break
                end
            end
        end


        if Reject_Shell[current_shell] == 0

            if shell_size < 6
                for i = 1:shell_size
                    Dmatrix[i,i] = 1/sigma[i]
                end
            else
                for i = 1:6
                    Dmatrix[i,i] = 1/sigma[i]
                end
            end

            for i = 1:6, j = 1:6
                Umatrix[i,j] = Umat[i,j]
            end
            

            for i = 1:shell_size, j = 1:shell_size
                # VTmatrix[i,j] = Vmat[i,j]
                VTmatrix[i,j] = Vmat[j,i]
            end

            

            for i = 1:shell_size
                for j = 1:6
                    Sum = 0.0
                    for k = 1:shell_size
                        Sum += VTmatrix[k,i]*Dmatrix[j,k]
                    end
                    Amatrix[j,i] = Sum
                end
            end

            for i = 1:shell_size
                for j = 1:6
                    Sum = 0.0
                    for k = 1:6
                        Sum += Amatrix[k,i]*Umatrix[j,k]
                    end
                    Dmatrix[j,i] = Sum
                end
            end

            for i = 1:shell_size
                Sum = 0.0
                for k = 1:6
                    Sum += Dmatrix[k,i]*qvector[k]
                end
                wb[i] = Sum
            end
        

            flag = 0
            for j = 1:6
                Sum = 0.0
                for i = 1:shell_size
                    Sum += wb[i]*copy_Amatrix[j,i]
                end

                if abs(Sum-qvector[j]) > 1.0e-5
                    flag = 1
                end
            end

            if flag == 0
                find_w = 1
                num_shell = shell_size
                searced_shell = current_shell
            else
                shell_size += 1
            end
        end

        current_shell += 1
    end


    return find_w, num_shell, searced_shell
end


function Set_bvectors(Recvecs, kmesh, shell_num, searched_shell, Reject_Shell, tmp_bvector, tmp_M_s, tmp_wb)

    tot_bvector = 0
    M_s = zeros(Int32, shell_num)
    j = 1
    for i = 1:searched_shell
        if Reject_Shell[i] == 0
            M_s[j] = tmp_M_s[i]
            tot_bvector += M_s[j]
            j += 1
        end
    end

    wb = zeros(Float64, tot_bvector)
    k = 1
    for i = 1:shell_num, j = 1:M_s[i]
        wb[k] = tmp_wb[i]
        k += 1
    end

    frac_bv_int = Vector{Vector{Int64}}(undef, tot_bvector)
    bvector = Vector{Vector{Float64}}(undef, tot_bvector)
    for ib = 1:tot_bvector
        frac_bv_int[ib] = zeros(Int64, 3)
        bvector[ib] = zeros(Float64, 3)
    end

    k = 1
    for i = 0:searched_shell-1
        if Reject_Shell[i+1] == 0
            startbv = 0
            if i == 0
                startbv = 0
            else
                for j = 0:i-1
                    startbv += tmp_M_s[j+1]
                end
            end

            for j = startbv:startbv+tmp_M_s[i+1]-1
                frac_bv_int[k][1] = tmp_bvector[j+1][1]
                frac_bv_int[k][2] = tmp_bvector[j+1][2]
                frac_bv_int[k][3] = tmp_bvector[j+1][3]
                bvector[k][1] = tmp_bvector[j+1][1]/kmesh[1]
                bvector[k][2] = tmp_bvector[j+1][2]/kmesh[2]
                bvector[k][3] = tmp_bvector[j+1][3]/kmesh[3]
                k += 1
            end
        end
    end

    
    println("There are $shell_num shells and total number of b vectors is $tot_bvector")

    frac_bv = Vector{Vector{Float64}}(undef, tot_bvector)
    for ib = 1:tot_bvector
        frac_bv[ib] = zeros(Float64, 3)
    end

    tot_bvector = 1
    wbtot = 0.0
    for i = 1:shell_num
        println("Shell $i has $(M_s[i]) b vectors:")
        println("No.|      Fractional Coordinate     || Cartesian Coordinate (Angs^-1)||Weight_b(Angs^2)||")

        for j = 1:M_s[i]
            frac_bv[tot_bvector][1] = bvector[tot_bvector][1]
            frac_bv[tot_bvector][2] = bvector[tot_bvector][2]
            frac_bv[tot_bvector][3] = bvector[tot_bvector][3]
            dkx = bvector[tot_bvector][1]*Recvecs[1,1] + bvector[tot_bvector][2]*Recvecs[2,1] + bvector[tot_bvector][3]*Recvecs[3,1]
            dky = bvector[tot_bvector][1]*Recvecs[1,2] + bvector[tot_bvector][2]*Recvecs[2,2] + bvector[tot_bvector][3]*Recvecs[3,2]
            dkz = bvector[tot_bvector][1]*Recvecs[1,3] + bvector[tot_bvector][2]*Recvecs[2,3] + bvector[tot_bvector][3]*Recvecs[3,3]

            @printf(" %2d|  (%8.5f,%8.5f,%8.5f)  ||  (%8.5f,%8.5f,%8.5f) || %13.5f  ||\n",
                    j, bvector[tot_bvector][1],bvector[tot_bvector][2],bvector[tot_bvector][3],
                    dkx*Ang_to_bohr,dky*Ang_to_bohr,dkz*Ang_to_bohr,
                    wb[tot_bvector]/Ang_to_bohr/Ang_to_bohr)

            bvector[tot_bvector][1] = dkx
            bvector[tot_bvector][2] = dky
            bvector[tot_bvector][3] = dkz
            wbtot += wb[tot_bvector]
            tot_bvector += 1
        end
    end


    return tot_bvector, bvector, frac_bv_int, frac_bv, wb
end


function Set_bdirection(tot_bvector, frac_bv)

    bdirection = zeros(Int32, tot_bvector)
    for i = 1:tot_bvector
        bdirection[i] = -1
    end
    nbdir = 0

    for bindx = 1:tot_bvector
        k = 0
        if nbdir == 0
            bdirection[nbdir+1] = bindx
            nbdir += 1
        else
            for i = 1:nbdir
                j = bdirection[i]
                tmp = (frac_bv[bindx][1]+frac_bv[j][1])^2 + (frac_bv[bindx][2]+frac_bv[j][2])^2 + (frac_bv[bindx][3]+frac_bv[j][3])^2

                if tmp < 1e-8
                    k = 1
                end
            end

            if k == 0
                bdirection[nbdir+1] = bindx
                nbdir += 1
            end
        end
    end


    return nbdir, bdirection
end


function Set_kplusb(Nkpt, tot_bvector, kg, frac_bv)

    ktp = zeros(Float64, 3)
    b = zeros(Float64, 3)
    kplusb = Vector{Vector{Int32}}(undef, Nkpt)
    for ik = 1:Nkpt
        kplusb[ik] = zeros(Int32, tot_bvector)
    end


    for ik = 1:Nkpt, bindx = 1:tot_bvector
        ktp[1] = kg[ik][1] + frac_bv[bindx][1]
        ktp[2] = kg[ik][2] + frac_bv[bindx][2]
        ktp[3] = kg[ik][3] + frac_bv[bindx][3]

        for i = 1:3
            b[i] = ktp[i]
            if ktp[i] >= 1.0
                b[i] = ktp[i] - 1.0
            end

            if ktp[i] <= 0.0
                b[i] = ktp[i] + 1.0
            end
        end

        kj = 0
        kk = -1
        for ki = 1:Nkpt
            tmpx = abs(b[1]-kg[ki][1])
            tmpy = abs(b[2]-kg[ki][2])
            tmpz = abs(b[3]-kg[ki][3])
            
            if tmpx<1e-6 && tmpy<1e-6 && tmpz<1e-6
                kj = 1
                kk = ki
                break
            end
        end

        if kj == 0
            println("***************************** Error *********************************")
            println("* Check Wannier.Kgrids, the equivalent k points for k+b not found.  *")
            println("***************************** Error *********************************")
            println("***********************************INFO**************************************")
            # println("Reciprocal Lattices lengths are:%10.6f %10.6f %10.6f\n",rtv[0][0],rtv[0][1],rtv[0][2])
            # println("The ratio among them are: b1:b2=%10.6f b1:b3=%10.6f b2:b3=%10.6f\n",rtv[0][0]/rtv[0][1],rtv[0][0]/rtv[0][2],rtv[0][1]/rtv[0][2])
            println("Message: Please try to set Wannier.Kgrid has the similar ratio as above.")
            println("************************************INFO*************************************")
        end

        kplusb[ik][bindx] = kk-1
    end


    return kplusb
end
