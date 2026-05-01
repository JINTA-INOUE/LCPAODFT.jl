@timeit timer "Set_ProExpn_VNA" function Set_ProExpn_VNA!(HVNA, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Set_ProExpn!(HVNA, pao, pspot, system_grid)

    HVNA2 = Set_VNA2(pao, pspot, system_grid)

    Natom = system_grid.Natom
	FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    maxTotal_NumOrbs = maximum(Total_NumOrbs)
    Htemp = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)


    # fill zero HVNA
    hvna_sum = 0
    for atom = 1:Natom

        hvna_count1 = 0
        for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            hvna_count1 += 1
            HVNA[hvna_sum+hvna_count1] = 0.0
        end

        hvna_count = 0
        for Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            hvna_count += 1
        end

        hvna_sum += hvna_count
    end


    # add HVNA2 to HVNA
    hvna2_count = 0
    hvna_sum = 0
    hvna2_sum = 0
    for atom = 1:Natom

        fill!(Htemp, 0.0)
        for Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            hvna2_count += 1
            Htemp[ist,jst] += HVNA2[hvna2_count]
        end

        hvna_count1 = 0
        for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            hvna_count1 += 1
            HVNA[hvna_sum+hvna_count1] = Htemp[ist,jst]
        end

        hvna_count = 0
        for Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            hvna_count += 1
        end

        hvna_sum += hvna_count
        hvna2_sum += hvna_count1
    end
end


function Set_ProExpn!(HVNA, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)


    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    Nspecies = length(pao)
    Total_NumOrbs = system_grid.Total_NumOrbs

    maxL = maximum([pao[spe].Spe_MaxL_Basis for spe = 1:Nspecies]) + BufferL_ProVNA

    Num_RVNA = maxM*(maxL + 1)
    VNA_List = zeros(Int64, Num_RVNA)
    VNA_List2 = zeros(Int64, Num_RVNA)
    lnum = 1
    for l = 0:maxL
        for m = 0:maxM-1
            VNA_List[lnum] = l
            VNA_List2[lnum] = m
            lnum += 1
        end
    end
    Lmax_Four_Int = 2*maxL
    VNATotal_Num = (maxL+1)^2 * maxM



    GL_Abscissae = Gauss_Legendre_x(GL_Mesh)
    
    Sk = Nkmax + Radial_kmin
    Dk = Nkmax - Radial_kmin
    k1 = zeros(Float64, GL_Mesh)
    for ik = 1:GL_Mesh
        k1[ik] = 0.5*(Dk*GL_Abscissae[ik] + Sk)
    end


    VNA_proj_ene = Vector{Vector{Vector{Float64}}}(undef, Nspecies)
    VNA_Bessel = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspecies)
    for spe = 1:Nspecies
        VNA_proj_ene[spe] = Vector{Vector{Float64}}(undef, maxL+1)
        VNA_Bessel[spe] = Vector{Vector{Vector{Float64}}}(undef, maxL+1)
        for L = 0:maxL
            VNA_proj_ene[spe][L+1] = zeros(Float64, maxM)
            VNA_Bessel[spe][L+1] = Vector{Vector{Float64}}(undef, maxM)
            for m = 1:maxM
                VNA_Bessel[spe][L+1][m] = zeros(Float64, GL_Mesh)
            end
        end
    end
    Calc_VNA_Bessel!( pao, pspot, maxL, VNA_proj_ene, VNA_Bessel )

    

    # calc Bessel_Pro00
    Bessel_Pro00 = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_MaxL_Basis = pao[spe].Spe_MaxL_Basis
        Spe_Num_Basis = pao[spe].Spe_Num_Basis
        Bessel_Pro00[spe] = Vector{Vector{Vector{Float64}}}(undef, Spe_MaxL_Basis+1)
        for l = 0:Spe_MaxL_Basis
            Bessel_Pro00[spe][l+1] = Vector{Vector{Float64}}(undef, Spe_Num_Basis[l+1])
            for p = 1:Spe_Num_Basis[l+1]
                Bessel_Pro00[spe][l+1][p] = zeros(Float64, GL_Mesh)
            end
        end
    end
    Calc_Bessel_Pro00!( pao, Bessel_Pro00 )




    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size


    myDS_VNAsize = 0
    for loop = 1:MPI_size, ist = 1:Total_NumOrbs[MPI_atom[loop]], jst = 1:VNATotal_Num
        myDS_VNAsize += 1
    end
    MPI_DS_VNA = zeros(Float64, myDS_VNAsize)



    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    
    SumNL0 = zeros(Float64, Num_RVNA, 4, 4)
    SphB = zeros(Float64, GL_Mesh)
    fsize = maximum(Total_NumOrbs)
    VNAiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!(Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis)
        conj!(Ciα[spe])
    end
    Cjβ = Set_VNAComp2Real(maxL)



    fact2 = zeros(Float64, 2*Lmax_Four_Int+1, 2*Lmax_Four_Int+1)
    for i = 0:2*Lmax_Four_Int, j = 0:2*Lmax_Four_Int
        tmp0 = sqrt(factorial(big(i)))
        tmp1 = sqrt(factorial(big(j)))
        fact2[i+1,j+1] = tmp0/tmp1
    end



    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    SphB = Vector{Vector{Float64}}(undef, Lmax_Four_Int+1)
    for l = 1:Lmax_Four_Int+1
        SphB[l] = zeros(Float64, GL_Mesh)
    end
    tmpH1 = Vector{Vector{ComplexF64}}(undef, maxL+1)
    for l = 0:maxL
        tmpH1[l+1] = zeros(ComplexF64, 2*l+1)
    end
    tmpH = zeros(ComplexF64, fsize)
    tmpL = zeros(Float64, GL_Mesh)

    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)
    

    SH = zeros(Float64, 2)



    MPI.Barrier(comm)
    hst = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        NO0 = Total_NumOrbs[atom]

        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        cell = MPI_ncn[loop]+1
        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]
        R = sqrt(x^2 + y^2 + z^2)
        R = ifelse(R < 1e-10, 1e-10, R)

    

        for ik = 1:GL_Mesh
            Calc_SphericalBesselj2!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
            end
        end


        fill!(VNAiαjβ, 0.0)

        for L = 0:Lmax_Four_Int
                
            for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:Num_RVNA
                ll = VNA_List[lnum]
                index = VNA_List2[lnum]+1
                @. tmpL = SphB[L+1]*Bessel_Pro00[ispe][l+1][p]
                SumNL0[lnum,p,l+1] = dot(tmpL, VNA_Bessel[jspe][ll+1][index])
            end

            for M = -L:L
                ist = 0
                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                    jst = 0
                    ist += 1
                    for lnum = 1:Num_RVNA, mm = -VNA_List[lnum]:VNA_List[lnum]
                        ll = VNA_List[lnum]
                        jst += 1
                        if abs(ll-L) <= l <= abs(ll+L) && iszero(m-mm-M) && abs(m-M) <= ll
                            indx0 = L-abs(M)+1
                            indx1 = L+abs(M)+1
                            Ylm_complex!(L,M,fact2[indx0,indx1],x,y,z,SH)
                            Ylm = ComplexF64(SH[1], SH[2])
                            Ls = Float64(L+ll-l)
                            iYC = (-im)^Ls * conj(Ylm) * Gaunt(f,l,m,ll,mm,L,M)

                            VNAiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                        end
                    end
                end
            end
        end


        # complex to real
        for ist = 1:NO0
            tot = 1
            for lnum = 1:Num_RVNA
                ll = VNA_List[lnum]
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAiαjβ[ist, tot:tot+2*ll])
                @views VNAiαjβ[ist, tot:tot+2*ll] = tmpH1[ll+1]
                tot += 2*ll+1
            end
        end

        for jst = 1:VNATotal_Num
            @views mul!(tmpH, Ciα[ispe], VNAiαjβ[:,jst])
            @views VNAiαjβ[:,jst] = tmpH
        end
            
        
        for ist = 1:NO0, jst = 1:VNATotal_Num
            hst += 1
            MPI_DS_VNA[hst] = 8*real(VNAiαjβ[ist,jst])
        end
    end
    MPI.Barrier(comm)




    Total_DS_VNAsize = MPI.Allreduce(myDS_VNAsize, MPI.SUM, comm) 

    _counts = zeros(Int64, nprocs)
    for id = 1:nprocs
        if id-1 == myrank
            _counts[id] = myDS_VNAsize
        end
        MPI.Barrier(comm)
    end
    MPI.Barrier(comm)

    watemp = similar(_counts)
    MPI.Allreduce!(_counts,watemp,nprocs,MPI.SUM,comm) 


    if myrank == 0
        All_DS_VNA = zeros(Float64, Total_DS_VNAsize)
        MPI.Gatherv!(MPI_DS_VNA,VBuffer(All_DS_VNA, watemp),comm)
    else
        All_DS_VNA = Vector{Float64}(undef, Total_DS_VNAsize)
        MPI.Gatherv!(MPI_DS_VNA,nothing,0,comm)
    end
    MPI.Barrier(comm)
    MPI.Bcast!(All_DS_VNA, 0, comm)


    
    FNAN = system_grid.FNAN
	natn = system_grid.natn
    
    DS_VNA = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
	for atom = 1:Natom
		DS_VNA[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
		for Rn = 1:FNAN[atom]+1
			DS_VNA[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
			for ist = 1:Total_NumOrbs[atom]
				DS_VNA[atom][Rn][ist] = Vector{Float64}(undef, VNATotal_Num)
			end
		end
	end

    count = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:VNATotal_Num
        count += 1
        DS_VNA[atom][Rn][ist][jst] = All_DS_VNA[count]
    end



    VNAE = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        VNAE[atom] = zeros(Float64, VNATotal_Num)

        count = 0
        for lnum = 1:Num_RVNA
            L = VNA_List[lnum]
            m = VNA_List2[lnum]
            ene = VNA_proj_ene[spe][L+1][m+1]
            for m = 1:2*L+1
                count += 1
                VNAE[atom][count] = ene
            end
        end
    end



    MPI_Hsize = system_grid.MPI_Hsize
    MPI_HVNA = zeros(Float64, MPI_Hsize[myrank+1])


    Atom_Cut1 = system_grid.Atom_Cut1
    MPI_Dis = system_grid.MPI_Dis
    MPI_RMI = system_grid.MPI_RMI
    HVNA_temp = zeros(Float64, fsize, fsize)
    tmp = zeros(Float64, VNATotal_Num)
    count = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        fill!(HVNA_temp, 0.0)

        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kl = MPI_RMI[loop][Rm]
            if kl >= 0
                for ist = 1:NO0, jst = 1:NO1
                    @. tmp = DS_VNA[jatom][kl+1][jst]*VNAE[kg]
                    HVNA_temp[ist,jst] += dot(DS_VNA[atom][Rm][ist], tmp)
                end
            end
        end


        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, MPI_Dis[loop])
        for ist = 1:NO0, jst = 1:NO1
            count += 1
            MPI_HVNA[count] = dmp * HVNA_temp[ist,jst]
        end
    end
    MPI.Barrier(comm)



    MPI.Allgatherv!(MPI_HVNA,VBuffer(HVNA,MPI_Hsize),comm)
end


function Set_VNA2(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    Total_NumOrbs = system_grid.Total_NumOrbs
    HVNA2size = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
        HVNA2size += 1
    end
    HVNA2 = zeros(Float64, HVNA2size)

    Set_VNA2!(HVNA2, pao, pspot, system_grid)

    return HVNA2
end


function Set_VNA2!(HVNA2, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nspecies = length(pao)
    atom2spe = system_grid.atom2spe
    Total_NumOrbs = system_grid.Total_NumOrbs
    

    Lmax_Four_Int = 0
    for spe = 1:Nspecies
        Spe_MaxL_Basis = pao[spe].Spe_MaxL_Basis
        Lmax_Four_Int = max(Spe_MaxL_Basis, Lmax_Four_Int)
    end

    

    GL_Abscissae, GL_Weight = Gauss_Legendre(GL_Mesh)    
    Sk = Nkmax + Radial_kmin
    Dk = Nkmax - Radial_kmin
    k1 = zeros(Float64, GL_Mesh)
    for ik = 1:GL_Mesh
        k1[ik] = 0.5*(Dk*GL_Abscissae[ik] + Sk)
    end

    

    Spe_CrudeVNA_Bessel = Vector{Vector{Float64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_CrudeVNA_Bessel[spe] = zeros(Float64, GL_Mesh)
    end
    FT_VNA!(pao, pspot, Spe_CrudeVNA_Bessel)

    for spe = 1:Nspecies
        @. Spe_CrudeVNA_Bessel[spe] = 0.5*Dk*GL_Weight*k1*k1*Spe_CrudeVNA_Bessel[spe]
    end

    
    
    Spe_ProductRF_Bessel = Vector{Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}}(undef, Nspecies)
    for spe = 1:Nspecies

        MaxL = pao[spe].Spe_MaxL_Basis
        Mu = pao[spe].Spe_Num_Basis

        Spe_ProductRF_Bessel[spe] = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, MaxL+1)
        for GL1 = 0:MaxL
            Spe_ProductRF_Bessel[spe][GL1+1] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Mu[GL1+1])
            for Mul1 = 1:Mu[GL1+1]
                Spe_ProductRF_Bessel[spe][GL1+1][Mul1] = Vector{Vector{Vector{Vector{Float64}}}}(undef, MaxL+1)
                for GL2 = 0:MaxL
                    Spe_ProductRF_Bessel[spe][GL1+1][Mul1][GL2+1] = Vector{Vector{Vector{Float64}}}(undef, Mu[GL2+1])
                    
                    Lmax = ifelse(GL1 <= GL2, 2*GL2, 1)
                    num = ifelse(GL1 <= GL2, GL_Mesh, 1)
                    
                    for Mul2 = 1:Mu[GL2+1]
                        Spe_ProductRF_Bessel[spe][GL1+1][Mul1][GL2+1][Mul2] = Vector{Vector{Float64}}(undef, Lmax+1)
                        for l = 0:Lmax
                            Spe_ProductRF_Bessel[spe][GL1+1][Mul1][GL2+1][Mul2][l+1] = zeros(Float64, num)
                        end
                    end
                end
            end
        end
    end

    for spe = 1:Nspecies
        FT_ProductPAO!(pao[spe], pspot[spe].Spe_VPS_RV[begin], Spe_ProductRF_Bessel[spe])
    end
    


    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size


    myHVNA2size = 0
    for loop = 1:MPI_size, _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:Total_NumOrbs[MPI_atom[loop]]
        myHVNA2size += 1
    end
    MPI_HVNA2 = zeros(Float64, myHVNA2size)


    
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
     
    # error check
    Max_Spe_MaxL_Basis = 3
    for spe = 1:Nspecies
        L = pao[spe].Spe_MaxL_Basis
        if L > Max_Spe_MaxL_Basis
            println("Spe_MaxL_Basis = $L")
            error("please check")
        end
    end

    SumHVNA2 = zeros(Float64, 9, 4, 4, 4, 4)


    fsize = maximum(Total_NumOrbs)
    VNA2iαjβ = zeros(ComplexF64, fsize, fsize)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    Cjβ = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!( Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis )
        Cjβ[spe] = deepcopy(Ciα[spe])
        conj!(Ciα[spe])
    end

    
    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    
    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    SphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    for l = 1:2*Lmax_Four_Int+3
        SphB[l] = zeros(Float64, GL_Mesh)
    end
    tempL = zeros(Float64, GL_Mesh)
    tmpH = zeros(ComplexF64, fsize)



    MPI.Barrier(comm)
    counts = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        MaxL_Basis = pao[ispe].Spe_MaxL_Basis
        Num_Basis = pao[ispe].Spe_Num_Basis

        jatom = MPI_natn[loop]
        cell = MPI_ncn[loop]+1
        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]
        jspe = atom2spe[jatom]

        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[atom]

        R = sqrt(x^2 + y^2 + z^2)
        R = ifelse(R < 1e-10, 1e-10, R)
        
        Lmax_Four_Int = 2*MaxL_Basis
        for ik = 1:GL_Mesh
            Calc_SphericalBesselj2!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
            end
        end


        for l = 0:MaxL_Basis, p = 1:Num_Basis[l+1], ll = 0:MaxL_Basis, pp = 1:Num_Basis[ll+1]
            if l <= ll
                Lmax_Four_Int = 2*ll
                for L = 0:Lmax_Four_Int
                    if abs(ll-L) <= l <= ll+L
                        @. tempL = SphB[L+1]*Spe_CrudeVNA_Bessel[jspe]
                        @views SumHVNA2[L+1,pp,ll+1,p,l+1] = dot(tempL, Spe_ProductRF_Bessel[ispe][l+1][p][ll+1][pp][L+1])
                    end
                end
            end
        end

        fill!(VNA2iαjβ, 0.0)

        ist = 0
        for l = 0:MaxL_Basis, p = 1:Num_Basis[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis, pp = 1:Num_Basis[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    Lmax_Four_Int = 2*ll
                    for L = 0:Lmax_Four_Int, M = -L:L
                        if abs(ll-L) <= l <= ll+L && iszero(m-mm+M)
                            gaunt = (-1.0)^abs(M)*Gaunt(f,l,m,ll,mm,L,-M)
                            Ylm = Ylm_complex(L,M,x,y,z)
                            YC = Ylm * gaunt
                            
                            VNA2iαjβ[ist,jst] += YC*SumHVNA2[L+1,pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end


        ist = 0
        for l = 0:MaxL_Basis, p = 1:Num_Basis[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis, pp = 1:Num_Basis[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    VNA2iαjβ[jst,ist] = conj(VNA2iαjβ[ist,jst])
                end
            end
        end
            
                
        # complex to real
        for ist = 1:NO0
            @views mul!(tmpH, Cjβ[ispe], VNA2iαjβ[ist,:])
            @views VNA2iαjβ[ist,:] = tmpH
        end

        for jst = 1:NO1
            @views mul!(tmpH, Ciα[ispe], VNA2iαjβ[:,jst])
            @views VNA2iαjβ[:,jst] = tmpH
        end

            
        for ist = 1:NO0, jst = 1:NO1
            counts += 1
            MPI_HVNA2[counts] = 8*real(VNA2iαjβ[ist,jst])
        end
    end
    MPI.Barrier(comm)





    _counts = zeros(Int64, nprocs)
    for id = 1:nprocs
        if id-1 == myrank
            _counts[id] = myHVNA2size
        end
        MPI.Barrier(comm)
    end
    MPI.Barrier(comm)

    MPI_HVNA2size = similar(_counts)
    MPI.Allreduce!(_counts, MPI_HVNA2size, nprocs, MPI.SUM, comm) 


    MPI.Allgatherv!(MPI_HVNA2, VBuffer(HVNA2, MPI_HVNA2size), comm)
end
