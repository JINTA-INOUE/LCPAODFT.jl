function Set_DS_VNAforce(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Natom = system_grid.Natom
    Nspecies = length(pao)
    FNAN = system_grid.FNAN
    Total_NumOrbs = system_grid.Total_NumOrbs

    maxL = maximum([pao[spe].Spe_MaxL_Basis for spe = 1:Nspecies]) + BufferL_ProVNA
    VNATotal_Num = (maxL+1)^2 * maxM


    maxFNAN = maximum(FNAN)
	maxTotal_NumOrbs = maximum(Total_NumOrbs)
	DS_VNAforce = Vector{Vector{Vector{Vector{Vector{Float32}}}}}(undef, 4)
	for xyz = 1:4
		DS_VNAforce[xyz] = Vector{Vector{Vector{Vector{Float32}}}}(undef, Natom+1)
		for atom = 1:Natom+1
			if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
			DS_VNAforce[xyz][atom] = Vector{Vector{Vector{Float32}}}(undef, fan)
			for Rn = 1:fan
				DS_VNAforce[xyz][atom][Rn] = Vector{Vector{Float32}}(undef, NO0)
				for ist = 1:NO0
					DS_VNAforce[xyz][atom][Rn][ist] = zeros(Float32, VNATotal_Num)
				end
			end
		end
	end

    Set_DS_VNAforce!(DS_VNAforce, pao, pspot, system_grid)


    return DS_VNAforce
end



function Set_DS_VNAforce!(DS_VNAforce, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
    
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
            VNA_proj_ene[spe][L+1] = Vector{Float64}(undef, maxM)
            VNA_Bessel[spe][L+1] = Vector{Vector{Float64}}(undef, maxM)
            for m = 1:maxM
                VNA_Bessel[spe][L+1][m] = zeros(Float64, GL_Mesh)
            end
        end
    end
    Calc_VNA_Bessel!(pao, pspot, maxL, VNA_proj_ene, VNA_Bessel)

    

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
    Calc_Bessel_Pro00!(pao, Bessel_Pro00)



    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size
    FNAN = system_grid.FNAN
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz


    myDS_VNAsize = 0
    for loop = 1:MPI_size, ist = 1:Total_NumOrbs[MPI_atom[loop]], jst = 1:VNATotal_Num
        myDS_VNAsize += 1
    end
    MPI_DS_VNA1D = zeros(Float32, myDS_VNAsize)
    MPI_DS_VNA1Dx = zeros(Float32, myDS_VNAsize)
    MPI_DS_VNA1Dy = zeros(Float32, myDS_VNAsize)
    MPI_DS_VNA1Dz = zeros(Float32, myDS_VNAsize)

    


    SumNL0 = zeros(Float64, Num_RVNA, 4, 4)
    SumNLr0 = zeros(Float64, Num_RVNA, 4, 4)
    fsize = maximum(Total_NumOrbs)
    VNAiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    VNAriαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    VNAtiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    VNApiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!(Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis)
        conj!(Ciα[spe])
    end
    Cjβ = Set_VNAComp2Real( maxL )


    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    fact2 = zeros(Float64, 2*Lmax_Four_Int+1, 2*Lmax_Four_Int+1)
    for i = 0:2*Lmax_Four_Int, j = 0:2*Lmax_Four_Int
        tmp0 = sqrt(factorial(big(i)))
        tmp1 = sqrt(factorial(big(j)))
        fact2[i+1,j+1] = tmp0/tmp1
    end


    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 30)
    dSphB_l = zeros(Float64, 30)
    SphB = Vector{Vector{Float64}}(undef, Lmax_Four_Int+1)
    dSphB = Vector{Vector{Float64}}(undef, Lmax_Four_Int+1)
    for l = 1:Lmax_Four_Int+1
        SphB[l] = zeros(Float64, GL_Mesh)
        dSphB[l] = zeros(Float64, GL_Mesh)
    end
    tmpH1 = Vector{Vector{ComplexF64}}(undef, maxL+1)
    for l = 0:maxL
        tmpH1[l+1] = zeros(ComplexF64, 2*l+1)
    end
    tmpH = zeros(ComplexF64, fsize)
    tmpL = zeros(Float64, GL_Mesh)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)




    MPI.Barrier(comm)
    counts = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        NO0 = Total_NumOrbs[atom]

        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        cell = MPI_ncn[loop]+1

        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]
        

        R, theta, phi = xyz_to_spherical(x, y, z)
        R = ifelse(R < 1.0e-10, 1.0e-10, R)
        siT = sin(theta)
        coT = cos(theta)
        siP = sin(phi)
        coP = cos(phi)
            

        for ik = 1:GL_Mesh
            Calc_SphericalBesselj2!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                dSphB[l][ik] = dSphB_l[l]*k1[ik]
            end
        end


        fill!(VNAiαjβ, 0.0)
        fill!(VNAriαjβ, 0.0)
        fill!(VNAtiαjβ, 0.0)
        fill!(VNApiαjβ, 0.0)


        for L = 0:Lmax_Four_Int

            for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:Num_RVNA
                ll = VNA_List[lnum]
                index = VNA_List2[lnum]+1
                @. tmpL = Bessel_Pro00[ispe][l+1][p]*VNA_Bessel[jspe][ll+1][index]
                SumNL0[lnum,p,l+1] = dot(SphB[L+1], tmpL)
                SumNLr0[lnum,p,l+1] = dot(dSphB[L+1], tmpL)
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
                            Ylm_complex!(L,M,fact2[indx0,indx1],theta,phi,SH,dSHt,dSHp)
                            Ls = Float64(L+ll-l)
                            gaunt = Gaunt(f,l,m,ll,mm,L,M)
                            tmp = (-im)^Ls
                                
                            Ylm = ComplexF64(SH[1], SH[2])
                            dYlmdtheta = ComplexF64(dSHt[1], dSHt[2])
                            dYlmdphi = ComplexF64(dSHp[1], dSHp[2])
                            
                            iYC = conj(Ylm) * tmp * gaunt
                            iYCt = conj(dYlmdtheta) * tmp * gaunt
                            iYCp = conj(dYlmdphi) * tmp * gaunt

                            VNAiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                            VNAriαjβ[ist,jst] += iYC*SumNLr0[lnum,p,l+1]
                            VNAtiαjβ[ist,jst] += iYCt*SumNL0[lnum,p,l+1]
                            VNApiαjβ[ist,jst] += iYCp*SumNL0[lnum,p,l+1]
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
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAiαjβ[ist,tot:tot+2*ll])
                @views VNAiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAriαjβ[ist,tot:tot+2*ll])
                @views VNAriαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAtiαjβ[ist,tot:tot+2*ll])
                @views VNAtiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNApiαjβ[ist,tot:tot+2*ll])
                @views VNApiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                tot += 2*ll+1
            end
        end


        for jst = 1:VNATotal_Num
            @views mul!(tmpH, Ciα[ispe], VNAiαjβ[:,jst])
            @views VNAiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], VNAriαjβ[:,jst])
            @views VNAriαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], VNAtiαjβ[:,jst])
            @views VNAtiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], VNApiαjβ[:,jst])
            @views VNApiαjβ[:,jst] = tmpH
        end


        if Rn ≠ 1
            if abs(siT) < 1.0e-13
                for ist = 1:NO0, jst = 1:VNATotal_Num
                    counts += 1
                    MPI_DS_VNA1D[counts] = 8*real(VNAiαjβ[ist,jst])
                    MPI_DS_VNA1Dx[counts] = -8*real(siT*coP*VNAriαjβ[ist,jst] + coT*coP/R*VNAtiαjβ[ist,jst])
                    MPI_DS_VNA1Dy[counts] = -8*real(siT*siP*VNAriαjβ[ist,jst] + coT*siP/R*VNAtiαjβ[ist,jst])
                    MPI_DS_VNA1Dz[counts] = -8*real(coT*VNAriαjβ[ist,jst] - siT/R*VNAtiαjβ[ist,jst])
                end
            else
                for ist = 1:NO0, jst = 1:VNATotal_Num
                    counts += 1
                    MPI_DS_VNA1D[counts] = 8*real(VNAiαjβ[ist,jst])
                    MPI_DS_VNA1Dx[counts] = -8*real(siT*coP*VNAriαjβ[ist,jst] + coT*coP/R*VNAtiαjβ[ist,jst] - siP/siT/R*VNApiαjβ[ist,jst])
                    MPI_DS_VNA1Dy[counts] = -8*real(siT*siP*VNAriαjβ[ist,jst] + coT*siP/R*VNAtiαjβ[ist,jst] + coP/siT/R*VNApiαjβ[ist,jst])
                    MPI_DS_VNA1Dz[counts] = -8*real(coT*VNAriαjβ[ist,jst] - siT/R*VNAtiαjβ[ist,jst])
                end
            end
        else
            for ist = 1:NO0, jst = 1:VNATotal_Num
                counts += 1
                MPI_DS_VNA1D[counts] = 8*real(VNAiαjβ[ist,jst])
                MPI_DS_VNA1Dx[counts] = 0.0
                MPI_DS_VNA1Dy[counts] = 0.0
                MPI_DS_VNA1Dz[counts] = 0.0
            end
        end
    end



    VNAE = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        VNAE[atom] = zeros(Float64, VNATotal_Num)

        counts = 0
        for lnum = 1:Num_RVNA
            L = VNA_List[lnum]
            m = VNA_List2[lnum]
            ene = VNA_proj_ene[spe][L+1][m+1]
            for m = 1:2*L+1
                counts += 1
                VNAE[atom][counts] = ene
            end
        end
    end

    
    counts = 0
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        for ist = 1:Total_NumOrbs[atom], jst = 1:VNATotal_Num
            ene = VNAE[jatom][jst]
            counts += 1
            MPI_DS_VNA1Dx[counts] = ene*MPI_DS_VNA1Dx[counts]
            MPI_DS_VNA1Dy[counts] = ene*MPI_DS_VNA1Dy[counts]
            MPI_DS_VNA1Dz[counts] = ene*MPI_DS_VNA1Dz[counts]
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

    MPI_DS_VNAsize = zeros(Int64, nprocs)
    MPI.Allreduce!(_counts, MPI_DS_VNAsize, nprocs, MPI.SUM, comm)  


    
    MPI_DS_VNAxyz = zeros(Float32, Total_DS_VNAsize)
    for i = 1:4
        if i == 1
            MPI.Allgatherv!(MPI_DS_VNA1D, VBuffer(MPI_DS_VNAxyz, MPI_DS_VNAsize), comm)
        elseif i == 2
            MPI.Allgatherv!(MPI_DS_VNA1Dx, VBuffer(MPI_DS_VNAxyz, MPI_DS_VNAsize), comm)
        elseif i == 3
            MPI.Allgatherv!(MPI_DS_VNA1Dy, VBuffer(MPI_DS_VNAxyz, MPI_DS_VNAsize), comm)
        elseif i == 4
            MPI.Allgatherv!(MPI_DS_VNA1Dz, VBuffer(MPI_DS_VNAxyz, MPI_DS_VNAsize), comm)
        end

        counts = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:VNATotal_Num
            counts += 1
            DS_VNAforce[i][atom][Rn][ist][jst] = MPI_DS_VNAxyz[counts]
        end
    end
end



function Set_HVNA2_3force(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Natom = system_grid.Natom
	FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HVNA2force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
    HVNA3force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
	for xyz = 1:3
		HVNA2force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		HVNA3force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			HVNA2force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			HVNA3force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				HVNA2force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					HVNA2force[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[atom])
				end

                HVNA3force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[natn[atom][Rn]])
                for ist = 1:Total_NumOrbs[natn[atom][Rn]]
					HVNA3force[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    Set_HVNA2_3force!(HVNA2force, HVNA3force, pao, pspot, system_grid)


    return HVNA2force, HVNA3force
end


function Set_HVNA2_3force!(HVNA2force, HVNA3force, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
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
    FT_VNA!( pao, pspot, Spe_CrudeVNA_Bessel )

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
        FT_ProductPAO!( pao[spe], pspot[spe].Spe_VPS_RV[begin], Spe_ProductRF_Bessel[spe] )
    end


    FNAN = system_grid.FNAN
    natn = system_grid.natn
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size


    myHVNA2size = 0
    for loop = 1:MPI_size, ist = 1:Total_NumOrbs[MPI_atom[loop]], jst = 1:Total_NumOrbs[MPI_atom[loop]]
        myHVNA2size += 1
    end

    myHVNA3size = 0
    for loop = 1:MPI_size, ist = 1:Total_NumOrbs[MPI_natn[loop]], jst = 1:Total_NumOrbs[MPI_natn[loop]]
        myHVNA3size += 1
    end

    MPI_HVNA21Dx = zeros(Float64, myHVNA2size)
    MPI_HVNA21Dy = zeros(Float64, myHVNA2size)
    MPI_HVNA21Dz = zeros(Float64, myHVNA2size)
    MPI_HVNA31Dx = zeros(Float64, myHVNA3size)
    MPI_HVNA31Dy = zeros(Float64, myHVNA3size)
    MPI_HVNA31Dz = zeros(Float64, myHVNA3size)
     


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
    SumHVNAr2 = zeros(Float64, 9, 4, 4, 4, 4)
    SumHVNA3 = zeros(Float64, 9, 4, 4, 4, 4)
    SumHVNAr3 = zeros(Float64, 9, 4, 4, 4, 4)


    fsize = maximum(Total_NumOrbs)
    VNA2riαjβ = zeros(ComplexF64, fsize, fsize)
    VNA2tiαjβ = zeros(ComplexF64, fsize, fsize)
    VNA2piαjβ = zeros(ComplexF64, fsize, fsize)
    VNA3riαjβ = zeros(ComplexF64, fsize, fsize)
    VNA3tiαjβ = zeros(ComplexF64, fsize, fsize)
    VNA3piαjβ = zeros(ComplexF64, fsize, fsize)
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

    fact2 = zeros(Float64, (2*Lmax_Four_Int+1)^2+1, (2*Lmax_Four_Int+1)^2+1)
    for i = 0:(2*Lmax_Four_Int+1)^2, j = 0:(2*Lmax_Four_Int+1)^2
        tmp0 = sqrt(factorial(big(i)))
        tmp1 = sqrt(factorial(big(j)))
        fact2[i+1,j+1] = tmp0/tmp1
    end


    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    dSphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    SphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    dSphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    for l = 1:2*Lmax_Four_Int+3
        SphB[l] = zeros(Float64, GL_Mesh)
        dSphB[l] = zeros(Float64, GL_Mesh)
    end
    tempL = zeros(Float64, GL_Mesh)
    tmpH = zeros(ComplexF64, fsize)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)


    MPI.Barrier(comm)
    counts1 = 0
    counts2 = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        NO2 = Total_NumOrbs[atom]
        MaxL_Basis2 = pao[ispe].Spe_MaxL_Basis
        Num_Basis2 = pao[ispe].Spe_Num_Basis

        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        cell = MPI_ncn[loop]+1
        NO3 = Total_NumOrbs[jatom]
        MaxL_Basis3 = pao[jspe].Spe_MaxL_Basis
        Num_Basis3 = pao[jspe].Spe_Num_Basis
        
        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]

        R2, theta2, phi2 = xyz_to_spherical(x, y, z)
        R2 = ifelse(R2 < 1.0e-10, 1.0e-10, R2)
        siT2 = sin(theta2)
        coT2 = cos(theta2)
        siP2 = sin(phi2)
        coP2 = cos(phi2)


        R3, theta3, phi3 = xyz_to_spherical(-x, -y, -z)
        R3 = ifelse(R3 < 1.0e-10, 1.0e-10, R3)
        siT3 = sin(theta3)
        coT3 = cos(theta3)
        siP3 = sin(phi3)
        coP3 = cos(phi3)
        R = R2
        

        Lmax_Four_Int = 2*max(MaxL_Basis2, MaxL_Basis3)
        for ik = 1:GL_Mesh
            Calc_SphericalBesselj2!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                dSphB[l][ik] = dSphB_l[l]*k1[ik]
            end
        end

            
        for l = 0:MaxL_Basis2, p = 1:Num_Basis2[l+1], ll = 0:MaxL_Basis2, pp = 1:Num_Basis2[ll+1]
            if l <= ll
                Lmax_Four_Int = 2*ll
                for L = 0:Lmax_Four_Int
                    if abs(ll-L) <= l <= ll+L
                        @. tempL = Spe_ProductRF_Bessel[ispe][l+1][p][ll+1][pp][L+1]*Spe_CrudeVNA_Bessel[jspe]
                        @views SumHVNA2[L+1,pp,ll+1,p,l+1] = dot(SphB[L+1], tempL)
                        @views SumHVNAr2[L+1,pp,ll+1,p,l+1] = dot(dSphB[L+1], tempL)
                    end
                end
            end
        end


        for l = 0:MaxL_Basis3, p = 1:Num_Basis3[l+1], ll = 0:MaxL_Basis3, pp = 1:Num_Basis3[ll+1]
            if l <= ll
                Lmax_Four_Int = 2*ll
                for L = 0:Lmax_Four_Int
                    if abs(ll-L) <= l <= ll+L
                        @. tempL = Spe_ProductRF_Bessel[jspe][l+1][p][ll+1][pp][L+1]*Spe_CrudeVNA_Bessel[ispe]
                        @views SumHVNA3[L+1,pp,ll+1,p,l+1] = dot(SphB[L+1], tempL)
                        @views SumHVNAr3[L+1,pp,ll+1,p,l+1] = dot(dSphB[L+1], tempL)
                    end
                end
            end
        end

        fill!(VNA2riαjβ, 0.0)
        fill!(VNA2tiαjβ, 0.0)
        fill!(VNA2piαjβ, 0.0)
        fill!(VNA3riαjβ, 0.0)
        fill!(VNA3tiαjβ, 0.0)
        fill!(VNA3piαjβ, 0.0)


        ist = 0
        for l = 0:MaxL_Basis2, p = 1:Num_Basis2[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis2, pp = 1:Num_Basis2[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    Lmax_Four_Int = 2*ll
                    for L = 0:Lmax_Four_Int, M = -L:L
                        if abs(ll-L) <= l <= ll+L && iszero(m-mm+M)
                                
                            gaunt = (-1.0)^abs(M)*Gaunt(f,l,m,ll,mm,L,-M)

                            indx0 = L-abs(M)+1
                            indx1 = L+abs(M)+1
                            Ylm_complex!(L,M,fact2[indx0,indx1],theta2,phi2,SH,dSHt,dSHp)
                            Ylm = ComplexF64(SH[1], SH[2])
                            dYlmdtheta = ComplexF64(dSHt[1], dSHt[2])
                            dYlmdphi = ComplexF64(dSHp[1], dSHp[2])
                                
                            YC = Ylm * gaunt
                            YCt = dYlmdtheta * gaunt
                            YCp = dYlmdphi * gaunt

                            VNA2riαjβ[ist,jst] += YC*SumHVNAr2[L+1,pp,ll+1,p,l+1]
                            VNA2tiαjβ[ist,jst] += YCt*SumHVNA2[L+1,pp,ll+1,p,l+1]
                            VNA2piαjβ[ist,jst] += YCp*SumHVNA2[L+1,pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end



        ist = 0
        for l = 0:MaxL_Basis3, p = 1:Num_Basis3[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis3, pp = 1:Num_Basis3[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    Lmax_Four_Int = 2*ll
                    for L = 0:Lmax_Four_Int, M = -L:L
                        if abs(ll-L) <= l <= ll+L && iszero(m-mm+M)
                                
                            gaunt = (-1.0)^abs(M)*Gaunt(f,l,m,ll,mm,L,-M)

                            indx0 = L-abs(M)+1
                            indx1 = L+abs(M)+1
                            Ylm_complex!(L,M,fact2[indx0,indx1],theta3,phi3,SH,dSHt,dSHp)
                                
                            Ylm = ComplexF64(SH[1], SH[2])
                            dYlmdtheta = ComplexF64(dSHt[1], dSHt[2])
                            dYlmdphi = ComplexF64(dSHp[1], dSHp[2])
                                
                            YC = Ylm * gaunt
                            YCt = dYlmdtheta * gaunt
                            YCp = dYlmdphi * gaunt

                            VNA3riαjβ[ist,jst] += YC*SumHVNAr3[L+1,pp,ll+1,p,l+1]
                            VNA3tiαjβ[ist,jst] += YCt*SumHVNA3[L+1,pp,ll+1,p,l+1]
                            VNA3piαjβ[ist,jst] += YCp*SumHVNA3[L+1,pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end


        ist = 0
        for l = 0:MaxL_Basis2, p = 1:Num_Basis2[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis2, pp = 1:Num_Basis2[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    VNA2riαjβ[jst,ist] = conj(VNA2riαjβ[ist,jst])
                    VNA2tiαjβ[jst,ist] = conj(VNA2tiαjβ[ist,jst])
                    VNA2piαjβ[jst,ist] = conj(VNA2piαjβ[ist,jst])
                end
            end
        end


        ist = 0
        for l = 0:MaxL_Basis3, p = 1:Num_Basis3[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis3, pp = 1:Num_Basis3[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    VNA3riαjβ[jst,ist] = conj(VNA3riαjβ[ist,jst])
                    VNA3tiαjβ[jst,ist] = conj(VNA3tiαjβ[ist,jst])
                    VNA3piαjβ[jst,ist] = conj(VNA3piαjβ[ist,jst])
                end
            end
        end


        # complex to real
        for ist = 1:NO2
            @views mul!(tmpH, Cjβ[ispe], VNA2riαjβ[ist,:])
            @views VNA2riαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[ispe], VNA2tiαjβ[ist,:])
            @views VNA2tiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[ispe], VNA2piαjβ[ist,:])
            @views VNA2piαjβ[ist,:] = tmpH
        end

        for jst = 1:NO2
            @views mul!(tmpH, Ciα[ispe], VNA2riαjβ[:,jst])
            @views VNA2riαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], VNA2tiαjβ[:,jst])
            @views VNA2tiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], VNA2piαjβ[:,jst])
            @views VNA2piαjβ[:,jst] = tmpH
        end

        for ist = 1:NO3
            @views mul!(tmpH, Cjβ[jspe], VNA3riαjβ[ist,:])
            @views VNA3riαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], VNA3tiαjβ[ist,:])
            @views VNA3tiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], VNA3piαjβ[ist,:])
            @views VNA3piαjβ[ist,:] = tmpH
        end
    
        for jst = 1:NO3
            @views mul!(tmpH, Ciα[jspe], VNA3riαjβ[:,jst])
            @views VNA3riαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[jspe], VNA3tiαjβ[:,jst])
            @views VNA3tiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[jspe], VNA3piαjβ[:,jst])
            @views VNA3piαjβ[:,jst] = tmpH
        end

            
        
            
        if Rn ≠ 1
            if abs(siT2) < 1.0e-13
                for ist = 1:NO2, jst = 1:NO2
                    counts1 += 1
                    MPI_HVNA21Dx[counts1] = -8*real(siT2*coP2*VNA2riαjβ[ist,jst] + coT2*coP2/R*VNA2tiαjβ[ist,jst])
                    MPI_HVNA21Dy[counts1] = -8*real(siT2*siP2*VNA2riαjβ[ist,jst] + coT2*siP2/R*VNA2tiαjβ[ist,jst])
                    MPI_HVNA21Dz[counts1] = -8*real(coT2*VNA2riαjβ[ist,jst] - siT2/R*VNA2tiαjβ[ist,jst])
                end
            else
                for ist = 1:NO2, jst = 1:NO2
                    counts1 += 1
                    MPI_HVNA21Dx[counts1] = -8*real(siT2*coP2*VNA2riαjβ[ist,jst] + coT2*coP2/R*VNA2tiαjβ[ist,jst] - siP2/siT2/R*VNA2piαjβ[ist,jst])
                    MPI_HVNA21Dy[counts1] = -8*real(siT2*siP2*VNA2riαjβ[ist,jst] + coT2*siP2/R*VNA2tiαjβ[ist,jst] + coP2/siT2/R*VNA2piαjβ[ist,jst])
                    MPI_HVNA21Dz[counts1] = -8*real(coT2*VNA2riαjβ[ist,jst] - siT2/R*VNA2tiαjβ[ist,jst])
                end
            end

            if abs(siT3) < 1.0e-13
                for ist = 1:NO3, jst = 1:NO3
                    counts2 += 1
                    MPI_HVNA31Dx[counts2] = -8*real(siT3*coP3*VNA3riαjβ[ist,jst] + coT3*coP3/R*VNA3tiαjβ[ist,jst])
                    MPI_HVNA31Dy[counts2] = -8*real(siT3*siP3*VNA3riαjβ[ist,jst] + coT3*siP3/R*VNA3tiαjβ[ist,jst])
                    MPI_HVNA31Dz[counts2] = -8*real(coT3*VNA3riαjβ[ist,jst] - siT3/R*VNA3tiαjβ[ist,jst])
                end
            else
                for ist = 1:NO3, jst = 1:NO3
                    counts2 += 1
                    MPI_HVNA31Dx[counts2] = -8*real(siT3*coP3*VNA3riαjβ[ist,jst] + coT3*coP3/R*VNA3tiαjβ[ist,jst] - siP3/siT3/R*VNA3piαjβ[ist,jst])
                    MPI_HVNA31Dy[counts2] = -8*real(siT3*siP3*VNA3riαjβ[ist,jst] + coT3*siP3/R*VNA3tiαjβ[ist,jst] + coP3/siT3/R*VNA3piαjβ[ist,jst])
                    MPI_HVNA31Dz[counts2] = -8*real(coT3*VNA3riαjβ[ist,jst] - siT3/R*VNA3tiαjβ[ist,jst])
                end
            end
        else
            for ist = 1:NO2, jst = 1:NO2
                counts1 += 1
                MPI_HVNA21Dx[counts1] = 0.0
                MPI_HVNA21Dy[counts1] = 0.0
                MPI_HVNA21Dz[counts1] = 0.0
            end

            for ist = 1:NO3, jst = 1:NO3
                counts2 += 1
                MPI_HVNA31Dx[counts2] = 0.0
                MPI_HVNA31Dy[counts2] = 0.0
                MPI_HVNA31Dz[counts2] = 0.0
            end
        end
    end
    MPI.Barrier(comm)



    Total_HVNA2size = MPI.Allreduce(myHVNA2size, MPI.SUM, comm) 
    Total_HVNA3size = MPI.Allreduce(myHVNA3size, MPI.SUM, comm) 

    _counts2 = zeros(Int64, nprocs)
    _counts3 = zeros(Int64, nprocs)
    for id = 1:nprocs
        if id-1 == myrank
            _counts2[id] = myHVNA2size
            _counts3[id] = myHVNA3size
        end
        MPI.Barrier(comm)
    end
    MPI.Barrier(comm)

    MPI_HVNA2size = zeros(Int64, nprocs)
    MPI_HVNA3size = zeros(Int64, nprocs)
    MPI.Allreduce!(_counts2, MPI_HVNA2size, nprocs, MPI.SUM, comm)  
    MPI.Allreduce!(_counts3, MPI_HVNA3size, nprocs, MPI.SUM, comm)  


    MPI_HVNA2x = zeros(Float64, Total_HVNA2size)
    MPI_HVNA2y = zeros(Float64, Total_HVNA2size)
    MPI_HVNA2z = zeros(Float64, Total_HVNA2size)
    MPI_HVNA3x = zeros(Float64, Total_HVNA3size)
    MPI_HVNA3y = zeros(Float64, Total_HVNA3size)
    MPI_HVNA3z = zeros(Float64, Total_HVNA3size)


    MPI.Allgatherv!(MPI_HVNA21Dx, VBuffer(MPI_HVNA2x, MPI_HVNA2size), comm)
    MPI.Allgatherv!(MPI_HVNA21Dy, VBuffer(MPI_HVNA2y, MPI_HVNA2size), comm)
    MPI.Allgatherv!(MPI_HVNA21Dz, VBuffer(MPI_HVNA2z, MPI_HVNA2size), comm)
    MPI.Allgatherv!(MPI_HVNA31Dx, VBuffer(MPI_HVNA3x, MPI_HVNA3size), comm)
    MPI.Allgatherv!(MPI_HVNA31Dy, VBuffer(MPI_HVNA3y, MPI_HVNA3size), comm)
    MPI.Allgatherv!(MPI_HVNA31Dz, VBuffer(MPI_HVNA3z, MPI_HVNA3size), comm)
    
    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
        counts += 1
        HVNA2force[1][atom][Rn][ist][jst] = MPI_HVNA2x[counts]
        HVNA2force[2][atom][Rn][ist][jst] = MPI_HVNA2y[counts]
        HVNA2force[3][atom][Rn][ist][jst] = MPI_HVNA2z[counts]
    end
    
    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[natn[atom][Rn]], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        counts += 1
        HVNA3force[1][atom][Rn][ist][jst] = MPI_HVNA3x[counts]
        HVNA3force[2][atom][Rn][ist][jst] = MPI_HVNA3y[counts]
        HVNA3force[3][atom][Rn][ist][jst] = MPI_HVNA3z[counts]
    end
    MPI.Barrier(comm)
end
