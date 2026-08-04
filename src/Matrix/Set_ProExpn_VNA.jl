# @timeit timer "Set_ProExpn_VNA" function Set_ProExpn_VNA!(DS_VNA, HVNA, HVNA2force, HVNA3force, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
@timeit timer "Set_ProExpn_VNA" function Set_ProExpn_VNA!(DS_VNA, HVNA, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    Total_NumOrbs = system_grid.Total_NumOrbs

    Set_ProExpn!(DS_VNA, HVNA, pao, pspot, system_grid)

    HVNA2 = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
	for atom = 1:Natom
		HVNA2[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
		for Rn = 1:FNAN[atom]+1
			HVNA2[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
			for ist = 1:Total_NumOrbs[atom]
				HVNA2[atom][Rn][ist] = zeros(Float64, Total_NumOrbs[atom])
			end
		end
	end
    # Set_HVNA2_3!(HVNA2, HVNA2force, HVNA3force, pao, pspot, system_grid)
    Set_HVNA2_3!(HVNA2, pao, pspot, system_grid)

    for atom = 1:Natom
        NO0 = Total_NumOrbs[atom]
        @inbounds for ist = 1:NO0, jst = 1:NO0
            HVNA[atom][1][ist][jst] = 0.0
        end

        @inbounds for Rn = 1:FNAN[atom]+1, ist = 1:NO0, jst = 1:NO0
            HVNA[atom][1][ist][jst] += HVNA2[atom][Rn][ist][jst]
        end
    end
end


function Set_ProExpn!(DS_VNA, HVNA, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

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


    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size
    



    SumNL0 = zeros(Float64, Num_RVNA, 4, 4)
    # SumNLr0 = zeros(Float64, Num_RVNA, 4, 4)
    fsize = maximum(Total_NumOrbs)
    VNAiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    # VNAriαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    # VNAtiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    # VNApiαjβ = zeros(ComplexF64, fsize, VNATotal_Num)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!(Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis)
        conj!(Ciα[spe])
    end
    Cjβ = Set_VNAComp2Real(maxL)


    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    fact2 = zeros(Float64, 2*Lmax_Four_Int+1, 2*Lmax_Four_Int+1)
    Set_SqrtFactorial_Ratio!(fact2)


    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 30)
    dSphB_l = zeros(Float64, 30)
    SphB = Vector{Vector{Float64}}(undef, Lmax_Four_Int+1)
    # dSphB = Vector{Vector{Float64}}(undef, Lmax_Four_Int+1)
    for l = 1:Lmax_Four_Int+1
        SphB[l] = zeros(Float64, GL_Mesh)
        # dSphB[l] = zeros(Float64, GL_Mesh)
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
        # siT = sin(theta)
        # coT = cos(theta)
        # siP = sin(phi)
        # coP = cos(phi)
            

        for ik = 1:GL_Mesh
            Calc_SphericalBesselj2!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            @inbounds for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                # dSphB[l][ik] = dSphB_l[l]*k1[ik]
            end
        end


        fill!(VNAiαjβ, 0.0)
        # fill!(VNAriαjβ, 0.0)
        # fill!(VNAtiαjβ, 0.0)
        # fill!(VNApiαjβ, 0.0)


        for L = 0:Lmax_Four_Int

            for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:Num_RVNA
                ll = VNA_List[lnum]
                index = VNA_List2[lnum]+1
                @. tmpL = Bessel_Pro00[ispe][l+1][p]*VNA_Bessel[jspe][ll+1][index]
                SumNL0[lnum,p,l+1] = dot(SphB[L+1], tmpL)
                # SumNLr0[lnum,p,l+1] = dot(dSphB[L+1], tmpL)
            end


            for M = -L:L
                ist = 0
                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                    jst = 0
                    ist += 1
                    @inbounds for lnum = 1:Num_RVNA, mm = -VNA_List[lnum]:VNA_List[lnum]
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
                            # dYlmdtheta = ComplexF64(dSHt[1], dSHt[2])
                            # dYlmdphi = ComplexF64(dSHp[1], dSHp[2])
                            
                            iYC = conj(Ylm) * tmp * gaunt
                            # iYCt = conj(dYlmdtheta) * tmp * gaunt
                            # iYCp = conj(dYlmdphi) * tmp * gaunt

                            VNAiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                            # VNAriαjβ[ist,jst] += iYC*SumNLr0[lnum,p,l+1]
                            # VNAtiαjβ[ist,jst] += iYCt*SumNL0[lnum,p,l+1]
                            # VNApiαjβ[ist,jst] += iYCp*SumNL0[lnum,p,l+1]
                        end
                    end
                end
            end
        end

            
            
        # complex to real
        for ist = 1:NO0
            tot = 1
            @inbounds for lnum = 1:Num_RVNA
                ll = VNA_List[lnum]
                @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAiαjβ[ist,tot:tot+2*ll])
                @views VNAiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                # @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAriαjβ[ist,tot:tot+2*ll])
                # @views VNAriαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                # @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNAtiαjβ[ist,tot:tot+2*ll])
                # @views VNAtiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                # @views mul!(tmpH1[ll+1], Cjβ[ll+1], VNApiαjβ[ist,tot:tot+2*ll])
                # @views VNApiαjβ[ist,tot:tot+2*ll] = tmpH1[ll+1]
                tot += 2*ll+1
            end
        end


        @inbounds for jst = 1:VNATotal_Num
            @views mul!(tmpH, Ciα[ispe], VNAiαjβ[:,jst])
            @views VNAiαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNAriαjβ[:,jst])
            # @views VNAriαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNAtiαjβ[:,jst])
            # @views VNAtiαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNApiαjβ[:,jst])
            # @views VNApiαjβ[:,jst] = tmpH
        end



        DS_VNA1 = DS_VNA[1][atom][Rn]
        # DS_VNA2 = DS_VNA[2][atom][Rn]
        # DS_VNA3 = DS_VNA[3][atom][Rn]
        # DS_VNA4 = DS_VNA[4][atom][Rn]

        @inbounds for ist = 1:NO0, jst = 1:VNATotal_Num
            DS_VNA1[ist,jst] = 8*real(VNAiαjβ[ist,jst])
        end

        #=
        if Rn ≠ 1
            if abs(siT) < 1.0e-13
                @inbounds for ist = 1:NO0, jst = 1:VNATotal_Num
                    DS_VNA2[ist,jst] = -8*real(siT*coP*VNAriαjβ[ist,jst] + coT*coP/R*VNAtiαjβ[ist,jst])
                    DS_VNA3[ist,jst] = -8*real(siT*siP*VNAriαjβ[ist,jst] + coT*siP/R*VNAtiαjβ[ist,jst])
                    DS_VNA4[ist,jst] = -8*real(coT*VNAriαjβ[ist,jst] - siT/R*VNAtiαjβ[ist,jst])
                end
            else
                @inbounds for ist = 1:NO0, jst = 1:VNATotal_Num
                    DS_VNA2[ist,jst] = -8*real(siT*coP*VNAriαjβ[ist,jst] + coT*coP/R*VNAtiαjβ[ist,jst] - siP/siT/R*VNApiαjβ[ist,jst])
                    DS_VNA3[ist,jst] = -8*real(siT*siP*VNAriαjβ[ist,jst] + coT*siP/R*VNAtiαjβ[ist,jst] + coP/siT/R*VNApiαjβ[ist,jst])
                    DS_VNA4[ist,jst] = -8*real(coT*VNAriαjβ[ist,jst] - siT/R*VNAtiαjβ[ist,jst])
                end
            end
        else
            @inbounds for ist = 1:NO0, jst = 1:VNATotal_Num
                DS_VNA2[ist,jst] = 0.0
                DS_VNA3[ist,jst] = 0.0
                DS_VNA4[ist,jst] = 0.0
            end
        end=#
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

    #=
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        @inbounds for ist = 1:Total_NumOrbs[atom], jst = 1:VNATotal_Num
            ene = VNAE[jatom][jst]
            DS_VNA[2][atom][Rn][ist,jst] = ene*DS_VNA[2][atom][Rn][ist,jst]
            DS_VNA[3][atom][Rn][ist,jst] = ene*DS_VNA[3][atom][Rn][ist,jst]
            DS_VNA[4][atom][Rn][ist,jst] = ene*DS_VNA[4][atom][Rn][ist,jst]
        end
    end
    =#


    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        MPI.Allreduce!(DS_VNA[1][atom][Rn], MPI.SUM, comm)
        # MPI.Allreduce!(DS_VNA[2][atom][Rn], MPI.SUM, comm)
        # MPI.Allreduce!(DS_VNA[3][atom][Rn], MPI.SUM, comm)
        # MPI.Allreduce!(DS_VNA[4][atom][Rn], MPI.SUM, comm)
    end


    Set_HVNA!(VNAE, DS_VNA[1], HVNA, system_grid)
end


function Set_HVNA!(VNAE, DS_VNA, HVNA, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    FNAN = system_grid.FNAN
	natn = system_grid.natn
    Dis = system_grid.Dis
    RMI = system_grid.RMI
    Atom_Cut1 = system_grid.Atom_Cut1
    VNATotal_Num = length(VNAE[begin])

    max_orbitals = maximum(Total_NumOrbs)
    HVNA_temp = zeros(Float64, max_orbitals, max_orbitals)
    weighted_projector = zeros(Float64, max_orbitals, VNATotal_Num)


    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        fill!(HVNA_temp, 0.0)
        C = @view HVNA_temp[1:NO0, 1:NO1]

        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            if kl >= 0
                A = DS_VNA[atom][Rm]
                B = DS_VNA[jatom][kl+1]
                energies = VNAE[kg]
                @inbounds for projector = 1:VNATotal_Num, jst = 1:NO1
                    weighted_projector[jst, projector] = B[jst, projector] * energies[projector]
                end
                Bweighted = @view weighted_projector[1:NO1, 1:VNATotal_Num]
                mul!(C, A, transpose(Bweighted), 1.0, 1.0)
            end
        end

        _HVNA = HVNA[atom][Rn]
        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, Dis[atom][Rn])
        @inbounds for ist = 1:NO0, jst = 1:NO1
            _HVNA[ist][jst] = dmp * HVNA_temp[ist, jst]
        end
    end

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(HVNA[atom][Rn][ist], MPI.SUM, comm)
    end
end


# function Set_HVNA2_3!(HVNA2, HVNA2force, HVNA3force, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
function Set_HVNA2_3!(HVNA2, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    
    Nspecies = length(pao)


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


    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size

     

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
    # SumHVNAr2 = zeros(Float64, 9, 4, 4, 4, 4)
    # SumHVNA3 = zeros(Float64, 9, 4, 4, 4, 4)
    # SumHVNAr3 = zeros(Float64, 9, 4, 4, 4, 4)


    fsize = maximum(Total_NumOrbs)
    VNA2iαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA2riαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA2tiαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA2piαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA3riαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA3tiαjβ = zeros(ComplexF64, fsize, fsize)
    # VNA3piαjβ = zeros(ComplexF64, fsize, fsize)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    Cjβ = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!(Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis)
        Cjβ[spe] = deepcopy(Ciα[spe])
        conj!(Ciα[spe])
    end


    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    fact2 = zeros(Float64, (2*Lmax_Four_Int+1)^2+1, (2*Lmax_Four_Int+1)^2+1)
    Set_SqrtFactorial_Ratio!(fact2)


    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    dSphB_l = zeros(Float64, 2*Lmax_Four_Int+3)
    SphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    # dSphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    for l = 1:2*Lmax_Four_Int+3
        SphB[l] = zeros(Float64, GL_Mesh)
        # dSphB[l] = zeros(Float64, GL_Mesh)
    end
    tempL = zeros(Float64, GL_Mesh)
    tmpH = zeros(ComplexF64, fsize)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)




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
                # dSphB[l][ik] = dSphB_l[l]*k1[ik]
            end
        end

            
        for l = 0:MaxL_Basis2, p = 1:Num_Basis2[l+1], ll = 0:MaxL_Basis2, pp = 1:Num_Basis2[ll+1]
            if l <= ll
                Lmax_Four_Int = 2*ll
                for L = 0:Lmax_Four_Int
                    if abs(ll-L) <= l <= ll+L
                        @. tempL = Spe_ProductRF_Bessel[ispe][l+1][p][ll+1][pp][L+1]*Spe_CrudeVNA_Bessel[jspe]
                        @views SumHVNA2[L+1,pp,ll+1,p,l+1] = dot(SphB[L+1], tempL)
                        # @views SumHVNAr2[L+1,pp,ll+1,p,l+1] = dot(dSphB[L+1], tempL)
                    end
                end
            end
        end


        #=
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
        =#

        fill!(VNA2iαjβ, 0.0)
        # fill!(VNA2riαjβ, 0.0)
        # fill!(VNA2tiαjβ, 0.0)
        # fill!(VNA2piαjβ, 0.0)
        # fill!(VNA3riαjβ, 0.0)
        # fill!(VNA3tiαjβ, 0.0)
        # fill!(VNA3piαjβ, 0.0)


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
                            # dYlmdtheta = ComplexF64(dSHt[1], dSHt[2])
                            # dYlmdphi = ComplexF64(dSHp[1], dSHp[2])
                                
                            YC = Ylm * gaunt
                            # YCt = dYlmdtheta * gaunt
                            # YCp = dYlmdphi * gaunt

                            VNA2iαjβ[ist,jst] += YC*SumHVNA2[L+1,pp,ll+1,p,l+1]
                            # VNA2riαjβ[ist,jst] += YC*SumHVNAr2[L+1,pp,ll+1,p,l+1]
                            # VNA2tiαjβ[ist,jst] += YCt*SumHVNA2[L+1,pp,ll+1,p,l+1]
                            # VNA2piαjβ[ist,jst] += YCp*SumHVNA2[L+1,pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end


        #=
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
        =#


        ist = 0
        for l = 0:MaxL_Basis2, p = 1:Num_Basis2[l+1], m = -l:l
            jst = 0
            ist += 1
            for ll = 0:MaxL_Basis2, pp = 1:Num_Basis2[ll+1], mm = -ll:ll
                jst += 1
                if l <= ll
                    VNA2iαjβ[jst,ist] = conj(VNA2iαjβ[ist,jst])
                    # VNA2riαjβ[jst,ist] = conj(VNA2riαjβ[ist,jst])
                    # VNA2tiαjβ[jst,ist] = conj(VNA2tiαjβ[ist,jst])
                    # VNA2piαjβ[jst,ist] = conj(VNA2piαjβ[ist,jst])
                end
            end
        end

        #=
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
        end=#


        # complex to real
        @inbounds for ist = 1:NO2
            @views mul!(tmpH, Cjβ[ispe], VNA2iαjβ[ist,:])
            @views VNA2iαjβ[ist,:] = tmpH
            # @views mul!(tmpH, Cjβ[ispe], VNA2riαjβ[ist,:])
            # @views VNA2riαjβ[ist,:] = tmpH
            # @views mul!(tmpH, Cjβ[ispe], VNA2tiαjβ[ist,:])
            # @views VNA2tiαjβ[ist,:] = tmpH
            # @views mul!(tmpH, Cjβ[ispe], VNA2piαjβ[ist,:])
            # @views VNA2piαjβ[ist,:] = tmpH
        end

        @inbounds for jst = 1:NO2
            @views mul!(tmpH, Ciα[ispe], VNA2iαjβ[:,jst])
            @views VNA2iαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNA2riαjβ[:,jst])
            # @views VNA2riαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNA2tiαjβ[:,jst])
            # @views VNA2tiαjβ[:,jst] = tmpH
            # @views mul!(tmpH, Ciα[ispe], VNA2piαjβ[:,jst])
            # @views VNA2piαjβ[:,jst] = tmpH
        end

        #=
        @inbounds for ist = 1:NO3
            @views mul!(tmpH, Cjβ[jspe], VNA3riαjβ[ist,:])
            @views VNA3riαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], VNA3tiαjβ[ist,:])
            @views VNA3tiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], VNA3piαjβ[ist,:])
            @views VNA3piαjβ[ist,:] = tmpH
        end
    
        @inbounds for jst = 1:NO3
            @views mul!(tmpH, Ciα[jspe], VNA3riαjβ[:,jst])
            @views VNA3riαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[jspe], VNA3tiαjβ[:,jst])
            @views VNA3tiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[jspe], VNA3piαjβ[:,jst])
            @views VNA3piαjβ[:,jst] = tmpH
        end
        =#

        


        HVNA20 = HVNA2[atom][Rn]
        # HVNA2force1 = HVNA2force[1][atom][Rn]
        # HVNA2force2 = HVNA2force[2][atom][Rn]
        # HVNA2force3 = HVNA2force[3][atom][Rn]
        # HVNA3force1 = HVNA3force[1][atom][Rn]
        # HVNA3force2 = HVNA3force[2][atom][Rn]
        # HVNA3force3 = HVNA3force[3][atom][Rn]


        @inbounds for ist = 1:NO2, jst = 1:NO2
            HVNA20[ist][jst] = 8*real(VNA2iαjβ[ist,jst])
        end

        #=
		if Rn ≠ 1
            if abs(siT2) < 1.0e-13
                @inbounds for ist = 1:NO2, jst = 1:NO2
                    HVNA2force1[ist][jst] = -8*real(siT2*coP2*VNA2riαjβ[ist,jst] + coT2*coP2/R*VNA2tiαjβ[ist,jst])
                    HVNA2force2[ist][jst] = -8*real(siT2*siP2*VNA2riαjβ[ist,jst] + coT2*siP2/R*VNA2tiαjβ[ist,jst])
                    HVNA2force3[ist][jst] = -8*real(coT2*VNA2riαjβ[ist,jst] - siT2/R*VNA2tiαjβ[ist,jst])
                end
            else
                @inbounds for ist = 1:NO2, jst = 1:NO2
                    HVNA2force1[ist][jst] = -8*real(siT2*coP2*VNA2riαjβ[ist,jst] + coT2*coP2/R*VNA2tiαjβ[ist,jst] - siP2/siT2/R*VNA2piαjβ[ist,jst])
                    HVNA2force2[ist][jst] = -8*real(siT2*siP2*VNA2riαjβ[ist,jst] + coT2*siP2/R*VNA2tiαjβ[ist,jst] + coP2/siT2/R*VNA2piαjβ[ist,jst])
                    HVNA2force3[ist][jst] = -8*real(coT2*VNA2riαjβ[ist,jst] - siT2/R*VNA2tiαjβ[ist,jst])
                end
            end

            if abs(siT3) < 1.0e-13
                @inbounds for ist = 1:NO3, jst = 1:NO3
                    HVNA3force1[ist][jst] = -8*real(siT3*coP3*VNA3riαjβ[ist,jst] + coT3*coP3/R*VNA3tiαjβ[ist,jst])
                    HVNA3force2[ist][jst] = -8*real(siT3*siP3*VNA3riαjβ[ist,jst] + coT3*siP3/R*VNA3tiαjβ[ist,jst])
                    HVNA3force3[ist][jst] = -8*real(coT3*VNA3riαjβ[ist,jst] - siT3/R*VNA3tiαjβ[ist,jst])
                end
            else
                @inbounds for ist = 1:NO3, jst = 1:NO3
                    HVNA3force1[ist][jst] = -8*real(siT3*coP3*VNA3riαjβ[ist,jst] + coT3*coP3/R*VNA3tiαjβ[ist,jst] - siP3/siT3/R*VNA3piαjβ[ist,jst])
                    HVNA3force2[ist][jst] = -8*real(siT3*siP3*VNA3riαjβ[ist,jst] + coT3*siP3/R*VNA3tiαjβ[ist,jst] + coP3/siT3/R*VNA3piαjβ[ist,jst])
                    HVNA3force3[ist][jst] = -8*real(coT3*VNA3riαjβ[ist,jst] - siT3/R*VNA3tiαjβ[ist,jst])
                end
            end
        else
            @inbounds for ist = 1:NO2, jst = 1:NO2
                HVNA2force1[ist][jst] = 0.0
                HVNA2force2[ist][jst] = 0.0
                HVNA2force3[ist][jst] = 0.0
            end

            @inbounds for ist = 1:NO3, jst = 1:NO3
                HVNA3force1[ist][jst] = 0.0
                HVNA3force2[ist][jst] = 0.0
                HVNA3force3[ist][jst] = 0.0
            end
        end
        =#
    end




    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(HVNA2[atom][Rn][ist], MPI.SUM, comm)
        # MPI.Allreduce!(HVNA2force[1][atom][Rn][ist], MPI.SUM, comm)
        # MPI.Allreduce!(HVNA2force[2][atom][Rn][ist], MPI.SUM, comm)
        # MPI.Allreduce!(HVNA2force[3][atom][Rn][ist], MPI.SUM, comm)
    end

    #=
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[natn[atom][Rn]]
        MPI.Allreduce!(HVNA3force[1][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(HVNA3force[2][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(HVNA3force[3][atom][Rn][ist], MPI.SUM, comm)
    end=#
end
