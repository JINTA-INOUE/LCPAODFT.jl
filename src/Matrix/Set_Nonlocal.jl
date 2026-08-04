@timeit timer "Set_Nonlocal" function Set_Nonlocal!(SpinPol::AbstractString, NLPforce, HNL, iHNL, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
   
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)


    Natom = system_grid.Natom
    Nspecies = length(pao)
    atom2spe = system_grid.atom2spe
    Total_NumOrbs = system_grid.Total_NumOrbs
    

    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot[spe].Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end

        NLTotal_Num[atom] = tot
    end

    Lmax_Four_Int = 0
    for spe = 1:Nspecies
        Lmax_Four_Int = max(Lmax_Four_Int, maximum(pspot[spe].Spe_VPS_List))
        Lmax_Four_Int = max(Lmax_Four_Int, pao[spe].Spe_MaxL_Basis)
    end


    NLRF_Bessel = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_Num_RVPS = pspot[spe].Spe_Num_RVPS
        Spe_Atom_Cut1 = pao[spe].Spe_Atom_Cut1

        VPS_j_dependency = pspot[spe].VPS_j_dependency
        NLRF_Bessel[spe] = Vector{Vector{Vector{Float64}}}(undef, VPS_j_dependency+1)
        for so = 1:VPS_j_dependency+1
            NLRF_Bessel[spe][so] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
            for l = 1:Spe_Num_RVPS
                NLRF_Bessel[spe][so][l] = zeros(Float64, NkGrid+1)
            end
        end
        FT_NLP!(Spe_Atom_Cut1, NLRF_Bessel[spe], pspot[spe])
    end



    k1 = zeros(Float64, NkGrid+1)
    k2 = zeros(Float64, NkGrid+1)
    k3 = zeros(Float64, NkGrid+1)
    dk = (Nkmax-Radial_kmin)/NkGrid
    for ik = 1:NkGrid+1
        k1[ik] = Radial_kmin + (ik-1)*dk
    end
    @. k2 = k1^2
    @. k3 = k1^3
    k2[begin] = 0.5*k2[begin]
    k2[end] = 0.5*k2[end]
    k3[begin] = 0.5*k3[begin]
    k3[end] = 0.5*k3[end]



    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz



    fsize = maximum(Total_NumOrbs)
    NLfsize = maximum(NLTotal_Num)
    NLPiαjβ = zeros(ComplexF64, fsize, NLfsize)
    # NLPriαjβ = zeros(ComplexF64, fsize, NLfsize)
    # NLPtiαjβ = zeros(ComplexF64, fsize, NLfsize)
    # NLPpiαjβ = zeros(ComplexF64, fsize, NLfsize)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    Cjβ = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Cjβ[spe] = zeros(ComplexF64, NLfsize, NLfsize)
        Set_Comp2Real!(Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis)
        Set_NLComp2Real!(Cjβ[spe], pspot[spe].Spe_VPS_List)
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
    dSphB = Vector{Vector{Float64}}(undef, 2*Lmax_Four_Int+3)
    for l = 1:2*Lmax_Four_Int+3
        SphB[l] = zeros(Float64, NkGrid+1)
        dSphB[l] = zeros(Float64, NkGrid+1)
    end
    SphB2 = zeros(Float64, NkGrid+1)
    dSphB3 = zeros(Float64, NkGrid+1)
    SumNL0 = zeros(Float64, 15, 4, 4)
    # SumNLr0 = zeros(Float64, 15, 4, 4)
    tmpL = zeros(Float64, NkGrid+1)
    tmpH1 = zeros(ComplexF64, fsize)
    tmpH2 = zeros(ComplexF64, NLfsize)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)



    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        iRF_Bessel = pao[ispe].Spe_RF_Bessel

        NO0 = Total_NumOrbs[atom]

        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        cell = MPI_ncn[loop]+1
        NO1 = NLTotal_Num[jatom]

        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]


        R, theta, phi = xyz_to_spherical(x, y, z)
        R = ifelse(R < 1.0e-10, 1.0e-10, R)
        # siT = sin(theta)
        # coT = cos(theta)
        # siP = sin(phi)
        # coP = cos(phi)
            

        jNum_RVPS = pspot[jspe].Spe_Num_RVPS
        jVPS_List = pspot[jspe].Spe_VPS_List
        VPS_j_dependency = pspot[jspe].VPS_j_dependency

        Lmax = maximum(jVPS_List)
        Lmax_Four_Int = 2*ifelse(Lmax>iMaxL_Basis, Lmax, iMaxL_Basis)

        for ik = 1:NkGrid+1
            Calc_SphericalBesselj!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                # dSphB[l][ik] = dSphB_l[l]
            end
        end
        

        for so = 1:VPS_j_dependency+1

            fill!(NLPiαjβ, 0.0)
            # fill!(NLPriαjβ, 0.0)
            # fill!(NLPtiαjβ, 0.0)
            # fill!(NLPpiαjβ, 0.0)

            # Σ_{L=0}^{Lmax_Four_Int}Sum_{M=-L}^{L}
            for L = 0:Lmax_Four_Int

                @. SphB2 = SphB[L+1]*k2
                @. dSphB3 = dSphB[L+1]*k3

                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:jNum_RVPS
                    @. tmpL = iRF_Bessel[l+1][p]*NLRF_Bessel[jspe][so][lnum]
                    SumNL0[lnum,p,l+1] = dot(SphB2, tmpL)*dk
                    # SumNLr0[lnum,p,l+1] = dot(dSphB3, tmpL)*dk
                end



                for M = -L:L
                    ist = 0
                    for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                        jst = 0
                        ist += 1
                        @inbounds for lnum = 1:jNum_RVPS, mm = -jVPS_List[lnum]:jVPS_List[lnum]
                            jst += 1
                            ll = jVPS_List[lnum]
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

                                NLPiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                                # NLPriαjβ[ist,jst] += iYC*SumNLr0[lnum,p,l+1]
                                # NLPtiαjβ[ist,jst] += iYCt*SumNL0[lnum,p,l+1]
                                # NLPpiαjβ[ist,jst] += iYCp*SumNL0[lnum,p,l+1]
                            end
                        end
                    end
                end
            end
            

            # complex to real        
            @inbounds for ist = 1:NO0
                @views mul!(tmpH2, Cjβ[jspe], NLPiαjβ[ist,:])
                @views NLPiαjβ[ist,:] = tmpH2
                # @views mul!(tmpH2, Cjβ[jspe], NLPriαjβ[ist,:])
                # @views NLPriαjβ[ist,:] = tmpH2
                # @views mul!(tmpH2, Cjβ[jspe], NLPtiαjβ[ist,:])
                # @views NLPtiαjβ[ist,:] = tmpH2
                # @views mul!(tmpH2, Cjβ[jspe], NLPpiαjβ[ist,:])
                # @views NLPpiαjβ[ist,:] = tmpH2
            end

            @inbounds for jst = 1:NO1
                @views mul!(tmpH1, Ciα[ispe], NLPiαjβ[:,jst])
                @views NLPiαjβ[:,jst] = tmpH1
                # @views mul!(tmpH1, Ciα[ispe], NLPriαjβ[:,jst])
                # @views NLPriαjβ[:,jst] = tmpH1
                # @views mul!(tmpH1, Ciα[ispe], NLPtiαjβ[:,jst])
                # @views NLPtiαjβ[:,jst] = tmpH1
                # @views mul!(tmpH1, Ciα[ispe], NLPpiαjβ[:,jst])
                # @views NLPpiαjβ[:,jst] = tmpH1
            end
            

            
            NLPforce1 = NLPforce[1][atom][Rn][so]
            # NLPforce2 = NLPforce[2][atom][Rn][so]
            # NLPforce3 = NLPforce[3][atom][Rn][so]
            # NLPforce4 = NLPforce[4][atom][Rn][so]

            @inbounds for ist = 1:NO0, jst = 1:NO1
                NLPforce1[ist,jst] = 8*real(NLPiαjβ[ist,jst])
            end

            #=
            if Rn ≠ 1
                if abs(siT) < 1.0e-13
                    @inbounds for ist = 1:NO0, jst = 1:NO1
                        NLPforce2[ist,jst] = -8*real(siT*coP*NLPriαjβ[ist,jst] + coT*coP/R*NLPtiαjβ[ist,jst])
                        NLPforce3[ist,jst] = -8*real(siT*siP*NLPriαjβ[ist,jst] + coT*siP/R*NLPtiαjβ[ist,jst])
                        NLPforce4[ist,jst] = -8*real(coT*NLPriαjβ[ist,jst] - siT/R*NLPtiαjβ[ist,jst])
                    end
                else
                    @inbounds for ist = 1:NO0, jst = 1:NO1
                        NLPforce2[ist,jst] = -8*real(siT*coP*NLPriαjβ[ist,jst] + coT*coP/R*NLPtiαjβ[ist,jst] - siP/siT/R*NLPpiαjβ[ist,jst])
                        NLPforce3[ist,jst] = -8*real(siT*siP*NLPriαjβ[ist,jst] + coT*siP/R*NLPtiαjβ[ist,jst] + coP/siT/R*NLPpiαjβ[ist,jst])
                        NLPforce4[ist,jst] = -8*real(coT*NLPriαjβ[ist,jst] - siT/R*NLPtiαjβ[ist,jst])
                    end
                end
            else
                @inbounds for ist = 1:NO0, jst = 1:NO1
                    NLPforce2[ist,jst] = 0.0
                    NLPforce3[ist,jst] = 0.0
                    NLPforce4[ist,jst] = 0.0
                end
            end=#
        end
    end


    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        jspe = atom2spe[jatom]
        VPS_j_dependency = pspot[jspe].VPS_j_dependency
        @inbounds for so = 1:VPS_j_dependency+1
            MPI.Allreduce!(NLPforce[1][atom][Rn][so], MPI.SUM, comm)
            # MPI.Allreduce!(NLPforce[2][atom][Rn][so], MPI.SUM, comm)
            # MPI.Allreduce!(NLPforce[3][atom][Rn][so], MPI.SUM, comm)
            # MPI.Allreduce!(NLPforce[4][atom][Rn][so], MPI.SUM, comm)
        end
    end
    
    
    if SpinPol ∈ ("off", "on")

        VNLE = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            spe = atom2spe[atom]
            VNLE[atom] = zeros(Float64, NLTotal_Num[atom])

            counts = 0
            for lst = 1:pspot[spe].Spe_Num_RVPS
                ene = pspot[spe].Spe_VNLE[1,lst]
                L = pspot[spe].Spe_VPS_List[lst]
                for _ = -L:L
                    counts += 1
                    VNLE[atom][counts] = ene
                end
            end
        end

        Set_Nonlocal_Col!(HNL, NLPforce[1], VNLE, NLTotal_Num, system_grid)

    elseif SpinPol == "nc"
        Set_Nonlocal_NonCol!(HNL, iHNL, NLPforce[1], pspot, system_grid)
    end
end


function Set_Nonlocal_Col!(HNL, NLP, VNLE, NLTotal_Num, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    FNAN = system_grid.FNAN
	natn = system_grid.natn
    Atom_Cut1 = system_grid.Atom_Cut1
    Dis = system_grid.Dis
    RMI = system_grid.RMI
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size


    tmpL = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        tmpL[atom] = zeros(Float64, NLTotal_Num[atom])
    end

    max_orbitals = Int(maximum(Total_NumOrbs))
    max_projectors = Int(maximum(NLTotal_Num))
    HNL_temp = zeros(Float64, max_orbitals, max_orbitals)
    weighted_projector = zeros(Float64, max_orbitals, max_projectors)


    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        fill!(HNL_temp, 0.0)
        C = @view HNL_temp[1:NO0, 1:NO1]
        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            if kl >= 0
                A = NLP[atom][Rm][1]
                B = NLP[jatom][kl+1][1]
                nprojectors = Int(NLTotal_Num[kg])
                energies = VNLE[kg]
                @inbounds for projector = 1:nprojectors, jst = 1:NO1
                    weighted_projector[jst, projector] = B[jst, projector] * energies[projector]
                end
                Bweighted = @view weighted_projector[1:NO1, 1:nprojectors]
                mul!(C, A, transpose(Bweighted), 1.0, 1.0)
            end
        end

        HNL1 = HNL[1][atom][Rn]
        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, Dis[atom][Rn])
        @inbounds for ist = 1:NO0, jst = 1:NO1
            HNL1[ist][jst] = dmp * HNL_temp[ist, jst]
        end
    end

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(HNL[1][atom][Rn][ist], MPI.SUM, comm)
    end
end


function Set_Nonlocal_NonCol!(HNL, iHNL, NLP, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    Total_NumOrbs = system_grid.Total_NumOrbs
    FNAN = system_grid.FNAN
	natn = system_grid.natn
    Atom_Cut1 = system_grid.Atom_Cut1
    Dis = system_grid.Dis
    RMI = system_grid.RMI
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    max_orbitals = Int(maximum(Total_NumOrbs))
    
    HNL_temp = zeros(Float64, max_orbitals, max_orbitals, 3)
    iHNL_temp = zeros(Float64, max_orbitals, max_orbitals, 3)


    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        fill!(HNL_temp, 0.0)
        fill!(iHNL_temp, 0.0)

        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kgspe = atom2spe[kg]
            kl = RMI[atom][Rn][Rm]

            Spe_Num_RVPS = pspot[kgspe].Spe_Num_RVPS
            Spe_VNLE = pspot[kgspe].Spe_VNLE
            Spe_VPS_List = pspot[kgspe].Spe_VPS_List 

            if kl >= 0

                NLP11 = NLP[atom][Rm][1]
                NLP12 = NLP[atom][Rm][2]
                NLP21 = NLP[jatom][kl+1][1]
                NLP22 = NLP[jatom][kl+1][2]

                for ist = 1:NO0, jst = 1:NO1

                    Sum0_r = 0.0
                    Sum0_i = 0.0
                    Sum1_r = 0.0
                    Sum1_i = 0.0
                    Sum2_r = 0.0
                    Sum2_i = 0.0
                    
                    L = 1
                    @inbounds for lnum = 1:Spe_Num_RVPS
                        ene_p = Spe_VNLE[1,lnum]
                        ene_m = Spe_VNLE[2,lnum]
                        if Spe_VPS_List[lnum] == 0
                            L2 = 0
                            PFp = 1.0
                            PFm = 0.0
                        elseif Spe_VPS_List[lnum] == 1
                            L2 = 2
                            PFp = 2/3
                            PFm = 1/3
                        elseif Spe_VPS_List[lnum] == 2
                            L2 = 4
                            PFp = 3/5
                            PFm = 2/5
                        elseif Spe_VPS_List[lnum] == 3
                            L2 = 6
                            PFp = 4/7
                            PFm = 3/7
                        else
                            error("not support Spe_VPS_List >= 4")
                        end
                        

                        # off-diagonal contribution on up-dn
                        if L2 == 2
                            Sum2_r += (ene_p/3 * NLP11[ist,L  ] * NLP21[jst,L+2]
                                      -ene_p/3 * NLP11[ist,L+2] * NLP21[jst,L  ])

                            Sum2_i += (-ene_p/3 * NLP11[ist,L+1] * NLP21[jst,L+2]
                                       +ene_p/3 * NLP11[ist,L+2] * NLP21[jst,L+1])


                            Sum2_r -= (ene_m/3 * NLP12[ist,L  ] * NLP22[jst,L+2]
                                      -ene_m/3 * NLP12[ist,L+2] * NLP22[jst,L  ])

                            Sum2_i -= (-ene_m/3 * NLP12[ist,L+1] * NLP22[jst,L+2]
                                       +ene_m/3 * NLP12[ist,L+2] * NLP22[jst,L+1])
                        elseif L2 == 4
                            tmp0 = sqrt(3)
                            tmp1 = ene_p/5
                            tmp2 = tmp0*tmp1

                            Sum2_r += (-tmp2 * NLP11[ist,L  ] * NLP21[jst,L+3]
                                       +tmp2 * NLP11[ist,L+3] * NLP21[jst,L  ]
                                       +tmp1 * NLP11[ist,L+1] * NLP21[jst,L+3]
                                       -tmp1 * NLP11[ist,L+3] * NLP21[jst,L+1]
                                       +tmp1 * NLP11[ist,L+2] * NLP21[jst,L+4]
                                       -tmp1 * NLP11[ist,L+4] * NLP21[jst,L+2])
                                
                            Sum2_i +=  (tmp2 * NLP11[ist,L  ] * NLP21[jst,L+4]
                                       -tmp2 * NLP11[ist,L+4] * NLP21[jst,L  ]
                                       +tmp1 * NLP11[ist,L+1] * NLP21[jst,L+4]
                                       -tmp1 * NLP11[ist,L+4] * NLP21[jst,L+1]
                                       -tmp1 * NLP11[ist,L+2] * NLP21[jst,L+3]
                                       +tmp1 * NLP11[ist,L+3] * NLP21[jst,L+2])

                            tmp1 = ene_m/5
                            tmp2 = tmp0*tmp1

                            Sum2_r -= (-tmp2 * NLP12[ist,L  ] * NLP22[jst,L+3]
                                       +tmp2 * NLP12[ist,L+3] * NLP22[jst,L  ]
                                       +tmp1 * NLP12[ist,L+1] * NLP22[jst,L+3]
                                       -tmp1 * NLP12[ist,L+3] * NLP22[jst,L+1]
                                       +tmp1 * NLP12[ist,L+2] * NLP22[jst,L+4]
                                       -tmp1 * NLP12[ist,L+4] * NLP22[jst,L+2])

                            Sum2_i -= ( tmp2 * NLP12[ist,L  ] * NLP22[jst,L+4]
                                       -tmp2 * NLP12[ist,L+4] * NLP22[jst,L  ]
                                       +tmp1 * NLP12[ist,L+1] * NLP22[jst,L+4]
                                       -tmp1 * NLP12[ist,L+4] * NLP22[jst,L+1]
                                       -tmp1 * NLP12[ist,L+2] * NLP22[jst,L+3]
                                       +tmp1 * NLP12[ist,L+3] * NLP22[jst,L+2])
                        elseif L2 == 6
                            tmp0 = sqrt(6)
                            tmp1 = sqrt(3/2)
                            tmp2 = sqrt(5/2)
                            tmp3 = ene_p/7
                            tmp4 = tmp1*tmp3
                            tmp5 = tmp2*tmp3
                            tmp6 = tmp0*tmp3

                            Sum2_r += (-tmp6*NLP11[ist,L  ] * NLP21[jst,L+1]
                                       +tmp6*NLP11[ist,L+1] * NLP21[jst,L  ]
                                       -tmp5*NLP11[ist,L+1] * NLP21[jst,L+3]
                                       +tmp5*NLP11[ist,L+3] * NLP21[jst,L+1]
                                       -tmp5*NLP11[ist,L+2] * NLP21[jst,L+4]
                                       +tmp5*NLP11[ist,L+4] * NLP21[jst,L+2]
                                       -tmp4*NLP11[ist,L+3] * NLP21[jst,L+5]
                                       +tmp4*NLP11[ist,L+5] * NLP21[jst,L+3]
                                       -tmp4*NLP11[ist,L+4] * NLP21[jst,L+6]
                                       +tmp4*NLP11[ist,L+6] * NLP21[jst,L+4])

                            Sum2_i += ( tmp6*NLP11[ist,L  ] * NLP21[jst,L+2]
                                       -tmp6*NLP11[ist,L+2] * NLP21[jst,L  ]
                                       +tmp5*NLP11[ist,L+1] * NLP21[jst,L+4]
                                       -tmp5*NLP11[ist,L+4] * NLP21[jst,L+1]
                                       -tmp5*NLP11[ist,L+2] * NLP21[jst,L+3]
                                       +tmp5*NLP11[ist,L+3] * NLP21[jst,L+2]
                                       +tmp4*NLP11[ist,L+3] * NLP21[jst,L+6]
                                       -tmp4*NLP11[ist,L+6] * NLP21[jst,L+3]
                                       -tmp4*NLP11[ist,L+4] * NLP21[jst,L+5]
                                       +tmp4*NLP11[ist,L+5] * NLP21[jst,L+4])            

                            tmp3 = ene_m/7
                            tmp4 = tmp1*tmp3
                            tmp5 = tmp2*tmp3
                            tmp6 = tmp0*tmp3

                            Sum2_r -= (-tmp6*NLP12[ist,L  ] * NLP22[jst,L+1]
                                       +tmp6*NLP12[ist,L+1] * NLP22[jst,L  ]
                                       -tmp5*NLP12[ist,L+1] * NLP22[jst,L+3]
                                       +tmp5*NLP12[ist,L+3] * NLP22[jst,L+1]
                                       -tmp5*NLP12[ist,L+2] * NLP22[jst,L+4]
                                       +tmp5*NLP12[ist,L+4] * NLP22[jst,L+2]
                                       -tmp4*NLP12[ist,L+3] * NLP22[jst,L+5]
                                       +tmp4*NLP12[ist,L+5] * NLP22[jst,L+3]
                                       -tmp4*NLP12[ist,L+4] * NLP22[jst,L+6]
                                       +tmp4*NLP12[ist,L+6] * NLP22[jst,L+4])

                            Sum2_i -= ( tmp6*NLP12[ist,L  ] * NLP22[jst,L+2]
                                       -tmp6*NLP12[ist,L+2] * NLP22[jst,L  ]
                                       +tmp5*NLP12[ist,L+1] * NLP22[jst,L+4]
                                       -tmp5*NLP12[ist,L+4] * NLP22[jst,L+1]
                                       -tmp5*NLP12[ist,L+2] * NLP22[jst,L+3]
                                       +tmp5*NLP12[ist,L+3] * NLP22[jst,L+2]
                                       +tmp4*NLP12[ist,L+3] * NLP22[jst,L+6]
                                       -tmp4*NLP12[ist,L+6] * NLP22[jst,L+3]
                                       -tmp4*NLP12[ist,L+4] * NLP22[jst,L+5]
                                       +tmp4*NLP12[ist,L+5] * NLP22[jst,L+4])
                        end

                        # off-diagonal contribution on up-up and dn-dn
                        if L2 == 2
                            tmp0 = (ene_p/3 * NLP11[ist,L  ] * NLP21[jst,L+1]
                                   -ene_p/3 * NLP11[ist,L+1] * NLP21[jst,L  ])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = (ene_m/3 * NLP12[ist,L  ] * NLP22[jst,L+1]
                                   -ene_m/3 * NLP12[ist,L+1] * NLP22[jst,L  ])

                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        elseif L2 == 4
                            tmp0 = ( ene_p*2/5 * NLP11[ist,L+1] * NLP21[jst,L+2]
                                    -ene_p*2/5 * NLP11[ist,L+2] * NLP21[jst,L+1]
                                    +ene_p*1/5 * NLP11[ist,L+3] * NLP21[jst,L+4]
                                    -ene_p*1/5 * NLP11[ist,L+4] * NLP21[jst,L+3])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = ( ene_m*2/5 * NLP12[ist,L+1] * NLP22[jst,L+2]
                                    -ene_m*2/5 * NLP12[ist,L+2] * NLP22[jst,L+1]
                                    +ene_m*1/5 * NLP12[ist,L+3] * NLP22[jst,L+4]
                                    -ene_m*1/5 * NLP12[ist,L+4] * NLP22[jst,L+3])

                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        elseif L2 == 6
                            tmp0 = ( ene_p*1/7 * NLP11[ist,L+1] * NLP21[jst,L+2]
                                    -ene_p*1/7 * NLP11[ist,L+2] * NLP21[jst,L+1]
                                    +ene_p*2/7 * NLP11[ist,L+3] * NLP21[jst,L+4]
                                    -ene_p*2/7 * NLP11[ist,L+4] * NLP21[jst,L+3]
                                    +ene_p*3/7 * NLP11[ist,L+5] * NLP21[jst,L+6]
                                    -ene_p*3/7 * NLP11[ist,L+6] * NLP21[jst,L+5])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = ( ene_m*1/7 * NLP12[ist,L+1] * NLP22[jst,L+2]
                                    -ene_m*1/7 * NLP12[ist,L+2] * NLP22[jst,L+1]
                                    +ene_m*2/7 * NLP12[ist,L+3] * NLP22[jst,L+4]
                                    -ene_m*2/7 * NLP12[ist,L+4] * NLP22[jst,L+3]
                                    +ene_m*3/7 * NLP12[ist,L+5] * NLP22[jst,L+6]
                                    -ene_m*3/7 * NLP12[ist,L+6] * NLP22[jst,L+5])
                                    
                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        end

                        # diagonal contribution on up-up and dn-dn
                        for _ = 0:L2
                            Sum0_r += PFp*ene_p * NLP11[ist,L] * NLP21[jst,L]
                            Sum1_r += PFp*ene_p * NLP11[ist,L] * NLP21[jst,L]

                            Sum0_r += PFm*ene_m * NLP12[ist,L] * NLP22[jst,L]
                            Sum1_r += PFm*ene_m * NLP12[ist,L] * NLP22[jst,L]

                            L += 1
                        end
                    end

                    HNL_temp[ist,jst,1] += Sum0_r
                    HNL_temp[ist,jst,2] += Sum1_r
                    HNL_temp[ist,jst,3] += Sum2_r
                    iHNL_temp[ist,jst,1] += Sum0_i
                    iHNL_temp[ist,jst,2] += Sum1_i
                    iHNL_temp[ist,jst,3] += Sum2_i
                end
            end
        end


        HNL1 = HNL[1][atom][Rn]
        HNL2 = HNL[2][atom][Rn]
        HNL3 = HNL[3][atom][Rn]
        iHNL1 = iHNL[1][atom][Rn]
        iHNL2 = iHNL[2][atom][Rn]
        iHNL3 = iHNL[3][atom][Rn]

        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, Dis[atom][Rn])
        for ist = 1:NO0, jst = 1:NO1
            HNL1[ist][jst] = dmp * HNL_temp[ist,jst,1]
            HNL2[ist][jst] = dmp * HNL_temp[ist,jst,2]
            HNL3[ist][jst] = dmp * HNL_temp[ist,jst,3]
            iHNL1[ist][jst] = dmp * iHNL_temp[ist,jst,1]
            iHNL2[ist][jst] = dmp * iHNL_temp[ist,jst,2]
            iHNL3[ist][jst] = dmp * iHNL_temp[ist,jst,3]
        end
    end


    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(HNL[1][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(HNL[2][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(HNL[3][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iHNL[1][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iHNL[2][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iHNL[3][atom][Rn][ist], MPI.SUM, comm)
    end
end
