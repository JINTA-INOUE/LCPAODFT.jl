@timeit timer "Set_Nonlocal" function Set_Nonlocal!(SpinPol::AbstractString, HNL, iHNL, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
   
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

    Lmax = 0
    for spe = 1:Nspecies
        Lmax = max(Lmax, maximum(pspot[spe].Spe_VPS_List))
    end
    for spe = 1:Nspecies
        Lmax = max(Lmax, pao[spe].Spe_MaxL_Basis)
    end
    Lmax = Lmax+1


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
    dk = (Nkmax-Radial_kmin)/NkGrid
    for ik = 1:NkGrid+1
        k1[ik] = Radial_kmin + (ik-1)*dk
    end
    @. k2 = k1^2
    k2[begin] = 0.5*k2[begin]
    k2[end] = 0.5*k2[end]

    


    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz


    MPI_NLPsize = zeros(Int64, nprocs)
    MP_NLP = zeros(Int64, nprocs)

    myNLPsize = 0
    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        VPS_j_dependency = pspot[jspe].VPS_j_dependency
        for so = 1:VPS_j_dependency+1, ist = 1:Total_NumOrbs[atom], jst = 1:NLTotal_Num[jatom]
            myNLPsize += 1
        end
    end

    MPI_NLPsize[myrank+1] = myNLPsize
    Total_NLPsize = MPI.Allreduce(myNLPsize, MPI.SUM, comm)
    MPI.Allreduce!(MPI_NLPsize, MPI.SUM, comm)

    Sum = 0
    for id = 1:nprocs
        MP_NLP[id] = Sum
        Sum += MPI_NLPsize[id]
    end
    NLP_Num = MP_NLP[myrank+1]
    MPI_NLP = zeros(Float64, Total_NLPsize)



    


    fsize = maximum(Total_NumOrbs)
    NLfsize = maximum(NLTotal_Num)
    NLPiαjβ = zeros(ComplexF64, fsize, NLfsize)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    Cjβ = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Cjβ[spe] = zeros(ComplexF64, NLfsize, NLfsize)
        Set_Comp2Real!( Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis )
        Set_NLComp2Real!( Cjβ[spe], pspot[spe].Spe_VPS_List )
        conj!(Ciα[spe])
    end


    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    
    SphB = zeros(Float64, NkGrid+1)
    SphB2 = zeros(Float64, NkGrid+1)
    SumNL0 = zeros(Float64, 15, 4, 4)
    tmpL = zeros(Float64, NkGrid+1)
    tmpH1 = zeros(ComplexF64, fsize)
    tmpH2 = zeros(ComplexF64, NLfsize)



    hst = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        ispe = atom2spe[atom]
        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        iRF_Bessel = pao[ispe].Spe_RF_Bessel

        jatom = MPI_natn[loop]
        cell = MPI_ncn[loop]+1
        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]
        R = sqrt(x^2 + y^2 + z^2)
        R = ifelse(R < 1e-10, 1e-10, R)

        NO0 = Total_NumOrbs[atom]
        NO1 = NLTotal_Num[jatom]

        jspe = atom2spe[jatom]
        jNum_RVPS = pspot[jspe].Spe_Num_RVPS
        jVPS_List = pspot[jspe].Spe_VPS_List

        VPS_j_dependency = pspot[jspe].VPS_j_dependency

        Lmax = maximum(jVPS_List)
        Lmax_Four_Int = 2*max(Lmax, iMaxL_Basis)

        for so = 1:VPS_j_dependency+1

            fill!(NLPiαjβ, 0.0)

            # Σ_{L=0}^{Lmax_Four_Int}Sum_{M=-L}^{L}
            for L = 0:Lmax_Four_Int

                @. SphB = SphericalBesselj(L, R*k1)
                @. SphB2 = SphB*k2

                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:jNum_RVPS
                    @. tmpL = iRF_Bessel[l+1][p]*NLRF_Bessel[jspe][so][lnum]
                    SumNL0[lnum,p,l+1] = dot(SphB2, tmpL)*dk
                end

                for M = -L:L
                    ist = 0
                    for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                        jst = 0
                        ist += 1
                        for lnum = 1:jNum_RVPS, mm = -jVPS_List[lnum]:jVPS_List[lnum]
                            jst += 1
                            ll = jVPS_List[lnum]
                            if abs(ll-L) <= l <= abs(ll+L) && iszero(m-mm-M) && abs(m-M) <= ll
                                Ylm = Ylm_complex(L,M,x,y,z)
                                Ls = Float64(L+ll-l)
                                iYC = (-im)^Ls * conj(Ylm) * Gaunt(f,l,m,ll,mm,L,M)
                                
                                NLPiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                            end
                        end
                    end
                end
            end
            

            # complex to real        
            for ist = 1:NO0
                @views mul!(tmpH2, Cjβ[jspe], NLPiαjβ[ist,:])
                @. @views NLPiαjβ[ist,:] = tmpH2
            end

            for jst = 1:NO1
                @views mul!(tmpH1, Ciα[ispe], NLPiαjβ[:,jst])
                @. @views NLPiαjβ[:,jst] = tmpH1
            end
            
            for ist = 1:NO0, jst = 1:NO1
                hst += 1
                MPI_NLP[NLP_Num+hst] = 8*real(NLPiαjβ[ist,jst])
            end
        end
    end

    
    MPI.Allreduce!(MPI_NLP, MPI.SUM, comm) 


    FNAN = system_grid.FNAN
    natn = system_grid.natn

    NLP = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Natom)
	for atom = 1:Natom
		NLP[atom] = Vector{Vector{Vector{Vector{Float64}}}}(undef, FNAN[atom]+1)
		for Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            jspe = atom2spe[jatom]
            VPS_j_dependency = pspot[jspe].VPS_j_dependency
			NLP[atom][Rn] = Vector{Vector{Vector{Float64}}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                NLP[atom][Rn][ist] = Vector{Vector{Float64}}(undef, VPS_j_dependency+1)
                for so = 1:VPS_j_dependency+1
				    NLP[atom][Rn][ist][so] = zeros(Float64, NLTotal_Num[jatom])
                end
			end
		end
	end

    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        jspe = atom2spe[jatom]
        VPS_j_dependency = pspot[jspe].VPS_j_dependency
        for so = 1:VPS_j_dependency+1, ist = 1:Total_NumOrbs[atom], jst = 1:NLTotal_Num[jatom]
            hst += 1
            NLP[atom][Rn][ist][so][jst] = MPI_NLP[hst]
        end
    end
    MPI.Barrier(comm)
    
    
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

        Set_Nonlocal_Col!(HNL, NLP, VNLE, NLTotal_Num, system_grid)

    elseif SpinPol == "nc"
        Set_Nonlocal_NonCol!(HNL, iHNL, NLP, pspot, system_grid)
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

    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    Dis = system_grid.Dis
    RMI = system_grid.RMI

    MPI_size = system_grid.MPI_size
    MPHks = system_grid.MPHks
    Hks_Num = MPHks[myrank+1]


    tmpL = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        tmpL[atom] = zeros(Float64, NLTotal_Num[atom])
    end


    HNL_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs))
    hst = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        fill!(HNL_temp, 0.0)

        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            if kl >= 0
                @inbounds for jst = 1:NO1, ist = 1:NO0
                    @. tmpL[kg] = NLP[jatom][kl+1][jst][1]*VNLE[kg]
                    HNL_temp[ist,jst] += dot(NLP[atom][Rm][ist][1], tmpL[kg])
                end
            end
        end

        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, Dis[atom][Rn])
        for ist = 1:NO0, jst = 1:NO1
            hst += 1
            HNL[1][Hks_Num+hst] = dmp * HNL_temp[ist,jst]
        end
    end

    MPI.Allreduce!(HNL[1], MPI.SUM, comm) 
end


function Set_Nonlocal_NonCol!(HNL, iHNL, NLP, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    
    atom2spe = system_grid.atom2spe
    Total_NumOrbs = system_grid.Total_NumOrbs
    FNAN = system_grid.FNAN
	natn = system_grid.natn
    Atom_Cut1 = system_grid.Atom_Cut1

    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    Dis = system_grid.Dis
    RMI = system_grid.RMI

    MPI_size = system_grid.MPI_size
    MPHks = system_grid.MPHks
    Hks_Num = MPHks[myrank+1]

    
    HNL_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 3)
    iHNL_temp = zeros(Float64, maximum(Total_NumOrbs), maximum(Total_NumOrbs), 3)

    hst = 0
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
                            Sum2_r += (ene_p/3 * NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+2]
                                      -ene_p/3 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L  ])

                            Sum2_i += (-ene_p/3 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+2]
                                       +ene_p/3 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+1])


                            Sum2_r -= (ene_m/3 * NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+2]
                                      -ene_m/3 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L  ])

                            Sum2_i -= (-ene_m/3 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+2]
                                       +ene_m/3 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+1])
                        elseif L2 == 4
                            tmp0 = sqrt(3)
                            tmp1 = ene_p/5
                            tmp2 = tmp0*tmp1

                            Sum2_r += (-tmp2 * NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+3]
                                       +tmp2 * NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L  ]
                                       +tmp1 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+3]
                                       -tmp1 * NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+1]
                                       +tmp1 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+4]
                                       -tmp1 * NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+2])
                                
                            Sum2_i +=  (tmp2 * NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+4]
                                       -tmp2 * NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L  ]
                                       +tmp1 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+4]
                                       -tmp1 * NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+1]
                                       -tmp1 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+3]
                                       +tmp1 * NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+2])

                            tmp1 = ene_m/5
                            tmp2 = tmp0*tmp1

                            Sum2_r -= (-tmp2 * NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+3]
                                       +tmp2 * NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L  ]
                                       +tmp1 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+3]
                                       -tmp1 * NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+1]
                                       +tmp1 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+4]
                                       -tmp1 * NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+2])

                            Sum2_i -= ( tmp2 * NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+4]
                                       -tmp2 * NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L  ]
                                       +tmp1 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+4]
                                       -tmp1 * NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+1]
                                       -tmp1 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+3]
                                       +tmp1 * NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+2])
                        elseif L2 == 6
                            tmp0 = sqrt(6)
                            tmp1 = sqrt(3/2)
                            tmp2 = sqrt(5/2)
                            tmp3 = ene_p/7
                            tmp4 = tmp1*tmp3
                            tmp5 = tmp2*tmp3
                            tmp6 = tmp0*tmp3

                            Sum2_r += (-tmp6*NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+1]
                                       +tmp6*NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L  ]
                                       -tmp5*NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+3]
                                       +tmp5*NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+1]
                                       -tmp5*NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+4]
                                       +tmp5*NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+2]
                                       -tmp4*NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+5]
                                       +tmp4*NLP[atom][Rm][ist][1][L+5] * NLP[jatom][kl+1][jst][1][L+3]
                                       -tmp4*NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+6]
                                       +tmp4*NLP[atom][Rm][ist][1][L+6] * NLP[jatom][kl+1][jst][1][L+4])

                            Sum2_i += ( tmp6*NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+2]
                                       -tmp6*NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L  ]
                                       +tmp5*NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+4]
                                       -tmp5*NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+1]
                                       -tmp5*NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+3]
                                       +tmp5*NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+2]
                                       +tmp4*NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+6]
                                       -tmp4*NLP[atom][Rm][ist][1][L+6] * NLP[jatom][kl+1][jst][1][L+3]
                                       -tmp4*NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+5]
                                       +tmp4*NLP[atom][Rm][ist][1][L+5] * NLP[jatom][kl+1][jst][1][L+4])            

                            tmp3 = ene_m/7
                            tmp4 = tmp1*tmp3
                            tmp5 = tmp2*tmp3
                            tmp6 = tmp0*tmp3

                            Sum2_r -= (-tmp6*NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+1]
                                       +tmp6*NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L  ]
                                       -tmp5*NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+3]
                                       +tmp5*NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+1]
                                       -tmp5*NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+4]
                                       +tmp5*NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+2]
                                       -tmp4*NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+5]
                                       +tmp4*NLP[atom][Rm][ist][2][L+5] * NLP[jatom][kl+1][jst][2][L+3]
                                       -tmp4*NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+6]
                                       +tmp4*NLP[atom][Rm][ist][2][L+6] * NLP[jatom][kl+1][jst][2][L+4])

                            Sum2_i -= ( tmp6*NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+2]
                                       -tmp6*NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L  ]
                                       +tmp5*NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+4]
                                       -tmp5*NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+1]
                                       -tmp5*NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+3]
                                       +tmp5*NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+2]
                                       +tmp4*NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+6]
                                       -tmp4*NLP[atom][Rm][ist][2][L+6] * NLP[jatom][kl+1][jst][2][L+3]
                                       -tmp4*NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+5]
                                       +tmp4*NLP[atom][Rm][ist][2][L+5] * NLP[jatom][kl+1][jst][2][L+4])
                        end

                        # off-diagonal contribution on up-up and dn-dn
                        if L2 == 2
                            tmp0 = (ene_p/3 * NLP[atom][Rm][ist][1][L  ] * NLP[jatom][kl+1][jst][1][L+1]
                                   -ene_p/3 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L  ])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = (ene_m/3 * NLP[atom][Rm][ist][2][L  ] * NLP[jatom][kl+1][jst][2][L+1]
                                   -ene_m/3 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L  ])

                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        elseif L2 == 4
                            tmp0 = ( ene_p*2/5 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+2]
                                    -ene_p*2/5 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+1]
                                    +ene_p*1/5 * NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+4]
                                    -ene_p*1/5 * NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+3])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = ( ene_m*2/5 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+2]
                                    -ene_m*2/5 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+1]
                                    +ene_m*1/5 * NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+4]
                                    -ene_m*1/5 * NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+3])

                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        elseif L2 == 6
                            tmp0 = ( ene_p*1/7 * NLP[atom][Rm][ist][1][L+1] * NLP[jatom][kl+1][jst][1][L+2]
                                    -ene_p*1/7 * NLP[atom][Rm][ist][1][L+2] * NLP[jatom][kl+1][jst][1][L+1]
                                    +ene_p*2/7 * NLP[atom][Rm][ist][1][L+3] * NLP[jatom][kl+1][jst][1][L+4]
                                    -ene_p*2/7 * NLP[atom][Rm][ist][1][L+4] * NLP[jatom][kl+1][jst][1][L+3]
                                    +ene_p*3/7 * NLP[atom][Rm][ist][1][L+5] * NLP[jatom][kl+1][jst][1][L+6]
                                    -ene_p*3/7 * NLP[atom][Rm][ist][1][L+6] * NLP[jatom][kl+1][jst][1][L+5])

                            Sum0_i += -tmp0
                            Sum1_i += tmp0

                            tmp0 = ( ene_m*1/7 * NLP[atom][Rm][ist][2][L+1] * NLP[jatom][kl+1][jst][2][L+2]
                                    -ene_m*1/7 * NLP[atom][Rm][ist][2][L+2] * NLP[jatom][kl+1][jst][2][L+1]
                                    +ene_m*2/7 * NLP[atom][Rm][ist][2][L+3] * NLP[jatom][kl+1][jst][2][L+4]
                                    -ene_m*2/7 * NLP[atom][Rm][ist][2][L+4] * NLP[jatom][kl+1][jst][2][L+3]
                                    +ene_m*3/7 * NLP[atom][Rm][ist][2][L+5] * NLP[jatom][kl+1][jst][2][L+6]
                                    -ene_m*3/7 * NLP[atom][Rm][ist][2][L+6] * NLP[jatom][kl+1][jst][2][L+5])
                                    
                            Sum0_i += tmp0
                            Sum1_i += -tmp0
                        end

                        # diagonal contribution on up-up and dn-dn
                        for _ = 0:L2
                            Sum0_r += PFp*ene_p * NLP[atom][Rm][ist][1][L] * NLP[jatom][kl+1][jst][1][L]
                            Sum1_r += PFp*ene_p * NLP[atom][Rm][ist][1][L] * NLP[jatom][kl+1][jst][1][L]

                            Sum0_r += PFm*ene_m * NLP[atom][Rm][ist][2][L] * NLP[jatom][kl+1][jst][2][L]
                            Sum1_r += PFm*ene_m * NLP[atom][Rm][ist][2][L] * NLP[jatom][kl+1][jst][2][L]

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

        rcut = Atom_Cut1[atom] + Atom_Cut1[jatom]
        dmp = dampingF(rcut, Dis[atom][Rn])
        for ist = 1:NO0, jst = 1:NO1
            hst += 1
            HNL[1][Hks_Num+hst] = dmp * HNL_temp[ist,jst,1]
            HNL[2][Hks_Num+hst] = dmp * HNL_temp[ist,jst,2]
            HNL[3][Hks_Num+hst] = dmp * HNL_temp[ist,jst,3]
            iHNL[1][Hks_Num+hst] = dmp * iHNL_temp[ist,jst,1]
            iHNL[2][Hks_Num+hst] = dmp * iHNL_temp[ist,jst,2]
            iHNL[3][Hks_Num+hst] = dmp * iHNL_temp[ist,jst,3]
        end
    end


    MPI.Allreduce!(HNL[1], MPI.SUM, comm)
    MPI.Allreduce!(HNL[2], MPI.SUM, comm)
    MPI.Allreduce!(HNL[3], MPI.SUM, comm)
    MPI.Allreduce!(iHNL[1], MPI.SUM, comm)
    MPI.Allreduce!(iHNL[2], MPI.SUM, comm)
    MPI.Allreduce!(iHNL[3], MPI.SUM, comm)
end
