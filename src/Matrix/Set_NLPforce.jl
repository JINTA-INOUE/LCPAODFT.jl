"""
    Nonlocal potentials matrix
"""
function Set_NLPforce(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs


    VPS_j_Num = zeros(Int64, Natom)
    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot[spe].Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end

        VPS_j_Num[atom] = pspot[spe].VPS_j_dependency
        NLTotal_Num[atom] = tot
    end



    maxFNAN = maximum(FNAN)
	maxTotal_NumOrbs = maximum(Total_NumOrbs)
    maxVPS_j_Num = maximum(VPS_j_Num)
	maxNLTotal_Num = maximum(NLTotal_Num)+2
    NLPforce = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, 4)
    for xyz = 1:4
        NLPforce[xyz] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Natom+1)
        for atom = 1:Natom+1
            if atom == Natom+1
				fan = maxFNAN+1
				NO0 = maxTotal_NumOrbs
			else
				fan = FNAN[atom]+1
				NO0 = Total_NumOrbs[atom]
			end
            
            NLPforce[xyz][atom] = Vector{Vector{Vector{Vector{Float64}}}}(undef, fan)
            for Rn = 1:fan
                if atom == Natom+1
                    VPS_j_dependency = maxVPS_j_Num
                    NO1 = maxNLTotal_Num
                else
                    jatom = natn[atom][Rn]
                    VPS_j_dependency = VPS_j_Num[jatom]
                    NO1 = NLTotal_Num[jatom]+2
                end
                
                NLPforce[xyz][atom][Rn] = Vector{Vector{Vector{Float64}}}(undef, NO0)
                for ist = 1:NO0
                    NLPforce[xyz][atom][Rn][ist] = Vector{Vector{Float64}}(undef, VPS_j_dependency+1)
                    for so = 1:VPS_j_dependency+1
                        NLPforce[xyz][atom][Rn][ist][so] = zeros(Float64, NO1)
                    end
                end
            end
        end
    end
    Set_NLPforce!(NLPforce, pao, pspot, system_grid)


    return NLPforce
end


function Set_NLPforce!(NLPforce, pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

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
    end
    for spe = 1:Nspecies
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
    myNloop = system_grid.MPI_size
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz


    myNLPsize = 0
    for loop = 1:myNloop
        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        VPS_j_dependency = pspot[jspe].VPS_j_dependency
        for so = 1:VPS_j_dependency+1, ist = 1:Total_NumOrbs[atom], jst = 1:NLTotal_Num[jatom]
            myNLPsize += 1
        end
    end

    MPI_NLP1D = zeros(Float64, myNLPsize)
    MPI_NLP1Dx = zeros(Float64, myNLPsize)
    MPI_NLP1Dy = zeros(Float64, myNLPsize)
    MPI_NLP1Dz = zeros(Float64, myNLPsize)




    fsize = maximum(Total_NumOrbs)
    NLfsize = maximum(NLTotal_Num)
    NLPiαjβ = zeros(ComplexF64, fsize, NLfsize)
    NLPriαjβ = zeros(ComplexF64, fsize, NLfsize)
    NLPtiαjβ = zeros(ComplexF64, fsize, NLfsize)
    NLPpiαjβ = zeros(ComplexF64, fsize, NLfsize)
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
        SphB[l] = zeros(Float64, NkGrid+1)
        dSphB[l] = zeros(Float64, NkGrid+1)
    end
    SphB2 = zeros(Float64, NkGrid+1)
    dSphB3 = zeros(Float64, NkGrid+1)
    SumNL0 = zeros(Float64, 15, 4, 4)
    SumNLr0 = zeros(Float64, 15, 4, 4)
    tmpL = zeros(Float64, NkGrid+1)
    tmpH1 = zeros(ComplexF64, fsize)
    tmpH2 = zeros(ComplexF64, NLfsize)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)




    MPI.Barrier(comm)
    counts = 0
    for loop = 1:myNloop

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
        siT = sin(theta)
        coT = cos(theta)
        siP = sin(phi)
        coP = cos(phi)
            

        jNum_RVPS = pspot[jspe].Spe_Num_RVPS
        jVPS_List = pspot[jspe].Spe_VPS_List
        VPS_j_dependency = pspot[jspe].VPS_j_dependency

        Lmax = maximum(jVPS_List)
        Lmax_Four_Int = 2*ifelse(Lmax>iMaxL_Basis, Lmax, iMaxL_Basis)

        for ik = 1:NkGrid+1
            Calc_SphericalBesselj!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                dSphB[l][ik] = dSphB_l[l]
            end
        end
        

        for so = 1:VPS_j_dependency+1

            fill!(NLPiαjβ, 0.0)
            fill!(NLPriαjβ, 0.0)
            fill!(NLPtiαjβ, 0.0)
            fill!(NLPpiαjβ, 0.0)

            # Σ_{L=0}^{Lmax_Four_Int}Sum_{M=-L}^{L}
            for L = 0:Lmax_Four_Int

                @. SphB2 = SphB[L+1]*k2
                @. dSphB3 = dSphB[L+1]*k3

                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], lnum = 1:jNum_RVPS
                    @. tmpL = iRF_Bessel[l+1][p]*NLRF_Bessel[jspe][so][lnum]
                    SumNL0[lnum,p,l+1] = dot(SphB2, tmpL)*dk
                    SumNLr0[lnum,p,l+1] = dot(dSphB3, tmpL)*dk
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

                                NLPiαjβ[ist,jst] += iYC*SumNL0[lnum,p,l+1]
                                NLPriαjβ[ist,jst] += iYC*SumNLr0[lnum,p,l+1]
                                NLPtiαjβ[ist,jst] += iYCt*SumNL0[lnum,p,l+1]
                                NLPpiαjβ[ist,jst] += iYCp*SumNL0[lnum,p,l+1]
                            end
                        end
                    end
                end
            end
            

            # complex to real        
            for ist = 1:NO0
                @views mul!(tmpH2, Cjβ[jspe], NLPiαjβ[ist,:])
                @views NLPiαjβ[ist,:] = tmpH2
                @views mul!(tmpH2, Cjβ[jspe], NLPriαjβ[ist,:])
                @views NLPriαjβ[ist,:] = tmpH2
                @views mul!(tmpH2, Cjβ[jspe], NLPtiαjβ[ist,:])
                @views NLPtiαjβ[ist,:] = tmpH2
                @views mul!(tmpH2, Cjβ[jspe], NLPpiαjβ[ist,:])
                @views NLPpiαjβ[ist,:] = tmpH2
            end

            for jst = 1:NO1
                @views mul!(tmpH1, Ciα[ispe], NLPiαjβ[:,jst])
                @views NLPiαjβ[:,jst] = tmpH1
                @views mul!(tmpH1, Ciα[ispe], NLPriαjβ[:,jst])
                @views NLPriαjβ[:,jst] = tmpH1
                @views mul!(tmpH1, Ciα[ispe], NLPtiαjβ[:,jst])
                @views NLPtiαjβ[:,jst] = tmpH1
                @views mul!(tmpH1, Ciα[ispe], NLPpiαjβ[:,jst])
                @views NLPpiαjβ[:,jst] = tmpH1
            end
            


            if Rn ≠ 1
                if abs(siT) < 1.0e-13
                    for ist = 1:NO0, jst = 1:NO1
                        counts += 1
                        MPI_NLP1D[counts] = 8*real(NLPiαjβ[ist,jst])
                        MPI_NLP1Dx[counts] = -8*real(siT*coP*NLPriαjβ[ist,jst] + coT*coP/R*NLPtiαjβ[ist,jst])
                        MPI_NLP1Dy[counts] = -8*real(siT*siP*NLPriαjβ[ist,jst] + coT*siP/R*NLPtiαjβ[ist,jst])
                        MPI_NLP1Dz[counts] = -8*real(coT*NLPriαjβ[ist,jst] - siT/R*NLPtiαjβ[ist,jst])
                    end
                else
                    for ist = 1:NO0, jst = 1:NO1
                        counts += 1
                        MPI_NLP1D[counts] = 8*real(NLPiαjβ[ist,jst])
                        MPI_NLP1Dx[counts] = -8*real(siT*coP*NLPriαjβ[ist,jst] + coT*coP/R*NLPtiαjβ[ist,jst] - siP/siT/R*NLPpiαjβ[ist,jst])
                        MPI_NLP1Dy[counts] = -8*real(siT*siP*NLPriαjβ[ist,jst] + coT*siP/R*NLPtiαjβ[ist,jst] + coP/siT/R*NLPpiαjβ[ist,jst])
                        MPI_NLP1Dz[counts] = -8*real(coT*NLPriαjβ[ist,jst] - siT/R*NLPtiαjβ[ist,jst])
                    end
                end
            else
                for ist = 1:NO0, jst = 1:NO1
                    counts += 1
                    MPI_NLP1D[counts] = 8*real(NLPiαjβ[ist,jst])
                    MPI_NLP1Dx[counts] = 0.0
                    MPI_NLP1Dy[counts] = 0.0
                    MPI_NLP1Dz[counts] = 0.0
                end
            end
        end
    end
    MPI.Barrier(comm)
    



    Total_NLPsize = MPI.Allreduce(myNLPsize, MPI.SUM, comm) 

    _counts = zeros(Int64, nprocs)
    for id = 1:nprocs
        if id-1 == myrank
            _counts[id] = myNLPsize
        end
        MPI.Barrier(comm)
    end
    MPI.Barrier(comm)

    MPI_NLPsize = zeros(Int64, nprocs)
    MPI.Allreduce!(_counts, MPI_NLPsize, nprocs, MPI.SUM, comm)  


    MPI_NLP = zeros(Float64, Total_NLPsize)
    MPI_NLPx = zeros(Float64, Total_NLPsize)
    MPI_NLPy = zeros(Float64, Total_NLPsize)
    MPI_NLPz = zeros(Float64, Total_NLPsize)

    MPI.Allgatherv!(MPI_NLP1D, VBuffer(MPI_NLP, MPI_NLPsize), comm)
    MPI.Allgatherv!(MPI_NLP1Dx, VBuffer(MPI_NLPx, MPI_NLPsize), comm)
    MPI.Allgatherv!(MPI_NLP1Dy, VBuffer(MPI_NLPy, MPI_NLPsize), comm)
    MPI.Allgatherv!(MPI_NLP1Dz, VBuffer(MPI_NLPz, MPI_NLPsize), comm)

    
    
    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        jspe = atom2spe[jatom]
        VPS_j_dependency = pspot[jspe].VPS_j_dependency
        for so = 1:VPS_j_dependency+1, ist = 1:Total_NumOrbs[atom], jst = 1:NLTotal_Num[jatom]
            counts += 1
            NLPforce[1][atom][Rn][ist][so][jst] = MPI_NLP[counts]
            NLPforce[2][atom][Rn][ist][so][jst] = MPI_NLPx[counts]
            NLPforce[3][atom][Rn][ist][so][jst] = MPI_NLPy[counts]
            NLPforce[4][atom][Rn][ist][so][jst] = MPI_NLPz[counts]
        end
    end
    MPI.Barrier(comm)

end
