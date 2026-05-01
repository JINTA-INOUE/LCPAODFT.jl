function Set_OLP_Kinforce(pao::Vector{PAO}, system_grid::System_Grid)

    Natom = system_grid.Natom
	FNAN = system_grid.FNAN
	natn = system_grid.natn
	Total_NumOrbs = system_grid.Total_NumOrbs
    Total_Hsize = system_grid.Total_Hsize


    Hkin_force = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 3)
    for xyz = 1:3
        Hkin_force[xyz] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
        for atom = 1:Natom
            Hkin_force[xyz][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                Hkin_force[xyz][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    Hkin_force[xyz][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end

    OLP_force = Vector{Vector{Float64}}(undef, 3)
    for xyz = 1:3
        OLP_force[xyz] = zeros(Float64, Total_Hsize)
    end

    Set_OLP_Kinforce!(OLP_force, Hkin_force, pao, system_grid)


    return OLP_force, Hkin_force
end


function Set_OLP_Kinforce!(OLP_force, Hkin_force, pao::Vector{PAO}, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    atom2spe = system_grid.atom2spe
    Nspecies = length(pao)

    Lmax = 0
    for spe = 1:Nspecies
        Spe_MaxL_Basis = pao[spe].Spe_MaxL_Basis
        Lmax = max(Spe_MaxL_Basis, Lmax)
    end
    Lmax = Lmax+1

    k1 = zeros(Float64, NkGrid+1)
    k2 = zeros(Float64, NkGrid+1)
    k3 = zeros(Float64, NkGrid+1)
    k4 = zeros(Float64, NkGrid+1)
    k5 = zeros(Float64, NkGrid+1)
    dk = (Nkmax-Radial_kmin)/NkGrid
    for ik = 1:NkGrid+1
        k1[ik] = Radial_kmin + (ik-1)*dk
    end
    @. k2 = k1^2
    @. k3 = k1^3
    @. k4 = k1^4
    @. k5 = k1^5
    k2[begin] = 0.5*k2[begin]
    k2[end] = 0.5*k2[end]
    k3[begin] = 0.5*k3[begin]
    k3[end] = 0.5*k3[end]
    k4[begin] = 0.5*k4[begin]
    k4[end] = 0.5*k4[end]
    k5[begin] = 0.5*k5[begin]
    k5[end] = 0.5*k5[end]
    

    Natom = system_grid.Natom
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    Total_Hsize = system_grid.Total_Hsize
    MPHks = system_grid.MPHks
    MPI_Hsize = system_grid.MPI_Hsize
    MPI_size = system_grid.MPI_size
    HksNum = MPHks[myrank+1]
    myHsize = MPI_Hsize[myrank+1]

    natn = system_grid.natn
    FNAN = system_grid.FNAN
    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    Total_NumOrbs = system_grid.Total_NumOrbs


    MPI_OLPx = zeros(Float64, myHsize)
    MPI_OLPy = zeros(Float64, myHsize)
    MPI_OLPz = zeros(Float64, myHsize)
    MPI_Hkinx = zeros(Float64, myHsize)
    MPI_Hkiny = zeros(Float64, myHsize)
    MPI_Hkinz = zeros(Float64, myHsize)


 
    fsize = maximum(Total_NumOrbs)
    OLPriαjβ = zeros(ComplexF64, fsize, fsize)
    OLPtiαjβ = zeros(ComplexF64, fsize, fsize)
    OLPpiαjβ = zeros(ComplexF64, fsize, fsize)
    Hkinriαjβ = zeros(ComplexF64, fsize, fsize)
    Hkintiαjβ = zeros(ComplexF64, fsize, fsize)
    Hkinpiαjβ = zeros(ComplexF64, fsize, fsize)
    Ciα = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    Cjβ = Vector{Matrix{ComplexF64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Ciα[spe] = zeros(ComplexF64, fsize, fsize)
        Set_Comp2Real!( Ciα[spe], pao[spe].Spe_MaxL_Basis, pao[spe].Spe_Num_Basis )
        Cjβ[spe] = deepcopy(Ciα[spe])
        conj!(Ciα[spe])
    end


    Lmax2 = Lmax*2
    fact2 = zeros(Float64, 2*Lmax2+1, 2*Lmax2+1)
    for i = 0:2*Lmax2, j = 0:2*Lmax2
        tmp0 = sqrt(factorial(big(i)))
        tmp1 = sqrt(factorial(big(j)))
        fact2[i+1,j+1] = tmp0/tmp1
    end


    f = zeros(Float64, S3J_MAX_FACT)
    _Set_f_for_Gaunt!(f)

    asize_lmax = 30
    tsb = zeros(Float64, asize_lmax+10)
    SphB_l = zeros(Float64, 2*Lmax)
    dSphB_l = zeros(Float64, 2*Lmax)
    SphB = Vector{Vector{Float64}}(undef, 2*Lmax)
    dSphB = Vector{Vector{Float64}}(undef, 2*Lmax)
    for l = 1:2*Lmax
        SphB[l] = zeros(Float64, NkGrid+1)
        dSphB[l] = zeros(Float64, NkGrid+1)
    end
    SphB2 = zeros(Float64, NkGrid+1)
    dSphB3 = zeros(Float64, NkGrid+1)
    SphB4 = zeros(Float64, NkGrid+1)
    dSphB5 = zeros(Float64, NkGrid+1)
    SumS0 = zeros(Float64, 4, 4, 4, 4)
    SumSr0 = zeros(Float64, 4, 4, 4, 4)
    SumK0 = zeros(Float64, 4, 4, 4, 4)
    SumKr0 = zeros(Float64, 4, 4, 4, 4)
    tmpL = zeros(Float64, NkGrid+1)
    tmpH = zeros(ComplexF64, fsize)


    SH = zeros(Float64, 2)
    dSHt = zeros(Float64, 2)
    dSHp = zeros(Float64, 2)



    MPI.Barrier(comm)
    hst = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
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

        ispe = atom2spe[atom]
        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        iRF_Bessel = pao[ispe].Spe_RF_Bessel

        jMaxL_Basis = pao[jspe].Spe_MaxL_Basis
        jNum_Basis = pao[jspe].Spe_Num_Basis
        jRF_Bessel = pao[jspe].Spe_RF_Bessel

        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
            
        Lmax_Four_Int = 2*max(iMaxL_Basis, jMaxL_Basis)

        for ik = 1:NkGrid+1
            Calc_SphericalBesselj!(Lmax_Four_Int, R*k1[ik], tsb, SphB_l, dSphB_l)
            for l = 1:Lmax_Four_Int+1
                SphB[l][ik] = SphB_l[l]
                dSphB[l][ik] = dSphB_l[l]
            end
        end

        
        fill!(OLPriαjβ, 0.0)
        fill!(OLPtiαjβ, 0.0)
        fill!(OLPpiαjβ, 0.0)
        fill!(Hkinriαjβ, 0.0)
        fill!(Hkintiαjβ, 0.0)
        fill!(Hkinpiαjβ, 0.0)


        # Σ_{L=0}^{Lmax_Four_Int}Sum_{M=-L}^{L}
        for L = 0:Lmax_Four_Int
                
            @. SphB2 = SphB[L+1]*k2
            @. SphB4 = SphB[L+1]*k4
            @. dSphB3 = dSphB[L+1]*k3
            @. dSphB5 = dSphB[L+1]*k5

            for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], ll = 0:jMaxL_Basis, pp = 1:jNum_Basis[ll+1]
                @. tmpL = iRF_Bessel[l+1][p]*jRF_Bessel[ll+1][pp]
                SumS0[pp,ll+1,p,l+1] = dot(SphB2, tmpL)*dk
                SumSr0[pp,ll+1,p,l+1] = dot(dSphB3, tmpL)*dk
                SumK0[pp,ll+1,p,l+1] = dot(SphB4, tmpL)*dk
                SumKr0[pp,ll+1,p,l+1] = dot(dSphB5, tmpL)*dk
            end


            for M = -L:L
                ist = 0
                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                    jst = 0
                    ist += 1
                    for ll = 0:jMaxL_Basis, pp = 1:jNum_Basis[ll+1], mm = -ll:ll
                        jst += 1
                        Ls = Float64(L+ll-l)
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
                                
                            OLPriαjβ[ist,jst] += iYC*SumSr0[pp,ll+1,p,l+1]
                            OLPtiαjβ[ist,jst] += iYCt*SumS0[pp,ll+1,p,l+1]
                            OLPpiαjβ[ist,jst] += iYCp*SumS0[pp,ll+1,p,l+1]
                            Hkinriαjβ[ist,jst] += iYC*SumKr0[pp,ll+1,p,l+1]
                            Hkintiαjβ[ist,jst] += iYCt*SumK0[pp,ll+1,p,l+1]
                            Hkinpiαjβ[ist,jst] += iYCp*SumK0[pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end


        # complex to real
        for ist = 1:NO0
            @views mul!(tmpH, Cjβ[jspe], OLPriαjβ[ist,:])
            @views OLPriαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], OLPtiαjβ[ist,:])
            @views OLPtiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], OLPpiαjβ[ist,:])
            @views OLPpiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], Hkinriαjβ[ist,:])
            @views Hkinriαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], Hkintiαjβ[ist,:])
            @views Hkintiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], Hkinpiαjβ[ist,:])
            @views Hkinpiαjβ[ist,:] = tmpH
        end

        for jst = 1:NO1
            @views mul!(tmpH, Ciα[ispe], OLPriαjβ[:,jst])
            @views OLPriαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], OLPtiαjβ[:,jst])
            @views OLPtiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], OLPpiαjβ[:,jst])
            @views OLPpiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], Hkinriαjβ[:,jst])
            @views Hkinriαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], Hkintiαjβ[:,jst])
            @views Hkintiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], Hkinpiαjβ[:,jst])
            @views Hkinpiαjβ[:,jst] = tmpH
        end

    
        if Rn ≠ 1
            if abs(siT) < 1.0e-13
                for ist = 1:NO0, jst = 1:NO1
                    hst += 1
                    MPI_OLPx[hst] = -8*real(siT*coP*OLPriαjβ[ist,jst] + coT*coP/R*OLPtiαjβ[ist,jst])
                    MPI_OLPy[hst] = -8*real(siT*siP*OLPriαjβ[ist,jst] + coT*siP/R*OLPtiαjβ[ist,jst])
                    MPI_OLPz[hst] = -8*real(coT*OLPriαjβ[ist,jst] - siT/R*OLPtiαjβ[ist,jst])
                    MPI_Hkinx[hst] = -4*real(siT*coP*Hkinriαjβ[ist,jst] + coT*coP/R*Hkintiαjβ[ist,jst])
                    MPI_Hkiny[hst] = -4*real(siT*siP*Hkinriαjβ[ist,jst] + coT*siP/R*Hkintiαjβ[ist,jst])
                    MPI_Hkinz[hst] = -4*real(coT*Hkinriαjβ[ist,jst] - siT/R*Hkintiαjβ[ist,jst])
                end
            else
                for ist = 1:NO0, jst = 1:NO1
                    hst += 1
                    MPI_OLPx[hst] = -8*real(siT*coP*OLPriαjβ[ist,jst] + coT*coP/R*OLPtiαjβ[ist,jst] - siP/siT/R*OLPpiαjβ[ist,jst])
                    MPI_OLPy[hst] = -8*real(siT*siP*OLPriαjβ[ist,jst] + coT*siP/R*OLPtiαjβ[ist,jst] + coP/siT/R*OLPpiαjβ[ist,jst])
                    MPI_OLPz[hst] = -8*real(coT*OLPriαjβ[ist,jst] - siT/R*OLPtiαjβ[ist,jst])
                    MPI_Hkinx[hst] = -4*real(siT*coP*Hkinriαjβ[ist,jst] + coT*coP/R*Hkintiαjβ[ist,jst] - siP/siT/R*Hkinpiαjβ[ist,jst])
                    MPI_Hkiny[hst] = -4*real(siT*siP*Hkinriαjβ[ist,jst] + coT*siP/R*Hkintiαjβ[ist,jst] + coP/siT/R*Hkinpiαjβ[ist,jst])
                    MPI_Hkinz[hst] = -4*real(coT*Hkinriαjβ[ist,jst] - siT/R*Hkintiαjβ[ist,jst])
                end
            end
        else
            for ist = 1:NO0, jst = 1:NO1
                hst += 1
                MPI_OLPx[hst] = 0.0
                MPI_OLPy[hst] = 0.0
                MPI_OLPz[hst] = 0.0
                MPI_Hkinx[hst] = 0.0
                MPI_Hkiny[hst] = 0.0
                MPI_Hkinz[hst] = 0.0
            end
        end
    end
    MPI.Barrier(comm)

    

    Hkinforce1Dx = zeros(Float64, Total_Hsize)
    Hkinforce1Dy = zeros(Float64, Total_Hsize)
    Hkinforce1Dz = zeros(Float64, Total_Hsize)
    
    MPI.Allgatherv!(MPI_OLPx, VBuffer(OLP_force[1], MPI_Hsize), comm)
    MPI.Allgatherv!(MPI_OLPy, VBuffer(OLP_force[2], MPI_Hsize), comm)
    MPI.Allgatherv!(MPI_OLPz, VBuffer(OLP_force[3], MPI_Hsize), comm)
    MPI.Allgatherv!(MPI_Hkinx, VBuffer(Hkinforce1Dx, MPI_Hsize), comm)
    MPI.Allgatherv!(MPI_Hkiny, VBuffer(Hkinforce1Dy, MPI_Hsize), comm)
    MPI.Allgatherv!(MPI_Hkinz, VBuffer(Hkinforce1Dz, MPI_Hsize), comm)



    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        counts += 1
        Hkin_force[1][atom][Rn][ist][jst] = Hkinforce1Dx[counts]
        Hkin_force[2][atom][Rn][ist][jst] = Hkinforce1Dy[counts]
        Hkin_force[3][atom][Rn][ist][jst] = Hkinforce1Dz[counts]
    end
    MPI.Barrier(comm)
end
