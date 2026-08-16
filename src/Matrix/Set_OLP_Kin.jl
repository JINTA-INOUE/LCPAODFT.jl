@timeit timer "Set_OLP_Kin" function Set_OLP_Kin!(OLP, MPI_Hkin, pao::Vector{PAO}, system_grid::System_Grid)
    
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
    dk = (Nkmax-Radial_kmin)/NkGrid
    for ik = 1:NkGrid+1
        k1[ik] = Radial_kmin + (ik-1)*dk
    end
    @. k2 = k1^2
    k2[begin] = 0.5*k2[begin]
    k2[end] = 0.5*k2[end]


    atv = system_grid.atv
	Gxyz = system_grid.Gxyz
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    MPI_size = system_grid.MPI_size
    MPHks = system_grid.MPHks
    Hks_Num = MPHks[myrank+1]
    



 
    fsize = maximum(Total_NumOrbs)
    OLPiαjβ = zeros(ComplexF64, fsize, fsize)
    Hkiniαjβ = zeros(ComplexF64, fsize, fsize)
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

    SphB = zeros(Float64, NkGrid+1)
    SphB2 = zeros(Float64, NkGrid+1)
    SphB4 = zeros(Float64, NkGrid+1)
    SumS0 = zeros(Float64, 4, 4, 4, 4)
    SumK0 = zeros(Float64, 4, 4, 4, 4)
    tmpL = zeros(Float64, NkGrid+1)
    tmpH = zeros(ComplexF64, fsize)


    hst = 0
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        ispe = atom2spe[atom]
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        cell = MPI_ncn[loop]+1
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        x = Gxyz[jatom][1] + atv[cell][1] - Gxyz[atom][1]
        y = Gxyz[jatom][2] + atv[cell][2] - Gxyz[atom][2]
        z = Gxyz[jatom][3] + atv[cell][3] - Gxyz[atom][3]
        
        R = sqrt(x^2 + y^2 + z^2)
        R = ifelse(R < 1e-10, 1e-10, R)

        iMaxL_Basis = pao[ispe].Spe_MaxL_Basis
        iNum_Basis = pao[ispe].Spe_Num_Basis
        iRF_Bessel = pao[ispe].Spe_RF_Bessel

        jMaxL_Basis = pao[jspe].Spe_MaxL_Basis
        jNum_Basis = pao[jspe].Spe_Num_Basis
        jRF_Bessel = pao[jspe].Spe_RF_Bessel
            
        Lmax_Four_Int = 2*max(iMaxL_Basis, jMaxL_Basis)

        fill!(OLPiαjβ, 0.0)
        fill!(Hkiniαjβ, 0.0)

        # Σ_{L=0}^{Lmax_Four_Int}Sum_{M=-L}^{L}
        for L = 0:Lmax_Four_Int
                
            @. SphB = SphericalBesselj(L, R*k1)
            @. SphB2 = SphB*k2
            @. SphB4 = SphB2*k2

            for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], ll = 0:jMaxL_Basis, pp = 1:jNum_Basis[ll+1]
                @. tmpL = iRF_Bessel[l+1][p]*jRF_Bessel[ll+1][pp]
                SumS0[pp,ll+1,p,l+1] = dot(SphB2, tmpL)*dk
                SumK0[pp,ll+1,p,l+1] = dot(SphB4, tmpL)*dk
            end

            for M = -L:L
                ist = 0
                for l = 0:iMaxL_Basis, p = 1:iNum_Basis[l+1], m = -l:l
                    jst = 0
                    ist += 1
                    @inbounds for ll = 0:jMaxL_Basis, pp = 1:jNum_Basis[ll+1], mm = -ll:ll
                        jst += 1
                        if abs(ll-L) <= l <= abs(ll+L) && iszero(m-mm-M) && abs(m-M) <= ll
                            Ylm = Ylm_complex(L,M,x,y,z)
                            Ls = Float64(L+ll-l)
                            iYC = (-im)^Ls * conj(Ylm) * Gaunt(f,l,m,ll,mm,L,M)
                            OLPiαjβ[ist,jst] += iYC*SumS0[pp,ll+1,p,l+1]
                            Hkiniαjβ[ist,jst] += iYC*SumK0[pp,ll+1,p,l+1]
                        end
                    end
                end
            end
        end


        # complex to real
        @inbounds for ist = 1:NO0
            @views mul!(tmpH, Cjβ[jspe], OLPiαjβ[ist,:])
            @views OLPiαjβ[ist,:] = tmpH
            @views mul!(tmpH, Cjβ[jspe], Hkiniαjβ[ist,:])
            @views Hkiniαjβ[ist,:] = tmpH
        end

        @inbounds for jst = 1:NO1
            @views mul!(tmpH, Ciα[ispe], OLPiαjβ[:,jst])
            @views OLPiαjβ[:,jst] = tmpH
            @views mul!(tmpH, Ciα[ispe], Hkiniαjβ[:,jst])
            @views Hkiniαjβ[:,jst] = tmpH
        end

        
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            OLP[Hks_Num+hst] = 8*real(OLPiαjβ[ist,jst])
            MPI_Hkin[hst] = 4*real(Hkiniαjβ[ist,jst])
        end
    end

    MPI.Allreduce!(OLP, MPI.SUM, comm)
end