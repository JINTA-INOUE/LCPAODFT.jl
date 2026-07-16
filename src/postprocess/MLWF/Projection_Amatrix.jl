function Projection_Amatrix(kg, Nk, Cnk, mlwf_setup::MLWF_Setup)

    MLWF_kmesh = mlwf_setup.MLWF_kmesh
    Nkpt = prod(MLWF_kmesh)
    WANNUM = mlwf_setup.MLWF_Func_Num
    MLWF_Num_Kinds_Projectors = mlwf_setup.MLWF_Num_Kinds_Projectors
    MLWF_ProName = mlwf_setup.MLWF_ProName
    MLWF_Z_Direction = mlwf_setup.MLWF_Z_Direction
    MLWF_X_Direction = mlwf_setup.MLWF_X_Direction
    MLWF_Cut1 = mlwf_setup.MLWF_Cut1
    MLWF_ProjOrbitals = mlwf_setup.MLWF_ProjOrbitals
    MLWF_Gxyz_AU = mlwf_setup.MLWF_Gxyz_AU

    material = mlwf_setup.material
    SpinPol = material.SpinPol
    SpinP_switch = material.SpinP_switch
    TCpyCell = material.TCpyCell
    Natom = material.Natom
    Nspecies = material.Nspecies
    atv = material.atv
    atv_ijk = material.atv_ijk
    Gxyz = material.Gxyz
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    AtomOrbital = material.AtomOrbital
    xc_type = material.xc_type

    if SpinPol in ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end



    pao, _ = PAO_Pspot(AtomOrbital, xc_type, false) 
    proj_pao, _ = PAO_Pspot(MLWF_ProjOrbitals, xc_type, false) 


    FNAN_WP, natn_WP, ncn_WP = Find_NN_Projectors_Basis(MLWF_Cut1, MLWF_Gxyz_AU, TCpyCell, Natom, Gxyz, atv, [7.0,7.0])
    @show FNAN_WP
    OLP_WP = Set_OLP_WP(pao, proj_pao, FNAN_WP, natn_WP, ncn_WP, mlwf_setup)


    BANDNUM = 12

    
    Amnk = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Amnk[spin] = Vector{Vector{Vector{ComplexF64}}}(undef, Nkpt)
        for ik = 1:Nkpt
            Amnk[spin][ik] = Vector{Vector{ComplexF64}}(undef, BANDNUM)
            for μ = 1:BANDNUM
                Amnk[spin][ik][μ] = zeros(ComplexF64, WANNUM)
            end
        end
    end



    Projection_Amatrix!(Amnk, kg, Nk, Cnk, FNAN_WP, natn_WP, ncn_WP, OLP_WP, mlwf_setup)


    # C_Projection_Amatrix()
end


function Projection_Amatrix!(Amnk, kg, Nk, Cnk, FNAN_WP, natn_WP, ncn_WP, OLP_WP, mlwf_setup::MLWF_Setup)

    material = mlwf_setup.material
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)

    MLWF_Num_Kinds_Projectors = mlwf_setup.MLWF_Num_Kinds_Projectors
    MLWF_ProName = mlwf_setup.MLWF_ProName
    MLWF_Z_Direction = mlwf_setup.MLWF_Z_Direction
    MLWF_X_Direction = mlwf_setup.MLWF_X_Direction
    MLWF_Total_NumOrbs = mlwf_setup.MLWF_Total_NumOrbs


    MLWF_Num_Pro, MLWF_NumL_Pro, MLWF_Select_Matrix, MLWF_Projector_Hybridize_Matrix = Analyze_Wannier_Projectors(MLWF_Num_Kinds_Projectors, MLWF_ProName)
    MLWF_Euler_Rotation_Angle = Get_Wannier_Euler_Rotation_Angle(MLWF_Num_Kinds_Projectors, MLWF_Z_Direction, MLWF_X_Direction)

    MLWF_RotMat_for_Real_Func = Vector{Vector{Matrix{Float64}}}(undef, MLWF_Num_Kinds_Projectors)
    for p = 1:MLWF_Num_Kinds_Projectors
        MLWF_RotMat_for_Real_Func[p] = Vector{Matrix{Float64}}(undef, 3)
        for L = 1:3
            MLWF_RotMat_for_Real_Func[p][L] = zeros(Float64, 2*L+1, 2*L+1)
            if MLWF_NumL_Pro[p,L+1] ≠ 0
                Get_Rotational_Matrix!(L, MLWF_Euler_Rotation_Angle[p,:], MLWF_RotMat_for_Real_Func[p][L])
            end
        end
    end

    psize = 0
    for proj = 1:MLWF_Num_Kinds_Projectors, L = 0:3
        psize += MLWF_NumL_Pro[proj,L+1]*(2*L+1)
    end

    MLWF_MP = zeros(Int32, MLWF_Num_Kinds_Projectors+1)
    Sum = 0
    for atom = 1:MLWF_Num_Kinds_Projectors
        MLWF_MP[atom+1] = Sum + MLWF_Total_NumOrbs[atom]
        Sum += MLWF_Total_NumOrbs[atom]
    end


    spinsize = length(Amnk)
    Nkpt = length(Amnk[1])
    BANDNUM = length(Amnk[1][1])
    WANNUM = length(Amnk[1][1][1])

    
    disentangle = ifelse(BANDNUM>WANNUM, true, false)

    tmpAmnk = zeros(ComplexF64, fsize, psize)
    tmpResults = zeros(ComplexF64, max(WANNUM,psize,7))

    for spin = 1:spinsize, ik = 1:Nkpt
        band_num = ifelse(disentangle, Nk[spin,ik,2]-Nk[spin,ik,1], WANNUM)
        band_shift = Nk[spin,ik,1]
        for μ = 1:band_num
            windx = 1
            for proj = 1:MLWF_Num_Kinds_Projectors, ist = 1:MLWF_Total_NumOrbs[proj]
                Sum = 0.0 + im*0.0
                for Rn = 1:FNAN_WP[proj]
                    jatom = natn_WP[proj][Rn]
                    cell = ncn_WP[proj][Rn]+1
                    Bnum = MLWF_MP[jatom]
                    NO1 = MLWF_Total_NumOrbs[jatom]
                    
                    kRn = -2*pi*dot(kg[ik], atv_ijk[cell])
                    ex = cispi(im*kRn)
                    @inbounds for jst = 1:NO1
                        Sum += Cnk[spin][ik][jst+Bnum,μ+band_shift]*OLP_WP[proj][Rn][ist][jst]
                    end

                    Sum = Sum*ex
                end

                tmpAmnk[μ,windx] = Sum/sqrt(abs(FNAN_WP[proj]))
                windx += 1
            end
        end


        for μ = 1:band_num
            windx = 0
            for proj = 1:MLWF_Num_Kinds_Projectors
                for L = 0:3
                    if MLWF_NumL_Pro[proj,L+1]≠0 && L≠0
                        for i = 1:2*L+1
                            Sum = 0.0 + im*0.0
                            @inbounds for j = 1:2*L+1
                                Sum += MLWF_RotMat_for_Real_Func[proj][L+1][i,j]*tmpAmnk[μ,windx+j]
                            end

                            tmpResults[i] = Sum
                        end

                        for i = 1:2*L+1
                            tmpAmnk[μ,i] = tmpResults[i]
                        end

                        windx += 2*L+1
                    elseif MLWF_NumL_Pro[proj,L+1]≠0 && L==0
                        windx += 1
                    end
                end
            end
        end

        
        for μ = 1:band_num
            nindx = 0
            for proj = 1:MLWF_Num_Kinds_Projectors
                for i = 1:MLWF_Num_Pro[proj]
                    Sum = 0.0 + im*0.0
                    @inbounds for j = 1:MLWF_Num_Pro[proj]
                        Sum += MLWF_Projector_Hybridize_Matrix[proj][i,j]*Amnk[spin][ik][μ][nindx+j]
                    end

                    tmpResults[i] = Sum
                end


                for i = 1:MLWF_Num_Pro[proj]
                    Amnk[spin][ik][μ][i] = tmpResults[i]
                end

                nindx += Wannier_Num_Pro[proj]
            end
        end
    end 
end