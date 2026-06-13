using LinearAlgebra
using Printf

struct DFTD3vdW_params
    DFTD3_damp_dftD::Int64
    k1_dftD::Float64
    k2_dftD::Float64
    k3_dftD::Float64
    sr8_dftD::Float64
    alp6_dftD::Float64
    alp8_dftD::Float64
    cncut_dftD::Float64
    rcut_dftD::Float64
    s6_dftD::Float64
    s8_dftD::Float64
    sr6_dftD::Float64
    a1_dftD::Float64
    a2_dftD::Float64
end


function DFTD3vdW_params()

    k1_dftD = 16.0
    k2_dftD = 4.0/3.0
    k3_dftD = 4.0
    sr8_dftD = 1.0
    alp6_dftD = 14.0
    alp8_dftD = alp6_dftD+2.0
    DFTD3_damp_dftD = 2
    cncut_dftD = 40.0
    rcut_dftD = 100.0
    s6_dftD = 1.0
    s8_dftD = 0.7875
    sr6_dftD = 1.217
    a1_dftD = 0.4289
    a2_dftD = 4.4407

    return DFTD3vdW_params(
        k1_dftD, k2_dftD, k3_dftD, sr8_dft,
        alp6_dftD, alp8_dftD, 
        DFTD3_damp_dftD,
        cncut_dftD, rcut_dftD,
        s6_dftD, s8_dftD, sr6_dftD, a1_dftD, a2_dftD
    )
end


function limit_c6(iat::Int)
    iadr = 0
    at = iat

    while at > 100
        at -= 100
        iadr += 1
    end

    return iadr, at
end


function DFTD3vdW_init(Nspecies, Latvecs, Spe_Znumber::Vector{<:Integer}, DFTD_IntDir; verbosity=1)

    k1_dftD = 16.0
    k2_dftD = 4.0/3.0
    k3_dftD = 4.0
    sr8_dftD = 1.0
    alp6_dftD = 14.0
    alp8_dftD = alp6_dftD+2.0

    c6ab_tmp = DFTD3vdW_params_c6ab_tmp()
    r0ab_tmp = DFTD3vdW_params_r0ab_tmp()
    rcov, r2r4 = DFTD3vdW_params_rcov_r2r4()

    maxci = zeros(Int64, 95)
    c6ab = zeros(Float64, 95, 95, 5, 5, 3)
    k = 0
    for i = 1:nlines
        iat = floor(Int64, c6ab_tmp[k+2])
        jat = floor(Int64, c6ab_tmp[k+3])
        iadr, iat = limit_c6(iat)
        jadr, jat = limit_c6(jat)

        if iat<=94 && jat<=94
            maxci[iat] = max(maxci[iat], iadr)
            maxci[jat] = max(maxci[jat], jadr)
            c6ab[iat,jat,iadr+1,jadr+1,1] = c6ab_tmp[k+1]
            c6ab[iat,jat,iadr+1,jadr+1,2] = c6ab_tmp[k+4]
            c6ab[iat,jat,iadr+1,jadr+1,3] = c6ab_tmp[k+5]
            c6ab[jat,iat,jadr+1,iadr+1,1] = c6ab_tmp[k+1]
            c6ab[jat,iat,jadr+1,iadr+1,2] = c6ab_tmp[k+5]
            c6ab[jat,iat,jadr+1,iadr+1,3] = c6ab_tmp[k+4]
        end
        k = i*5
    end

    for iat = 1:95
        maxci[iat] += 1
    end


    k = 0
    r0ab = zeros(Float64, 95, 95)
    for i = 1:94, j = 1:i
        k += 1
        r0ab[i,j] = r0ab_tmp[k]/0.52917726
        r0ab[j,i] = r0ab_tmp[k]/0.52917726
    end
    
    
    
    for i = 1:94
        rcov[i] = k2_dftD*rcov[i]/0.52917726
    end

    Nspecies = 1
    Spe_Znumber = [6]
    DFTD_IntDir = [0, 0, 0]
    
    rcov_dftD = zeros(Float64, Nspecies)
    r2r4_dftD = zeros(Float64, Nspecies)
    maxcn_dftD = zeros(Int64, Nspecies)
    rcovab_dftD = zeros(Float64, Nspecies, Nspecies)
    r2r4ab_dftD = zeros(Float64, Nspecies, Nspecies)
    r0ab_dftD = zeros(Float64, Nspecies, Nspecies)

    C6ab_dftD = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspecies)
    for spe = 1:Nspecies
        C6ab_dftD[spe] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspecies)
        for jspe = 1:Nspecies
            C6ab_dftD[spe][jspe] = Vector{Vector{Vector{Vector{Float64}}}}(undef, 5)
            for i = 1:5
                C6ab_dftD[spe][jspe][i] = Vector{Vector{Float64}}(undef, 5)
                for j = 1:5
                    C6ab_dftD[spe][jspe][i][j] = zeros(Float64, 3)
                end
            end
        end
    end


    for spe = 1:Nspecies
        maxcn_dftD[spe] = -1
        rcov_dftD[spe] = -1.0
        r2r4_dftD[spe] = -1.0
    end

    for spe = 1:Nspecies
        iZ = Spe_Znumber[spe]
        if !(0 < iZ < 94)
            error("<DFTD3vdW_init> error in DFTD3vdw_init: DFTD3 only for atomic number supported 1 through 94.")
        end

        rcov_dftD[spe] = rcov[iZ]     
        r2r4_dftD[spe] = r2r4[iZ]
        maxcn_dftD[spe] = maxci[iZ]
    end


    for spe = 1:Nspecies, jspe = 1:Nspecies
        iZ = Spe_Znumber[spe]
        jZ = Spe_Znumber[jspe]

        r0ab_dftD[spe,jspe] = r0ab[iZ,jZ]
        r2r4ab_dftD[spe,jspe] = r2r4_dftD[spe]*r2r4_dftD[jspe]
        rcovab_dftD[spe,jspe] = rcov_dftD[spe]+rcov_dftD[jspe]

        for i = 1:maxci[iZ], j = 1:maxci[jZ]
            if c6ab[iZ,jZ,i,j,1] >= 0.0
                C6ab_dftD[spe][jspe][i][j][1] = c6ab[iZ,jZ,i,j,1]
                C6ab_dftD[spe][jspe][i][j][2] = c6ab[iZ,jZ,i,j,2]
                C6ab_dftD[spe][jspe][i][j][3] = c6ab[iZ,jZ,i,j,3]
            end
        end
    end

    verbosity = 1
    DFTD3_damp_dftD = 2
    cncut_dftD = 40.0
    rcut_dftD = 100.0
    s6_dftD = 1.0
    s8_dftD = 0.7875
    sr6_dftD = 1.217
    a1_dftD = 0.4289
    a2_dftD = 4.4407
    if verbosity >= 1
        println("<DFTD3vdW_init> DFT-D3 parameter")
        @printf("<DFTD3vdW_init> damping=%d (1 for ZERO, 2 for BJ)\n", DFTD3_damp_dftD)
        @printf("<DFTD3vdW_init> CN_TH=%f, r_cut=%f\n", cncut_dftD, rcut_dftD)
        @printf("<DFTD3vdW_init> k1=%f, k2=%f, k3=%f\n", k1_dftD, k2_dftD, k3_dftD)
        @printf("<DFTD3vdW_init> s6=%f, s8=%f\n", s6_dftD, s8_dftD)
        if DFTD3_damp_dftD == 1
            @printf("<DFTD3vdW_init> alp6=%f, alp8=%f, sr6=%f\n", alp6_dftD, alp8_dftD, sr6_dftD)
        elseif DFTD3_damp_dftD == 2
            @printf("<DFTD3vdW_init> a1=%f, a2=%f\n", a1_dftD, a2_dftD)
        end
        
        for spe = 1:Nspecies
            @printf("<DFTD3vdW_init> %s species rcov=%f\n", "Ga", rcov_dftD[spe])
        end

        for spe = 1:Nspecies
            @printf("<DFTD3vdW_init> %s species r2r4=%f\n", "Ga", r2r4_dftD[spe])
        end

        for spe = 1:Nspecies, jspe = 1:spe
            @printf("<DFTD3vdW_init> %s %s species r0ab=%f\n", "Ga", "Ga", r0ab_dftD[spe, jspe])
        end

        for spe = 1:Nspecies, jspe = 1:spe, i = 1:maxcn_dftD[spe], j = 1:maxcn_dftD[jspe]
            if C6ab_dftD[spe][jspe][i][j][1] >= 0.0
                @printf("<DFTD3vdW_init> CN_%s[%d]=%f, CN_%s[%d]=%f, C6ab=%f\n", 
                            "Ga", i, C6ab_dftD[spe][jspe][i][j][2], 
                            "Ga", j, C6ab_dftD[spe][jspe][i][j][3], C6ab_dftD[spe][jspe][i][j][1])
            end
        end
    end

  
    n = zeros(Float64, 3)
    for xyz = 1:3
        if xyz == 1
            i = 2
            j = 3
        elseif xyz == 2
            i = 3
            j = 1
        elseif xyz == 3
            i = 1
            j = 2
        end

        a = Latvecs[xyz,:]
        b = Latvecs[i,:]
        c = Latvecs[j,:]

        v = cross(b, c)
        coef = 1/norm(v)
        v = coef*v

        n[xyz] = abs(cncut_dftD/dot(v,a))
    end


    n_CN_DFT_D = zeros(Int64, 3)
    for i = 1:3
        n_CN_DFT_D[i] = ceil(Int64, n[i])
        if DFTD_IntDir[i] == 0
            n_CN_DFT_D[i] = 0
        end
    end


    for xyz = 1:3
        if xyz == 1
            i = 2
            j = 3
        elseif xyz == 2
            i = 3
            j = 1
        elseif xyz == 3
            i = 1
            j = 2
        end

        a = Latvecs[xyz,:]
        b = Latvecs[i,:]
        c = Latvecs[j,:]

        v = cross(b, c)
        coef = 1/norm(v)
        v = coef*v

        n[xyz] = abs(rcut_dftD/dot(v,a))
    end

    n_DFT_D = zeros(Int64, 3)
    for i = 1:3
        n_DFT_D[i] = ceil(Int64, n[i])
        if DFTD_IntDir[i] == 0
            n_DFT_D[i] = 0
        end
    end

    if verbosity >= 1
        @printf("<DFTD3vdW_init> n1=%2d n2=%2d n3=%2d\n", n_DFT_D[1], n_DFT_D[2], n_DFT_D[3])
        @printf("<DFTD3vdW_init> n1_CN=%2d n2_CN=%2d n3_CN=%2d\n", n_CN_DFT_D[1], n_CN_DFT_D[2], n_CN_DFT_D[3])
    end


    return 
end



# n2_CN_DFT_D, rcut_dftD, cncut_dftD, rcovab_dftD, C6ab_dftD, maxcn_dftD, r0ab_dftD, r2r4ab_dftD
function Calc_Energy_DFTD3(Natom, Latvecs, Gxyz, atom2spe, Atoms_vdW_dir, dftd3vdw_params::DFTD3vdW_params)
    
    DFTD3_damp_dftD = dftd3vdw_params.DFTD3_damp_dftD
    sr8_dftD = dftd3vdw_params.sr8_dftD
    alp6_dftD = dftd3vdw_params.alp6_dftD
    alp8_dftD = dftd3vdw_params.alp8_dftD
    cncut_dftD = dftd3vdw_params.cncut_dftD
    rcut_dftD = dftd3vdw_params.rcut_dftD
    s6_dftD = dftd3vdw_params.s6_dftD
    s8_dftD = dftd3vdw_params.s8_dftD
    sr6_dftD = dftd3vdw_params.sr6_dftD
    a1_dftD = dftd3vdw_params.a1_dftD
    a2_dftD = dftd3vdw_params.a2_dftD
    rcovab_dftD = dftd3vdw_params.rcovab_dftD
    C6ab_dftD = dftd3vdw_params.C6ab_dftD
    maxcn_dftD = dftd3vdw_params.maxcn_dftD
    r0ab_dftD = dftd3vdw_params.r0ab_dftD
    r2r4ab_dftD = dftd3vdw_params.r2r4ab_dftD
    n_CN_DFT_D = dftd3vdw_params.n_CN_DFT_D
    n1_CN_DFT_D, n2_CN_DFT_D, n3_CN_DFT_D = n_CN_DFT_D
    rcut2  = rcut_dftD^2
    cncut2 = cncut_dftD^2


    
    
    EdftD    = 0.0
    CN    = zeros(Float64, Natom+1)
    dC6ij = zeros(Float64, Natom+1, Natom+1)
    dEC0  = zeros(Float64, Natom+1, Natom+1)
    dCN = zeros(Float64, Natom+1, Natom+1, 2*n1_CN_DFT_D+1, 2*n2_CN_DFT_D+1, 2*n3_CN_DFT_D+1)


    DFTD3_Force = zeros(Float64, Natom, 3)


    # Compute coordination numbers CN_A and derivative dCN_AB/dr_AB
    for atom = 1:Natom
        spe  = atom2spe[atom]
        xn = 0.0
        for jatom = 1:Natom
            jspe = atom2spe[jatom]
            per_flagB = Atoms_vdW_dir[jatom]

            rij0_1 = Gxyz[atom][1] - Gxyz[jatom][1]
            rij0_2 = Gxyz[atom][2] - Gxyz[jatom][2]
            rij0_3 = Gxyz[atom][3] - Gxyz[jatom][3]

            if per_flagB == 0
                n1_max = 0
                n2_max = 0
                n3_max = 0
            else
                n1_max = n1_CN_DFT_D
                n2_max = n2_CN_DFT_D
                n3_max = n3_CN_DFT_D
            end

            for n1 = -n1_max:n1_max, n2 = -n2_max:n2_max, n3 = -n3_max:n3_max

                rij1 = rij0_1 - (n1*Latvecs[1,1] + n2*Latvecs[2,1] + n3*Latvecs[3,1])
                rij2 = rij0_2 - (n1*Latvecs[1,2] + n2*Latvecs[2,2] + n3*Latvecs[3,2])
                rij3 = rij0_3 - (n1*Latvecs[1,3] + n2*Latvecs[2,3] + n3*Latvecs[3,3])
                dist2 = rij1^2 + rij2^2 + rij3^2

                i1 = n1 + n1_CN_DFT_D + 1
                i2 = n2 + n2_CN_DFT_D + 1
                i3 = n3 + n3_CN_DFT_D + 1

                if dist2 < cncut2 && dist2 > 0.1

                    dist = sqrt(dist2)

                    exparg = -k1_dftD*(rcovab_dftD[spe,jspe]/dist - 1.0)
                    expval = exp(exparg)
                    fdamp  = 1/(1 + expval)
                    xn += fdamp

                    dCN[atom][jatom][i1][i2][i3] = -fdamp^2*expval*k1_dftD*rcovab_dftD[spe,jspe]/dist2
                else
                    dCN[atom][jatom][i1][i2][i3] = 0.0
                end
            end
        end

        CN[atom] = xn
    end



    # ------------------------------------------------------------
    # Calculate energy and collect gradients of two-body terms
    # ------------------------------------------------------------

    for atom = 1:Natom
        spe  = atom2spe[atom]
        dblcnt_factor = 0.5

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0
        Ene = 0.0

        per_flagA = Atoms_vdW_dir[atom]

        for jatom = 1:Natom
            jspe = atom2spe[jatom]

            dEC0[atom,jatom] = 0.0

            per_flagB = Atoms_vdW_dir[jatom]

            rij0_1 = Gxyz[atom][1] - Gxyz[jatom][1]
            rij0_2 = Gxyz[atom][2] - Gxyz[jatom][2]
            rij0_3 = Gxyz[atom][3] - Gxyz[jatom][3]

            # Calculate C6, C8 coefficient and derivatives

            Z   = 0.0
            W   = 0.0
            dZi = 0.0
            dZj = 0.0
            dWi = 0.0
            dWj = 0.0

            for i = 1:maxcn_dftD[spe],  j = 1:maxcn_dftD[jspe]

                C6_ref = C6ab_dftD[spe][jspe][i][j][1]

                if C6_ref > 1.0e-12

                    dAi = CN[atom] - C6ab_dftD[spe][jspe][i][j][2]
                    dBj = CN[jatom] - C6ab_dftD[spe][jspe][i][j][3]

                    exparg = -k3_dftD*(dAi^2 + dBj^2)
                    Lij = exp(exparg)

                    Z += C6_ref*Lij
                    W += Lij

                    dZi += C6_ref*Lij*2.0*k3_dftD*dAi
                    dZj += C6_ref*Lij*2.0*k3_dftD*dBj

                    dWi += Lij*2.0*k3_dftD*dAi
                    dWj += Lij*2.0*k3_dftD*dBj
                end

                if W > 1.0e-12
                    C6 = Z/W
                    C8 = 3.0*C6*r2r4ab_dftD[spe,jspe]
                    C8C6 = 3.0*r2r4ab_dftD[spe,jspe]
                    dC6ij[atom, Gc_BN] = (dZi*W - dWi*Z)/W^2
                else
                    C6 = 0.0
                    C8 = 0.0
                    C8C6 = 3.0*r2r4ab_dftD[spe,jspe]
                    dC6ij[atom, Gc_BN] = 0.0
                end


                # Calculate energy and first part of gradients
                if per_flagB == 0
                    n1_max = 0
                    n2_max = 0
                    n3_max = 0
                else
                    n1_max = n1_DFT_D
                    n2_max = n2_DFT_D
                    n3_max = n3_DFT_D
                end

                for n1 = -n1_max:n1_max, n2 = -n2_max:n2_max, n3 = -n3_max:n3_max

                    # double counting factor
                    if (abs(n1) + abs(n2) + abs(n3)) ≠ 0 && per_flagA == 0 && per_flagB == 1
                        dblcnt_factor = 1.0
                    else
                        dblcnt_factor = 0.5
                    end

                    rij1 = rij0_1 - (n1*Latvecs[1,1] + n2*Latvecs[2,1] + n3*Latvecs[3,1])
                    rij2 = rij0_2 - (n1*Latvecs[1,2] + n2*Latvecs[2,2] + n3*Latvecs[3,2])
                    rij3 = rij0_3 - (n1*Latvecs[1,3] + n2*Latvecs[2,3] + n3*Latvecs[3,3])
                    dist2 = rij1^2 + rij2^2 + rij3^2

                    if 0.1 < dist2 && dist2 < rcut2

                        dist  = sqrt(dist2)
                        dist5 = dist2*dist2*dist
                        dist6 = dist2*dist2*dist2
                        dist7 = dist6*dist
                        dist8 = dist6*dist2

                        dE6 = 0.0
                        dE8 = 0.0

                        # DFT-D3 zero damping
                        if DFTD3_damp_dftD == 1
                            # E6
                            powarg = dist/(sr6_dftD*r0ab_dftD[jspe,spe])
                            powval = powarg^(-alp6_dftD)

                            fdamp6 = 1.0/(1.0 + 6.0*powval)

                            E -= dblcnt_factor*s6_dftD*C6*fdamp6/dist6

                            dE6 = (s6_dftD*C6*fdamp6/dist6)*(6.0/dist)*(-1.0 + alp6_dftD * powval * fdamp6)

                            # E8
                            powarg = dist/(sr8_dftD*r0ab_dftD[jspe,spe])
                            powval = powarg^(-alp8_dftD)

                            fdamp8 = 1.0/(1.0 + 6.0*powval)

                            E -= dblcnt_factor*s8_dftD*C8*fdamp8/dist8
                            dE8 = (s8_dftD*C8*fdamp8/dist8)*(2.0/dist)*(-4.0 + 3.0*alp8_dftD*powval*fdamp8)

                            dEC0[atom, Gc_BN] += s6_dftD*fdamp6/dist6 + s8_dftD*3.0*r2r4ab_dftD[spe,jspe]*fdamp8/dist8

                        elseif DFTD3_damp_dftD == 2
                        # DFT-D3 BJ damping
                                        
                            fdamp = a1_dftD*sqrt(C8C6) + a2_dftD
                            fdamp6 = fdamp^6
                            fdamp8 = fdamp^8

                            t6  = dist6 + fdamp6
                            t62 = t6^2

                            t8  = dist8 + fdamp8
                            t82 = t8^2

                            E -= dblcnt_factor*s6_dftD*C6/t6
                            dE6 = -s6_dftD*C6*6.0*dist5/t62

                            E -= dblcnt_factor*s8_dftD*C8/t8
                            dE8 = -s8_dftD*C8*8.0*dist7/t82

                            dEC0[atom,jatom] += s6_dftD/t6 + s8_dftD*3.0*r2r4ab_dftD[spe,jspe]/t8
                        end

                        dEtot = dE6 + dE8

                        dEx -= dEtot*rij1/dist
                        dEy -= dEtot*rij2/dist
                        dEz -= dEtot*rij3/dist
                    end
                end
            end

            EdftD += E


            # Gradients from two-body terms
            DFTD3_Force[atom,1] += dEx
            DFTD3_Force[atom,2] += dEy
            DFTD3_Force[atom,3] += dEz
        end
    end



    for atom = 1:Natom
        spe  = atom2spe[atom]

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for jatom = 1:Natom

            jspe = atom2spe[jatom]
            per_flagB = Atoms_vdW_dir[jatom]

            rij0_1 = Gxyz[atom][1] - Gxyz[jatom][1]
            rij0_2 = Gxyz[atom][2] - Gxyz[jatom][2]
            rij0_3 = Gxyz[atom][3] - Gxyz[jatom][3]

            if per_flagB == 0
                n1_max = 0
                n2_max = 0
                n3_max = 0
            else
                n1_max = n1_CN_DFT_D
                n2_max = n2_CN_DFT_D
                n3_max = n3_CN_DFT_D
            end

            for n1 = -n1_max:n1_max, n2 = -n2_max:n2_max, n3 = -n3_max:n3_max

                dEC = 0.0
                rij1 = rij0_1 - (n1*Latvecs[1,1] + n2*Latvecs[2,1] + n3*Latvecs[3,1])
                rij2 = rij0_2 - (n1*Latvecs[1,2] + n2*Latvecs[2,2] + n3*Latvecs[3,2])
                rij3 = rij0_3 - (n1*Latvecs[1,3] + n2*Latvecs[2,3] + n3*Latvecs[3,3])
                dist2 = rij1^2 + rij2^2 + rij3^2

                if 0.1 < dist2 && dist2 < cncut2

                    dist = sqrt(dist2)

                    i1 = n1 + n1_CN_DFT_D + 1
                    i2 = n2 + n2_CN_DFT_D + 1
                    i3 = n3 + n3_CN_DFT_D + 1

                    dCN_AB = dCN[atom,jatom,i1,i2,i3]

                    for katom = 1:Natom
                        dEC += dEC0[atom,Gc_CN]*dC6ij[atom,Gc_CN]*dCN_AB
                        dEC += dEC0[Gc_BN,Gc_CN]*dC6ij[Gc_BN,Gc_CN]*dCN_AB
                    end

                    dEx += dEC*rij1/dist
                    dEy += dEC*rij2/dist
                    dEz += dEC*rij3/dist
                end
            end

            DFTD3_Force[atom,1] += dEx
            DFTD3_Force[atom,2] += dEy
            DFTD3_Force[atom,3] += dEz
        end
    end


    return EdftD, DFTD3_Force
end


# Graphite2
Natom = 4
Nspecies = 1
atom2spe = [1,1,1,1]
Latvecs = [0.0 2.48901587 0.0; 2.15555097384263 1.244507935 0.0; 0.0 0.0 20.0]*1.8897259886
Spe_Znumber = [4]
DFTD_IntDir = [0,0,1]
Gxyz = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
Atoms_vdW_dir = ones(Int64, Natom)
dftd3vdw_params = DFTD3vdW_params()
n2_CN_DFT_D, rcut_dftD, cncut_dftD, rcovab_dftD, C6ab_dftD, maxcn_dftD, r0ab_dftD, r2r4ab_dftD = DFTD3vdW_init(Nspecies, Latvecs, Spe_Znumber, DFTD_IntDir)
Calc_Energy_DFTD3(Natom, Latvecs, Gxyz, atom2spe, Atoms_vdW_dir, dftd3vdw_params)