function Generate_HermOrbs(material, ucell, Orbs_Grid, kmesh_c::Integer, kpts_ab::Vector{Float64})

    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    SpinPol = material.SpinPol
    Latvecs = material.Latvecs
    Natom = material.Natom
    Gxyz = material.Gxyz
    atv_ijk = material.atv_ijk
    ChemP = material.ChemP
    Total_SpinS = material.Total_SpinS
    Valence_Electrons = material.Valence_Electrons
    CellVolume = material.cellVolume
    E_Temp = material.E_Temp
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks

    
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    Grid_Origin = ucell.system_grid.Grid_Origin
    Ngrid = ucell.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    Nkpt = 1

    if SpinPol == "off"
        spinsize = 1
        Spindeg = 1
        Nfsize = fsize
        fsize3 = fsize + 2
        fsize4 = Int64(div(Valence_Electrons,2))
        fsize5 = fsize4*kmesh_c
    elseif SpinPol == "on"
        spinsize = 2
        Spindeg = 2
        Nfsize = fsize
        fsize3 = fsize + 2
        fsize4 = Int64(div(Valence_Electrons,2)) + Int64(fabs(floor(Total_SpinS)))*2 + 1
        fsize5 = fsize4*kmesh_c
    elseif SpinPol == "nc"
        spinsize = 1
        Spindeg = 2
        Nfsize = 2*fsize
        fsize3 = 2*fsize + 2
        fsize4 = Int64(Valence_Electrons)
        fsize5 = fsize4*kmesh_c
    end

    NHOMO = fsize4


    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    i = 0
    list_ik_mu = zeros(Int32, fsize5, 2)
    for mu = 1:NHOMO, ikc = 1:kmesh_c
        i += 1
        list_ik_mu[i,1] = mu
        list_ik_mu[i,2] = ikc
    end

    Smatrix = zeros(ComplexF64, fsize5, fsize5)



    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, kmesh_c)
        for ik = 1:kmesh_c
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end


    if SpinPol ∈ ("off", "on")
        S = zeros(ComplexF64, fsize, fsize)
        H = zeros(ComplexF64, fsize, fsize)
        for spin = 1:spinsize, ik = 1:kmesh_c
            HS_matrix_scfout!(S, OLP, material, [0.0,0.0,(ik-1)/kmesh_c])
            HS_matrix_scfout!(H, Hks[spin], material, [0.0,0.0,(ik-1)/kmesh_c])
            Cnk[spin][ik] = eigvecs(Hermitian(H), Hermitian(S))
        end
    elseif SpinPol == "nc"
        tmpH = zeros(ComplexF64, fsize, fsize)
        S = zeros(ComplexF64, 2*fsize, 2*fsize)
        H = zeros(ComplexF64, 2*fsize, 2*fsize)
        @inbounds for ik = 1:kmesh_c
            HS_matrix_NC_scfout!(tmpH, H, Hks, iHks, material, [kpts_ab[1],kpts_ab[2],(ik-1)/kmesh_c])
            HS_matrix_scfout!(tmpH, OLP, material, [kpts_ab[1],kpts_ab[2],(ik-1)/kmesh_c])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Cnk[1][ik] = eigvecs(Hermitian(H), Hermitian(S))
        end
    end

    


    
    HOMOs_Coef = Vector{Vector{Vector{Vector{Vector{ComplexF64}}}}}(undef, kmesh_c)
    for ik = 1:kmesh_c
        HOMOs_Coef[ik] = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, Spindeg)
        for spin = 1:Spindeg
            HOMOs_Coef[ik][spin] = Vector{Vector{Vector{ComplexF64}}}(undef, NHOMO)
            for μ = 1:NHOMO
                HOMOs_Coef[ik][spin][μ] = Vector{Vector{ComplexF64}}(undef, Natom)
                for atom = 1:Natom
                    HOMOs_Coef[ik][spin][μ][atom] = zeros(ComplexF64, Total_NumOrbs[atom])
                end
            end
        end
    end


    if SpinPol ∈ ("off", "on")
        for ik = 1:kmesh_c, spin = 1:Spindeg, μ = 1:NHOMO
            for atom = 1:Natom
                Anum = MP[atom]
                for ist = 1:Total_NumOrbs[atom]
                    HOMOs_Coef[ik][spin][μ][atom][ist] = Cnk[spin][ik][Anum+ist,μ]
                end
            end
        end
    else
        for ik = 1:kmesh_c, μ = 1:NHOMO
            for atom = 1:Natom
                Anum = MP[atom]
                for ist = 1:Total_NumOrbs[atom]
                    HOMOs_Coef[ik][1][μ][atom][ist] = Cnk[1][ik][Anum+ist,μ]
                    HOMOs_Coef[ik][2][μ][atom][ist] = Cnk[1][ik][Anum+ist+fsize,μ]
                end
            end
        end
    end



    unk_r1 = zeros(ComplexF64, NN)
    unk_r2 = zeros(ComplexF64, NN)


    kpts1 = zeros(Float64, 3)
    kpts2 = zeros(Float64, 3)

    kpts1[2] = kpts_ab[1]
    kpts1[3] = kpts_ab[2]

    kpts2[2] = kpts_ab[1]
    kpts2[3] = kpts_ab[2]

    for ist = 1:fsize5, jst = 1:fsize5

        mu1 = list_ik_mu[ist,1]
        ik1 = list_ik_mu[ist,2]
        mu2 = list_ik_mu[jst,1]
        ik2 = list_ik_mu[jst,2]

        kpts1[1] = (ik1-1)/kmesh_c
        kpts2[1] = (ik2-1)/kmesh_c

        if ik2-ik1 == 1

            k1_x = kpts1[1]*Recvecs[1,1] + kpts1[2]*Recvecs[2,1] + kpts1[3]*Recvecs[3,1]
            k1_y = kpts1[1]*Recvecs[1,2] + kpts1[2]*Recvecs[2,2] + kpts1[3]*Recvecs[3,2]
            k1_z = kpts1[1]*Recvecs[1,3] + kpts1[2]*Recvecs[2,3] + kpts1[3]*Recvecs[3,3]

            k2_x = kpts2[1]*Recvecs[1,1] + kpts2[2]*Recvecs[2,1] + kpts2[3]*Recvecs[3,1]
            k2_y = kpts2[1]*Recvecs[1,2] + kpts2[2]*Recvecs[2,2] + kpts2[3]*Recvecs[3,2]
            k2_z = kpts2[1]*Recvecs[1,3] + kpts2[2]*Recvecs[2,3] + kpts2[3]*Recvecs[3,3]
            
            # Overlap = Calc_Overlap_unk!(HOMOs_Coef[ik1][1][mu1], kpts1, Natom, GridN_Atom, GridListAtom, CellListAtom, atv_ijk, Total_NumOrbs, Orbs_Grid)
            Calc_unk_r!(unk_r1, HOMOs_Coef[ik1][1][mu1], kpts1, k1_x, k1_y, k1_z, Natom, Grid_Origin, GridN_Atom, GridListAtom, CellListAtom, atv_ijk, Total_NumOrbs, gLatvecs, Ngrid, Orbs_Grid)
            Calc_unk_r!(unk_r2, HOMOs_Coef[ik2][1][mu2], kpts2, k2_x, k2_y, k2_z, Natom, Grid_Origin, GridN_Atom, GridListAtom, CellListAtom, atv_ijk, Total_NumOrbs, gLatvecs, Ngrid, Orbs_Grid)

            # Smatrix[ist,jst] = Overlap/CellVolume
            Smatrix[ist,jst] = dot(unk_r1,unk_r2)/CellVolume
        end
    end


    zjk, wjk = eigen(Hermitian(Smatrix))


    return zjk, wjk
end


function Calc_unk_r!(unk_r, HOMOs_Coef, kpts, kx, ky, kz, Natom, Grid_Origin, GridN_Atom, GridListAtom, CellListAtom, atv_ijk, Total_NumOrbs, gLatvecs, Ngrid, Orbs_Grid)
    fill!(unk_r, 0.0)
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    for atom = 1:Natom, Nc = 1:GridN_Atom[atom]
        GN = GridListAtom[atom][Nc]
        Rn = CellListAtom[atom][Nc]+1

        n1 = div(GN, Ngrid2*Ngrid3)
        n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
        n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

        Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
        Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
        Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

        l1 = -atv_ijk[Rn][1]
        l2 = -atv_ijk[Rn][2]
        l3 = -atv_ijk[Rn][3]

        kRn = kpts[1]*l1 + kpts[2]*l2 + kpts[3]*l3
        ex = exp(im*2*pi*kRn-im*(kx*Cx + ky*Cy + kz*Cz))

        for ist = 1:Total_NumOrbs[atom]
            unk_r[GN+1] += ex*HOMOs_Coef[atom][ist]*Orbs_Grid[atom][ist][Nc]
        end
    end
end


function Calc_Overlap_unk(HOMOs_Coef1, HOMOs_Coef2, kpts, Natom, GridN_Atom, GridListAtom, CellListAtom, atv_ijk, Total_NumOrbs, Orbs_Grid)
    for atom = 1:Natom, Nc = 1:GridN_Atom[atom]
        GN = GridListAtom[atom][Nc]
        Rn = CellListAtom[atom][Nc]+1

        n1 = div(GN, Ngrid2*Ngrid3)
        n2 = div(GN - n1*Ngrid2*Ngrid3, Ngrid3)
        n3 = GN - n1*Ngrid2*Ngrid3 - n2*Ngrid3

        Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
        Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
        Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

        x = Cx + atv[cell][1] - Gxyz[atom][1]
        y = Cy + atv[cell][2] - Gxyz[atom][2]
        z = Cz + atv[cell][3] - Gxyz[atom][3]

        l1 = -atv_ijk[Rn][1]
        l2 = -atv_ijk[Rn][2]
        l3 = -atv_ijk[Rn][3]

        kRn = kpts[1]*l1 + kpts[2]*l2 + kpts[3]*l3
        ex = exp(im*2*pi*kRn)

        ex = exp(-im*(bvector[ib][1]*x + bvector[ib][2]*y + bvector[ib][3]*z))

        for ist = 1:Total_NumOrbs[atom]
            unk_r[GN] += ex*HOMOs_Coef[atom][ist]*Orbs_Grid[atom][ist][Nc]
        end
    end
end