function Calc_Diff_Coef(gLatvecs::Matrix{Float64})

    N = 27
    Mat = zeros(Float64, N, N)
    p = 1
    for i1 = -1:1, j1 = -1:1, k1 = -1:1
        x = i1*gLatvecs[1,1] + j1*gLatvecs[2,1] + k1*gLatvecs[3,1]
        y = i1*gLatvecs[1,2] + j1*gLatvecs[2,2] + k1*gLatvecs[3,2]
        z = i1*gLatvecs[1,3] + j1*gLatvecs[2,3] + k1*gLatvecs[3,3]

        q = 1
        for i2 = 0:2, j2 = 0:2, k2 = 0:2
            x1 = x^i2
            y1 = y^j2
            z1 = z^k2

            Mat[p,q] = x1*y1*z1
            q += 1
        end
        p += 1
    end


    invMat = inv(Mat)

    Diff_Coef = Vector{Vector{Float64}}(undef, 3)
    for i = 1:3
        Diff_Coef[i] = zeros(Float64, N)        
    end

    for j = 1:N
        Diff_Coef[1][j] = invMat[10,j]
        Diff_Coef[2][j] = invMat[4,j]
        Diff_Coef[3][j] = invMat[2,j]
    end


    return Diff_Coef
end


function transform2primcell(Ngrid, i)
    if 0 <= i <= Ngrid-1
        return i
    elseif i <= -1
        return i + Ngrid
    elseif i >= Ngrid
        return i - Ngrid
    end
end


function Set_dDensity_Grid!(Ngrid, SpinPol, Den, pccDen, Diff_Coef, dDensity_Grid)

    if SpinPol == "off"
        _Set_dDensity_Grid_woSpin!(Ngrid, Den, pccDen, Diff_Coef, dDensity_Grid)
    elseif SpinPol ∈ ("on", "nc")
        _Set_dDensity_Grid_wSpin!(Ngrid, Den, pccDen, Diff_Coef, dDensity_Grid)
    else
        println("SpinPol = $(SpinPol) in XC/GGA_PBE.jl")
        error("please check SpinPol")
    end
end


function _Set_dDensity_Grid_woSpin!(Ngrid, Den, pccDen, Diff_Coef, dDensity_Grid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    fden = zeros(Float64, 27)

    xyz = 0
    for i = 0:Ngrid1-1, j = 0:Ngrid2-1, k = 0:Ngrid3-1
        xyz += 1
        if den_min<(2*Den[1][xyz])

            p = 0
            for i1 = -1:1
                i2 = transform2primcell(Ngrid1, i+i1)
                for j1 = -1:1
                    j2 = transform2primcell(Ngrid2, j+j1)
                    for k1 = -1:1
                        k2 = transform2primcell(Ngrid3, k+k1)
                        
                        p += 1
                        xyz2 = i2*Ngrid2*Ngrid3 + j2*Ngrid3 + k2 + 1
                        fden[p] = Den[1][xyz2] + pccDen[xyz2]
                    end
                end
            end

            up_x = dot(fden, Diff_Coef[1])
            up_y = dot(fden, Diff_Coef[2])
            up_z = dot(fden, Diff_Coef[3])

            dDensity_Grid[1][1][xyz] = up_x
            dDensity_Grid[1][2][xyz] = up_y
            dDensity_Grid[1][3][xyz] = up_z
        else
            dDensity_Grid[1][1][xyz] = 0.0
            dDensity_Grid[1][2][xyz] = 0.0
            dDensity_Grid[1][3][xyz] = 0.0
        end
    end
end


function _Set_dDensity_Grid_wSpin!(Ngrid, Den, pccDen, Diff_Coef, dDensity_Grid)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    fden_up = zeros(Float64, 27)
    fden_dn = zeros(Float64, 27)

    xyz = 0
    for i = 0:Ngrid1-1, j = 0:Ngrid2-1, k = 0:Ngrid3-1
        xyz += 1
        if den_min<(Den[1][xyz]+Den[2][xyz])

            p = 0
            for i1 = -1:1
                i2 = transform2primcell(Ngrid1, i+i1)
                for j1 = -1:1
                    j2 = transform2primcell(Ngrid2, j+j1)
                    for k1 = -1:1
                        k2 = transform2primcell(Ngrid3, k+k1)
                        
                        p += 1
                        xyz2 = i2*Ngrid2*Ngrid3 + j2*Ngrid3 + k2 + 1
                        fden_up[p] = Den[1][xyz2] + pccDen[xyz2]
                        fden_dn[p] = Den[2][xyz2] + pccDen[xyz2]
                    end
                end
            end

            up_x = dot(fden_up, Diff_Coef[1])
            up_y = dot(fden_up, Diff_Coef[2])
            up_z = dot(fden_up, Diff_Coef[3])

            dn_x = dot(fden_dn, Diff_Coef[1])
            dn_y = dot(fden_dn, Diff_Coef[2])
            dn_z = dot(fden_dn, Diff_Coef[3])

            dDensity_Grid[1][1][xyz] = up_x
            dDensity_Grid[1][2][xyz] = up_y
            dDensity_Grid[1][3][xyz] = up_z

            dDensity_Grid[2][1][xyz] = dn_x
            dDensity_Grid[2][2][xyz] = dn_y
            dDensity_Grid[2][3][xyz] = dn_z
        else
            dDensity_Grid[1][1][xyz] = 0.0
            dDensity_Grid[1][2][xyz] = 0.0
            dDensity_Grid[1][3][xyz] = 0.0
            dDensity_Grid[2][1][xyz] = 0.0
            dDensity_Grid[2][2][xyz] = 0.0
            dDensity_Grid[2][3][xyz] = 0.0
        end
    end
end


function Set_Vxc_GGA_PBE!(Ngrid, SpinPol, Den, pccDen, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc)
    
    if SpinPol == "off"
        _Set_Vxc_GGA_PBE_woSpin!(Ngrid, Den, pccDen, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc)
    elseif SpinPol ∈ ("on", "nc")
        _Set_Vxc_GGA_PBE_wSpin!(Ngrid, Den, pccDen, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc)
    else
        println("SpinPol = $(SpinPol) in XC/GGA_PBE.jl")
        error("please check SpinPol")
    end


    if SpinPol == "nc"
        NN = prod(Ngrid)
        for xyz = 1:NN
            tmp0 = 0.5*(Vxc[1][xyz] + Vxc[2][xyz])
            tmp1 = 0.5*(Vxc[1][xyz] - Vxc[2][xyz])
            theta = Den[3][xyz]
            phi = Den[4][xyz]

            sit = sin(theta)
            cot = cos(theta)
            sip = sin(phi)
            cop = cos(phi)

            Vxc[4][xyz] = -tmp1*sit*sip
            Vxc[3][xyz] = -tmp1*sit*cop
            Vxc[2][xyz] = tmp0 - cot*tmp1
            Vxc[1][xyz] = tmp0 + cot*tmp1
        end
    end
end


function _Set_Vxc_GGA_PBE_woSpin!(Ngrid, Den, pccDen, dDensity_Grid, Vxc)

    NN = prod(Ngrid)
    Exc = zeros(Float64, 2)
    ED = zeros(Float64, 2)
    DEXDD = zeros(Float64, 2)
    DECDD = zeros(Float64, 2)
    GDENS = zeros(Float64, 6)
    DEXDGD = zeros(Float64, 6)
    DECDGD = zeros(Float64, 6)

    fill!(Vxc[1], 0.0)

    for xyz = 1:NN

        ED[1] = Den[1][xyz]
        ED[2] = Den[1][xyz]

        if 2*ED[1] > den_min

            GDENS[1] = dDensity_Grid[1][1][xyz]
            GDENS[2] = dDensity_Grid[1][2][xyz]
            GDENS[3] = dDensity_Grid[1][3][xyz]
            GDENS[4] = dDensity_Grid[1][1][xyz]
            GDENS[5] = dDensity_Grid[1][2][xyz]
            GDENS[6] = dDensity_Grid[1][3][xyz]

            ED[1] += pccDen[xyz]
	        ED[2] += pccDen[xyz]

            XC_PBE!(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD)

            Vxc[1][xyz] = Exc[1] + Exc[2]
        end
    end
end


function _Set_Vxc_GGA_PBE_wSpin!(Ngrid, Den, pccDen, dDensity_Grid, Vxc)

    NN = prod(Ngrid)
    Exc = zeros(Float64, 2)
    ED = zeros(Float64, 2)
    DEXDD = zeros(Float64, 2)
    DECDD = zeros(Float64, 2)
    GDENS = zeros(Float64, 6)
    DEXDGD = zeros(Float64, 6)
    DECDGD = zeros(Float64, 6)

    fill!(Vxc[1], 0.0)
    fill!(Vxc[2], 0.0)

    for xyz = 1:NN

        ED[1] = Den[1][xyz]
        ED[2] = Den[2][xyz]

        if (ED[1] + ED[2]) > den_min

            GDENS[1] = dDensity_Grid[1][1][xyz]
            GDENS[2] = dDensity_Grid[1][2][xyz]
            GDENS[3] = dDensity_Grid[1][3][xyz]
            GDENS[4] = dDensity_Grid[2][1][xyz]
            GDENS[5] = dDensity_Grid[2][2][xyz]
            GDENS[6] = dDensity_Grid[2][3][xyz]

            ED[1] += pccDen[xyz]
	        ED[2] += pccDen[xyz]

            XC_PBE!(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD)

            Vxc[1][xyz] = Exc[1] + Exc[2]
            Vxc[2][xyz] = Exc[1] + Exc[2]
        end
    end
end


function _Set_Vxc_GGA_PBE_woSpin!(Ngrid, Den, pccDen, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    
    fden_x = zeros(Float64, 27)
    fden_y = zeros(Float64, 27)
    fden_z = zeros(Float64, 27)
    Exc = zeros(Float64, 2)
    ED = zeros(Float64, 2)
    DEXDD = zeros(Float64, 2)
    DECDD = zeros(Float64, 2)
    GDENS = zeros(Float64, 6)
    DEXDGD = zeros(Float64, 6)
    DECDGD = zeros(Float64, 6)

    fill!(Vxc[1], 0.0)
    fill!(dEXC_dGD[1][1], 0.0)
    fill!(dEXC_dGD[1][2], 0.0)
    fill!(dEXC_dGD[1][3], 0.0)


    for xyz = 1:NN

        ED[1] = Den[1][xyz]
        ED[2] = Den[1][xyz]

        if 2*ED[1] >= den_min

            GDENS[1] = dDensity_Grid[1][1][xyz]
            GDENS[2] = dDensity_Grid[1][2][xyz]
            GDENS[3] = dDensity_Grid[1][3][xyz]
            GDENS[4] = dDensity_Grid[1][1][xyz]
            GDENS[5] = dDensity_Grid[1][2][xyz]
            GDENS[6] = dDensity_Grid[1][3][xyz]

            ED[1] += pccDen[xyz]
	        ED[2] += pccDen[xyz]

            XC_PBE!(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD)

            Vxc[1][xyz] = DEXDD[1] + DECDD[1]

            # up spin
            dEXC_dGD[1][1][xyz] = DEXDGD[1] + DECDGD[1]
            dEXC_dGD[1][2][xyz] = DEXDGD[2] + DECDGD[2]
            dEXC_dGD[1][3][xyz] = DEXDGD[3] + DECDGD[3]
        end
    end


    xyz = 0
    for i = 0:Ngrid1-1, j = 0:Ngrid2-1, k = 0:Ngrid3-1
        xyz += 1
        if 2*Den[1][xyz] > den_min

            p = 0
            for i1 = -1:1
                i2 = transform2primcell(Ngrid1, i+i1)
                for j1 = -1:1
                    j2 = transform2primcell(Ngrid2, j+j1)
                    for k1 = -1:1
                        k2 = transform2primcell(Ngrid3, k+k1)
                        
                        p += 1
                        xyz2 = i2*Ngrid2*Ngrid3 + j2*Ngrid3 + k2 + 1
                        fden_x[p] = dEXC_dGD[1][1][xyz2]
                        fden_y[p] = dEXC_dGD[1][2][xyz2]
                        fden_z[p] = dEXC_dGD[1][3][xyz2]
                    end
                end
            end

            tmp  = dot(fden_x, Diff_Coef[1])
            tmp += dot(fden_y, Diff_Coef[2])
            tmp += dot(fden_z, Diff_Coef[3])

            # XC potential
            Vxc[1][xyz] -= tmp
        end
    end
end


function _Set_Vxc_GGA_PBE_wSpin!(Ngrid, Den, pccDen, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc)

    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    
    fden_x_up = zeros(Float64, 27)
    fden_y_up = zeros(Float64, 27)
    fden_z_up = zeros(Float64, 27)
    fden_x_dn = zeros(Float64, 27)
    fden_y_dn = zeros(Float64, 27)
    fden_z_dn = zeros(Float64, 27)

    Exc = zeros(Float64, 2)
    ED = zeros(Float64, 2)
    DEXDD = zeros(Float64, 2)
    DECDD = zeros(Float64, 2)
    GDENS = zeros(Float64, 6)
    DEXDGD = zeros(Float64, 6)
    DECDGD = zeros(Float64, 6)

    fill!(Vxc[1], 0.0)
    fill!(Vxc[2], 0.0)
    fill!(dEXC_dGD[1][1], 0.0)
    fill!(dEXC_dGD[1][2], 0.0)
    fill!(dEXC_dGD[1][3], 0.0)
    fill!(dEXC_dGD[2][1], 0.0)
    fill!(dEXC_dGD[2][2], 0.0)
    fill!(dEXC_dGD[2][3], 0.0)


    for xyz = 1:NN

        ED[1] = Den[1][xyz]
        ED[2] = Den[2][xyz]

        if ED[1] + ED[2] > den_min

            GDENS[1] = dDensity_Grid[1][1][xyz]
            GDENS[2] = dDensity_Grid[1][2][xyz]
            GDENS[3] = dDensity_Grid[1][3][xyz]
            GDENS[4] = dDensity_Grid[2][1][xyz]
            GDENS[5] = dDensity_Grid[2][2][xyz]
            GDENS[6] = dDensity_Grid[2][3][xyz]

            ED[1] += pccDen[xyz]
	        ED[2] += pccDen[xyz]

            XC_PBE!(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD)

            Vxc[1][xyz] = DEXDD[1] + DECDD[1]
            Vxc[2][xyz] = DEXDD[2] + DECDD[2]

            # up spin
            dEXC_dGD[1][1][xyz] = DEXDGD[1] + DECDGD[1]
            dEXC_dGD[1][2][xyz] = DEXDGD[2] + DECDGD[2]
            dEXC_dGD[1][3][xyz] = DEXDGD[3] + DECDGD[3]

            # down spin
            dEXC_dGD[2][1][xyz] = DEXDGD[4] + DECDGD[4]
            dEXC_dGD[2][2][xyz] = DEXDGD[5] + DECDGD[5]
            dEXC_dGD[2][3][xyz] = DEXDGD[6] + DECDGD[6]
        end
    end


    xyz = 0
    for i = 0:Ngrid1-1, j = 0:Ngrid2-1, k = 0:Ngrid3-1
        xyz += 1
        if (Den[1][xyz] + Den[2][xyz]) > den_min

            p = 0
            for i1 = -1:1
                i2 = transform2primcell(Ngrid1, i+i1)
                for j1 = -1:1
                    j2 = transform2primcell(Ngrid2, j+j1)
                    for k1 = -1:1
                        k2 = transform2primcell(Ngrid3, k+k1)
                        
                        p += 1
                        xyz2 = i2*Ngrid2*Ngrid3 + j2*Ngrid3 + k2 + 1
                        fden_x_up[p] = dEXC_dGD[1][1][xyz2]
                        fden_y_up[p] = dEXC_dGD[1][2][xyz2]
                        fden_z_up[p] = dEXC_dGD[1][3][xyz2]

                        fden_x_dn[p] = dEXC_dGD[2][1][xyz2]
                        fden_y_dn[p] = dEXC_dGD[2][2][xyz2]
                        fden_z_dn[p] = dEXC_dGD[2][3][xyz2]
                    end
                end
            end

            tmp_up  = dot(fden_x_up, Diff_Coef[1])
            tmp_up += dot(fden_y_up, Diff_Coef[2])
            tmp_up += dot(fden_z_up, Diff_Coef[3])

            tmp_dn  = dot(fden_x_dn, Diff_Coef[1])
            tmp_dn += dot(fden_y_dn, Diff_Coef[2])
            tmp_dn += dot(fden_z_dn, Diff_Coef[3])

            # XC potential
            Vxc[1][xyz] -= tmp_up
            Vxc[2][xyz] -= tmp_dn
        end
    end
end
