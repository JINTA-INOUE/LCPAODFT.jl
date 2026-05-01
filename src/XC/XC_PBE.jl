function XC_PW92C(dens)

    den_min = 1e-14
    den_min_half = 0.5e-14
    COE1 = 0.6203504908994

    A1,A2,A3 = 0.0310910,0.0155450,0.0168870
    alpha1_1,alpha1_2,alpha1_3 = 0.2137000,0.2054800,0.1112500

    beta1_1,beta1_2,beta1_3 = 7.5957000,14.1189000,10.3570000
    beta2_1,beta2_2,beta2_3 = 3.5876000,6.1977000,3.6231000
    beta3_1,beta3_2,beta3_3 = 1.6382000,3.3662000,0.8802600
    beta4_1,beta4_2,beta4_3 = 0.4929400,0.6251700,0.4967100

    dtot = dens[1] + dens[2]

    if dtot <= den_min
        rs = 6203.504908994
        dens[1] = den_min_half
        dens[2] = den_min_half
        dtot = den_min
    else
        rs = COE1 * dtot^(-1/3)
    end

    tmp0 = 1.0/dtot
    zeta = tmp0*(dens[1]-dens[2])
    if zeta > 0.99
        zeta = 0.99
    end
    if zeta < -1.0
        zeta = -0.99
    end

    drsdd = -(1/3)*rs*tmp0
    dzdd1 = tmp0*(1-zeta)
    dzdd2 = tmp0*(-1-zeta)

    srs = sqrt(rs)

    ####################################
    # G(0)
    ####################################

    b = beta1_1*srs + rs*(beta2_1 + beta3_1*srs + beta4_1*rs)

    dbdrs = beta1_1*0.5/srs + beta2_1 + beta3_1*1.5*srs + beta4_1*2*rs

    c = 1 + 1/(2*A1*b)
    dcdrs = -(c-1)*dbdrs/b

    dum = log(c)
    dum1 = 1 + alpha1_1*rs

    G0 = -2*A1*dum1*dum
    dG0 = -2*A1*(alpha1_1*dum + dum1*dcdrs/c)

    ####################################
    # G(1)
    ####################################

    b = beta1_2*srs + rs*(beta2_2 + beta3_2*srs + beta4_2*rs)

    dbdrs = beta1_2*0.5/srs + beta2_2 + beta3_2*1.5*srs + beta4_2*2*rs

    c = 1 + 1/(2*A2*b)
    dcdrs = -(c-1)*dbdrs/b

    dum = log(c)
    dum1 = 1 + alpha1_2*rs

    G1 = -2*A2*dum1*dum
    dG1 = -2*A2*(alpha1_2*dum + dum1*dcdrs/c)

    ####################################
    # G(2)
    ####################################

    b = beta1_3*srs + rs*(beta2_3 + beta3_3*srs + beta4_3*rs)

    dbdrs = beta1_3*0.5/srs + beta2_3 + beta3_3*1.5*srs + beta4_3*2*rs

    c = 1 + 1/(2*A3*b)
    dcdrs = -(c-1)*dbdrs/b

    dum = log(c)
    dum1 = 1 + alpha1_3*rs

    G2 = -2*A3*dum1*dum
    dG2 = -2*A3*(alpha1_3*dum + dum1*dcdrs/c)

    ####################################
    # spin interpolation
    ####################################

    c = 1.92366105093154
    fpp0 = 1.70992093416137

    dum1 = 1 + zeta
    dum2 = 1 - zeta

    tmp1 = cbrt(dum1)
    tmp2 = cbrt(dum2)

    f = (tmp1^4 + tmp2^4 - 2)*c
    dfdz = (4/3)*(tmp1 - tmp2)*c

    dum1 = zeta^3
    dum = dum1*zeta

    Ec = G0 - G2*f/fpp0*(1-dum) + (G1-G0)*f*dum

    dEcdrs = dG0 - dG2*f/fpp0*(1-dum) + (dG1-dG0)*f*dum

    dEcdz = -G2/fpp0*(dfdz*(1-dum)-f*4*dum1) + (G1-G0)*(dfdz*dum + f*4*dum1)

    dum = dEcdrs*drsdd

    dEcdd1 = dum + dEcdz*dzdd1
    dEcdd2 = dum + dEcdz*dzdd2

    Vc1 = Ec + dtot*dEcdd1
    Vc2 = Ec + dtot*dEcdd2

    return Ec, Vc1, Vc2
end


############################################################
# LDA exchange
############################################################
function XC_EX(NSP, DS0, DS1, DS2)

    den_min = 1e-14
    den_min_half = 0.5e-14
    FTRD = 1.333333333333333
    TRD = 0.333333333333333
    TFTM = 0.519842099789746
    COE1 = 0.6203504908994
    COE2 = 0.610887057710857

    if NSP == 2
        D = DS1 + DS2
        if D <= den_min
            RS = COE1 * den_min^(-1/3)
            D0 = den_min_half
            D1 = den_min_half
            D = den_min
        else
            RS = COE1 * D^(-TRD)
            D0 = DS1
            D1 = DS2
        end
        Z = (D0-D1)/D

        FZ = ((1+Z)^FTRD + (1-Z)^FTRD - 2)/TFTM
        FZP = FTRD*(cbrt(1+Z) - cbrt(1-Z))/TFTM
    else
        D = max(DS0,den_min)
        RS = COE1 * D^(-TRD)

        Z = 0
        FZ = 0
        FZP = 0
    end

    VXP = -COE2/RS
    EXP = 0.75*VXP
    VXF = 0.111111111111111*VXP
    EXF = 0.111111111111111*EXP

    if NSP == 2
        VX1 = VXP + FZ*(VXF-VXP) + (1-Z)*FZP*(EXF-EXP)
        VX2 = VXP + FZ*(VXF-VXP) - (1+Z)*FZP*(EXF-EXP)

        EX1 = EXP + FZ*(EXF-EXP)
    else
        VX1 = VXP
        VX2 = 0.0
        EX1 = EXP
    end

    return EX1, VX1, VX2
end


function XC_PBE!(
    dens::Vector{Float64},
    GDENS::Vector{Float64},
    Exc::Vector{Float64},
    DEXDD::Vector{Float64},
    DECDD::Vector{Float64},
    DEXDGD::Vector{Float64},
    DECDGD::Vector{Float64})

    THD = 1/3
    den_min = 1e-14
    beta = 0.06672455060314922
    gamma = (1 - log(2))/(pi*pi)
    kappa = 0.8040
    mu = beta*pi*pi/3

    GD_min = 1e-14

    ############################################################
    # density
    ############################################################

    dens[1] = max(0.5*den_min, dens[1])
    dens[2] = max(0.5*den_min, dens[2])

    dt = dens[1] + dens[2]

    ############################################################
    # gradient
    ############################################################

    GD00 = GDENS[1]; GD10 = GDENS[2]; GD20 = GDENS[3]
    GD01 = GDENS[4]; GD11 = GDENS[5]; GD21 = GDENS[6]

    GDT0 = GD00 + GD01
    GDT1 = GD10 + GD11
    GDT2 = GD20 + GD21

    GDM0 = sqrt(GD00*GD00 + GD10*GD10 + GD20*GD20)
    GDM1 = sqrt(GD01*GD01 + GD11*GD11 + GD21*GD21)

    GDMT = sqrt(GDT0*GDT0 + GDT1*GDT1 + GDT2*GDT2)
    GDMT = max(GD_min,GDMT)

    ############################################################
    # Local correlation
    ############################################################

    Ec_unif, Vc1_unif, Vc2_unif = XC_PW92C(dens)

    ############################################################
    # correlation energy
    ############################################################

    rs = cbrt(3/(4*pi*dt))
    kF = cbrt(3*pi*pi*dt)
    ks = sqrt(4*kF/pi)

    zeta = (dens[1]-dens[2])/dt
    if zeta > 0.99
        zeta = 0.99
    end
    if zeta < -1.0
        zeta = -0.99
    end

    phi = 0.5*((1+zeta)^(2/3)+(1-zeta)^(2/3))

    t = GDMT/(2*phi*ks*dt)

    f1 = Ec_unif/(gamma*phi^3)
    f2 = exp(-f1)

    A = beta/gamma/(f2-1)

    f3 = t*t + A*t^4
    f4 = beta/gamma * f3/(1 + A*f3)

    H = gamma*phi^3*log(1+f4)

    Fc = Ec_unif + H

    ############################################################
    # derivatives
    ############################################################

    DKFDD = 1/3*kF/dt
    DKSDD = 0.5*ks*DKFDD/kF

    DZDD0 = 1/dt - zeta/dt
    DZDD1 = -1/dt - zeta/dt

    DPDZ = 0.5*2/3*(1/cbrt(1+zeta) - 1/cbrt(1-zeta))

    ############################################################
    # spin loop (IS=0)
    ############################################################

    DECUDD = (Vc1_unif-Ec_unif)/dt
    DPDD = DPDZ*DZDD0

    DTDD = (-t)*(DPDD/phi + DKSDD/ks + 1/dt)
    DF1DD = f1*(DECUDD/Ec_unif - 3*DPDD/phi)
    DF2DD = -f2*DF1DD
    DADD = -A*DF2DD/(f2-1)

    DF3DD = (2t + 4A*t^3)*DTDD + DADD*t^4

    DF4DD = f4*(DF3DD/f3 - (DADD*f3 + A*DF3DD)/(1 + A*f3))

    DHDD = 3H*DPDD/phi + gamma*phi^3*DF4DD/(1+f4)

    DFCDD0 = Vc1_unif + H + dt*DHDD

    DTDGD = (t/GDMT)*GDT0/GDMT
    DF3DGD = DTDGD*(2.0*t + 4.0*A*t*t*t)
    DF4DGD = f4*DF3DGD*(1.0/f3 - A/(1.0 + A*f3))
    DHDGD = gamma*phi*phi*phi*DF4DGD/(1.0 + f4)
    DFCDGD00 = dt*DHDGD
    DFCDGD01 = dt*DHDGD

    DTDGD = (t/GDMT)*GDT1/GDMT
    DF3DGD = DTDGD*(2.0*t + 4.0*A*t*t*t)
    DF4DGD = f4*DF3DGD*(1.0/f3 - A/(1.0 + A*f3))
    DHDGD = gamma*phi*phi*phi*DF4DGD/(1.0 + f4)
    DFCDGD10 = dt*DHDGD
    DFCDGD11 = dt*DHDGD

    DTDGD = (t/GDMT)*GDT2/GDMT
    DF3DGD = DTDGD*(2.0*t + 4.0*A*t*t*t)
    DF4DGD = f4*DF3DGD*(1.0/f3 - A/(1.0 + A*f3))
    DHDGD = gamma*phi*phi*phi*DF4DGD/(1.0 + f4)
    DFCDGD20 = dt*DHDGD
    DFCDGD21 = dt*DHDGD

    


    ############################################################
    # spin loop (IS=1)
    ############################################################

    DECUDD = (Vc2_unif-Ec_unif)/dt
    DPDD = DPDZ*DZDD1

    DTDD = (-t)*(DPDD/phi + DKSDD/ks + 1/dt)
    DF1DD = f1*(DECUDD/Ec_unif - 3*DPDD/phi)
    DF2DD = -f2*DF1DD
    DADD = -A*DF2DD/(f2-1)

    DF3DD = (2t + 4A*t^3)*DTDD + DADD*t^4
    DF4DD = f4*(DF3DD/f3 - (DADD*f3 + A*DF3DD)/(1 + A*f3))

    DHDD = 3H*DPDD/phi + gamma*phi^3*DF4DD/(1+f4)

    DFCDD1 = Vc2_unif + H + dt*DHDD

    ############################################################
    # exchange
    ############################################################

    Fx = 0.0

    DFXDD0=0.0
    DFXDD1=0.0

    DFXDGD00=0.0;DFXDGD10=0.0;DFXDGD20=0.0
    DFXDGD01=0.0;DFXDGD11=0.0;DFXDGD21=0.0

    ############################################################
    # spin loop (IS=0)
    ############################################################
    DS = max(den_min,2*dens[1])
    GDMS = max(GD_min,2*GDM0)

    KFS = cbrt(3*pi*pi*DS)

    s = GDMS/(2*KFS*DS)

    f1 = 1 + mu*s*s/kappa
    f = 1 + kappa - kappa/f1

    Ex_unif, Vx1_unif, Vx2_unif = XC_EX(1,DS,DS,DS)

    Fx += DS*Ex_unif*f

    DKFDD = THD*KFS/DS
    DSDD = s*(-(DKFDD/KFS) - 1/DS)

    DF1DD = 2*(f1-1)*DSDD/s
    DFDD = kappa*DF1DD/(f1*f1)

    DFXDD0 = Vx1_unif*f + DS*Ex_unif*DFDD

    # IX=0
    GDS = 2*GD00
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD00 = DS*Ex_unif*DFDGD

    # IX=1
    GDS = 2*GD10
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD10 = DS*Ex_unif*DFDGD

    # IX=2
    GDS = 2*GD20
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD20 = DS*Ex_unif*DFDGD



    ############################################################
    # spin loop (IS=1)
    ############################################################
    DS = max(den_min,2*dens[2])
    GDMS = max(GD_min,2*GDM1)

    KFS = cbrt(3*pi*pi*DS)

    s = GDMS/(2*KFS*DS)

    f1 = 1 + mu*s*s/kappa
    f = 1 + kappa - kappa/f1

    Ex_unif, Vx1_unif, Vx2_unif = XC_EX(1,DS,DS,DS)

    Fx += DS*Ex_unif*f

    DKFDD = THD*KFS/DS
    DSDD = s*(-(DKFDD/KFS) - 1/DS)

    DF1DD = 2*(f1-1)*DSDD/s
    DFDD = kappa*DF1DD/(f1*f1)

    DFXDD1 = Vx1_unif*f + DS*Ex_unif*DFDD

    # IX=0
    GDS = 2*GD01
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD01 = DS*Ex_unif*DFDGD

    # IX=1
    GDS = 2*GD11
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD11 = DS*Ex_unif*DFDGD

    # IX=2
    GDS = 2*GD21
    DSDGD = (s/GDMS)*GDS/GDMS
    DF1DGD = 2*mu*s*DSDGD/kappa
    DFDGD = kappa*DF1DGD/(f1*f1)
    DFXDGD21 = DS*Ex_unif*DFDGD


    Fx = 0.5*Fx/dt

    ############################################################
    # outputs
    ############################################################

    Exc[1] = Fx
    Exc[2] = Fc

    DEXDD[1] = DFXDD0
    DEXDD[2] = DFXDD1
    DECDD[1] = DFCDD0
    DECDD[2] = DFCDD1
    DEXDGD[1] = DFXDGD00
    DEXDGD[2] = DFXDGD10
    DEXDGD[3] = DFXDGD20
    DEXDGD[4] = DFXDGD01
    DEXDGD[5] = DFXDGD11
    DEXDGD[6] = DFXDGD21
    DECDGD[1] = DFCDGD00
    DECDGD[2] = DFCDGD10
    DECDGD[3] = DFCDGD20
    DECDGD[4] = DFCDGD01
    DECDGD[5] = DFCDGD11
    DECDGD[6] = DFCDGD21
end