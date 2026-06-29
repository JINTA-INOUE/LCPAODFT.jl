function Calc_Seebeck_decomp(boltz_setup::Boltz_Setup, TDF_Energy, TDF_decomp)

    plane_type = boltz_setup.plane_type
    if plane_type
        Calc_Seebeck_decomp_2D(boltz_setup, TDF_Energy, TDF_decomp)
    else
        Calc_Seebeck_decomp_3D(boltz_setup, TDF_Energy, TDF_decomp)
    end
end


function Calc_Seebeck_decomp_2D(boltz_setup::Boltz_Setup, TDF_Energy, TDF_decomp)

    material = boltz_setup.material
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    filename = boltz_setup.filename
    mat_type = boltz_setup.mat_type
    TDF_Erange = boltz_setup.TDF_Erange
    muE = boltz_setup.muE
    Temp = boltz_setup.Temp
    NTemp = length(Temp)

    TDF_EneNum = length(TDF_Energy)
    dTDF = (TDF_Erange[2]-TDF_Erange[1])/(TDF_EneNum-1)
    Nmu = length(muE)


    kBT = zeros(Float64, NTemp)
    for iTemp = 1:NTemp
        kBT[iTemp] = Temp[iTemp] * k_B_SI/elem_charge_SI
    end
    fermi_T_dE = zeros(Float64, TDF_EneNum)



    TDF = zeros(Float64,TDF_EneNum,spinsize,3)
    for spin = 1:spinsize, ie = 1:TDF_EneNum 

        TDF_Sumxx = 0.0
        TDF_Sumxy = 0.0
        TDF_Sumyy = 0.0
        for ist = 1:Nwann
            TDF_Sumxx += TDF_decomp[ist][ie,spin,1]
            TDF_Sumxy += TDF_decomp[ist][ie,spin,2]
            TDF_Sumyy += TDF_decomp[ist][ie,spin,3]
        end

        TDF[ie,spin,1] = TDF_Sumxx
        TDF[ie,spin,2] = TDF_Sumxy
        TDF[ie,spin,3] = TDF_Sumyy
    end



    sigma22 = zeros(Float64, 2, 2)
    inv_sigma22 = zeros(Float64, 2, 2)
    sigmaS22 = zeros(Float64, 2, 2)
    Seebeck22 = zeros(Float64, 2, 2)
    Seebeck = zeros(Float64, Nwann, spinsize, 4)

    for iTemp = 1:NTemp, imu = 1:Nmu

        mu = muE[imu]
        fill!(Seebeck, 0.0)

        for spin = 1:spinsize, ist = 1:Nwann
            
            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            sigma_Sumxx = 0.0
            sigma_Sumxy = 0.0
            sigma_Sumyy = 0.0

            sigmaS_Sumxx = 0.0
            sigmaS_Sumxy = 0.0
            sigmaS_Sumyy = 0.0

            @inbounds for ie = 1:TDF_EneNum
                sigma_Sumxx += fermi_T_dE[ie] * TDF[ie,spin,1]  # xx
                sigma_Sumxy += fermi_T_dE[ie] * TDF[ie,spin,2]  # xy
                sigma_Sumyy += fermi_T_dE[ie] * TDF[ie,spin,3]  # yy

                sigmaS_Sumxx += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,1] * (TDF_Energy[ie] - mu)
                sigmaS_Sumxy += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,2] * (TDF_Energy[ie] - mu)
                sigmaS_Sumyy += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,3] * (TDF_Energy[ie] - mu)
            end


            sigma22[1,1] = sigma_Sumxx*dTDF
            sigma22[1,2] = sigma_Sumxy*dTDF
            sigma22[2,2] = sigma_Sumyy*dTDF
            sigma22[2,1] = sigma22[1,2]

            sigmaS22[1,1] = sigmaS_Sumxx*dTDF/Temp[iTemp]
            sigmaS22[1,2] = sigmaS_Sumxy*dTDF/Temp[iTemp]
            sigmaS22[2,2] = sigmaS_Sumyy*dTDF/Temp[iTemp]
            sigmaS22[2,1] = sigmaS22[1,2]

            determinant = Seebeck_inv2!(sigma22, inv_sigma22)
            if iszero(abs(determinant))
                fill!(inv_sigma22, 0.0)
            else
                @. inv_sigma22 = inv_sigma22/determinant
            end

            LinearAlgebra.matmul2x2!(Seebeck22, 'N', 'N', inv_sigma22, sigmaS22)
            
            Seebeck[ist,spin,1] = -Seebeck22[1,1]    # xx
            Seebeck[ist,spin,2] = -Seebeck22[1,2]    # xy
            Seebeck[ist,spin,3] = -Seebeck22[2,1]    # yx
            Seebeck[ist,spin,4] = -Seebeck22[2,2]    # yy
        end

        Write_Seebeck_decomp_2D(filename, mat_type, SpinPol, Nwann, mu, Temp[iTemp], Seebeck)
    end
end


function Calc_Seebeck_decomp_3D(boltz_setup::Boltz_Setup, TDF_Energy, TDF_decomp)

    material = boltz_setup.material
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    filename = boltz_setup.filename
    mat_type = boltz_setup.mat_type
    TDF_Erange = boltz_setup.TDF_Erange
    muE = boltz_setup.muE
    Temp = boltz_setup.Temp
    NTemp = length(Temp)

    TDF_EneNum = length(TDF_Energy)
    dTDF = (TDF_Erange[2]-TDF_Erange[1])/(TDF_EneNum-1)
    Nmu = length(muE)


    kBT = zeros(Float64, NTemp)
    for iTemp = 1:NTemp
        kBT[iTemp] = Temp[iTemp] * k_B_SI/elem_charge_SI
    end
    fermi_T_dE = zeros(Float64, TDF_EneNum)


    TDF = zeros(Float64,TDF_EneNum,spinsize,6)
    for spin = 1:spinsize, ie = 1:TDF_EneNum 

        TDF_Sumxx = 0.0
        TDF_Sumxy = 0.0
        TDF_Sumyy = 0.0
        TDF_Sumxz = 0.0
        TDF_Sumyz = 0.0
        TDF_Sumzz = 0.0
        for ist = 1:Nwann
            TDF_Sumxx += TDF_decomp[ist][ie,spin,1]
            TDF_Sumxy += TDF_decomp[ist][ie,spin,2]
            TDF_Sumyy += TDF_decomp[ist][ie,spin,3]
            TDF_Sumxz += TDF_decomp[ist][ie,spin,4]
            TDF_Sumyz += TDF_decomp[ist][ie,spin,5]
            TDF_Sumzz += TDF_decomp[ist][ie,spin,6]
        end

        TDF[ie,spin,1] = TDF_Sumxx
        TDF[ie,spin,2] = TDF_Sumxy
        TDF[ie,spin,3] = TDF_Sumyy
        TDF[ie,spin,4] = TDF_Sumxz
        TDF[ie,spin,5] = TDF_Sumyz
        TDF[ie,spin,6] = TDF_Sumzz
    end


    sigma33 = zeros(Float64, 3, 3)
    inv_sigma33 = zeros(Float64, 3, 3)
    sigmaS33 = zeros(Float64, 3, 3)
    Seebeck33 = zeros(Float64, 3, 3)
    Seebeck = zeros(Float64, Nwann, spinsize, 9)

    for iTemp = 1:NTemp, imu = 1:Nmu

        mu = muE[imu]
        fill!(Seebeck, 0.0)
        for spin = 1:spinsize, ist = 1:Nwann

            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            sigma_Sumxx = 0.0
            sigma_Sumxy = 0.0
            sigma_Sumxz = 0.0
            sigma_Sumyy = 0.0
            sigma_Sumyz = 0.0
            sigma_Sumzz = 0.0

            sigmaS_Sumxx = 0.0
            sigmaS_Sumxy = 0.0
            sigmaS_Sumxz = 0.0
            sigmaS_Sumyy = 0.0
            sigmaS_Sumyz = 0.0
            sigmaS_Sumzz = 0.0

            @inbounds for ie = 1:TDF_EneNum
                sigma_Sumxx += fermi_T_dE[ie] * TDF[ie,spin,1]
                sigma_Sumxy += fermi_T_dE[ie] * TDF[ie,spin,2]
                sigma_Sumyy += fermi_T_dE[ie] * TDF[ie,spin,3]
                sigma_Sumxz += fermi_T_dE[ie] * TDF[ie,spin,4]
                sigma_Sumyz += fermi_T_dE[ie] * TDF[ie,spin,5]
                sigma_Sumzz += fermi_T_dE[ie] * TDF[ie,spin,6]

                sigmaS_Sumxx += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,1] * (TDF_Energy[ie] - mu)
                sigmaS_Sumxy += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,2] * (TDF_Energy[ie] - mu)
                sigmaS_Sumyy += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,3] * (TDF_Energy[ie] - mu)
                sigmaS_Sumxz += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,4] * (TDF_Energy[ie] - mu)
                sigmaS_Sumyz += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,5] * (TDF_Energy[ie] - mu)
                sigmaS_Sumzz += fermi_T_dE[ie] * TDF_decomp[ist][ie,spin,6] * (TDF_Energy[ie] - mu)
            end


            sigma33[1,1] = sigma_Sumxx*dTDF
            sigma33[1,2] = sigma_Sumxy*dTDF
            sigma33[2,2] = sigma_Sumyy*dTDF
            sigma33[1,3] = sigma_Sumxz*dTDF
            sigma33[2,3] = sigma_Sumyz*dTDF
            sigma33[3,3] = sigma_Sumzz*dTDF
            sigma33[2,1] = sigma33[1,2]
            sigma33[3,1] = sigma33[1,3]
            sigma33[3,2] = sigma33[2,3]

            sigmaS33[1,1] = sigmaS_Sumxx*dTDF/Temp[iTemp]
            sigmaS33[1,2] = sigmaS_Sumxy*dTDF/Temp[iTemp]
            sigmaS33[2,2] = sigmaS_Sumyy*dTDF/Temp[iTemp]
            sigmaS33[1,3] = sigmaS_Sumxz*dTDF/Temp[iTemp]
            sigmaS33[2,3] = sigmaS_Sumyz*dTDF/Temp[iTemp]
            sigmaS33[3,3] = sigmaS_Sumzz*dTDF/Temp[iTemp]
            sigmaS33[2,1] = sigmaS33[1,2]
            sigmaS33[3,1] = sigmaS33[1,3]
            sigmaS33[3,2] = sigmaS33[2,3]

            determinant = Seebeck_inv3!(sigma33, inv_sigma33)
            if iszero(abs(determinant))
                fill!(inv_sigma33, 0.0)
            else
                @. inv_sigma33 = inv_sigma33/determinant
            end

            LinearAlgebra.matmul3x3!(Seebeck33, 'N', 'N', inv_sigma33, sigmaS33)
            
            Seebeck[ist,spin,1] = -Seebeck33[1,1]    # xx
            Seebeck[ist,spin,2] = -Seebeck33[1,2]    # xy
            Seebeck[ist,spin,3] = -Seebeck33[1,3]    # xz
            Seebeck[ist,spin,4] = -Seebeck33[2,1]    # yx
            Seebeck[ist,spin,5] = -Seebeck33[2,2]    # yy
            Seebeck[ist,spin,6] = -Seebeck33[2,3]    # yz
            Seebeck[ist,spin,7] = -Seebeck33[3,1]    # zx
            Seebeck[ist,spin,8] = -Seebeck33[3,2]    # zy
            Seebeck[ist,spin,9] = -Seebeck33[3,3]    # zz
        end

        Write_Seebeck_decomp_3D(filename, mat_type, SpinPol, Nwann, mu, Temp[iTemp], Seebeck)
    end
end

