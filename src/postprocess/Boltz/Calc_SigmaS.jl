function Calc_SigmaS(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    plane_type = boltz_setup.plane_type
    if plane_type
        Calc_SigmaS_2D(boltz_setup, TDF_Energy, TDF)
    else
        Calc_SigmaS_3D(boltz_setup, TDF_Energy, TDF)
    end
end


function Calc_SigmaS_2D(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    material = boltz_setup.material
    filename = boltz_setup.filename
    SpinPol = material.SpinPol
    TDF_dE = boltz_setup.TDF_dE
    TDF_Erange = boltz_setup.TDF_Erange
    Temp = boltz_setup.Temp
    NTemp = length(Temp)
    cal_type = "CWF"

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end


    TDF_Emin, TDF_Emax = TDF_Erange
    TDF_EneNum = floor(Int, (TDF_Emax-TDF_Emin)/TDF_dE)  # length of TDF
    dTDF = (TDF_Emax-TDF_Emin)/(TDF_EneNum-1)

    Nmu = TDF_EneNum
    muE = TDF_Energy


    kBT = zeros(Float64, NTemp)
    for iTemp = 1:NTemp
        kBT[iTemp] = Temp[iTemp] * k_B_SI/elem_charge_SI
    end
    fermi_T_dE = zeros(Float64, TDF_EneNum)


    SigmaS = zeros(Float64, Nmu, spinsize, 3)

    for iTemp = 1:NTemp

        fill!(SigmaS, 0.0)
        for spin = 1:spinsize, imu = 1:Nmu

            mu = muE[imu]
            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            sigmaS_Sumxx = 0.0
            sigmaS_Sumxy = 0.0
            sigmaS_Sumyy = 0.0

            @inbounds for ie = 1:TDF_EneNum
                sigmaS_Sumxx += fermi_T_dE[ie]*TDF[ie,spin,1]*(TDF_Energy[ie] - mu)
                sigmaS_Sumxy += fermi_T_dE[ie]*TDF[ie,spin,2]*(TDF_Energy[ie] - mu)
                sigmaS_Sumyy += fermi_T_dE[ie]*TDF[ie,spin,3]*(TDF_Energy[ie] - mu)
            end

            SigmaS[imu,spin,1] = sigmaS_Sumxx*dTDF/Temp[iTemp]
            SigmaS[imu,spin,2] = sigmaS_Sumxy*dTDF/Temp[iTemp]
            SigmaS[imu,spin,3] = sigmaS_Sumyy*dTDF/Temp[iTemp]
        end

        Write_SigmaS_2D(filename, cal_type, SpinPol, muE, Temp[iTemp], SigmaS)
    end
end


function Calc_SigmaS_3D(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    material = boltz_setup.material
    filename = boltz_setup.filename
    SpinPol = material.SpinPol
    TDF_dE = boltz_setup.TDF_dE
    TDF_Erange = boltz_setup.TDF_Erange
    Temp = boltz_setup.Temp
    NTemp = length(Temp)
    cal_type = "CWF"

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    end


    TDF_Emin, TDF_Emax = TDF_Erange
    TDF_EneNum = floor(Int, (TDF_Emax-TDF_Emin)/TDF_dE)  # length of TDF
    dTDF = (TDF_Emax-TDF_Emin)/(TDF_EneNum-1)

    Nmu = TDF_EneNum
    muE = TDF_Energy


    kBT = zeros(Float64, NTemp)
    for iTemp = 1:NTemp
        kBT[iTemp] = Temp[iTemp] * k_B_SI/elem_charge_SI
    end
    fermi_T_dE = zeros(Float64, TDF_EneNum)


    SigmaS = zeros(Float64, Nmu, spinsize, 6)

    for iTemp = 1:NTemp

        fill!(SigmaS, 0.0)
        for spin = 1:spinsize, imu = 1:Nmu

            mu = muE[imu]
            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            sigmaS_Sumxx = 0.0
            sigmaS_Sumxy = 0.0
            sigmaS_Sumxz = 0.0
            sigmaS_Sumyy = 0.0
            sigmaS_Sumyz = 0.0
            sigmaS_Sumzz = 0.0

            @inbounds for ie = 1:TDF_EneNum
                sigmaS_Sumxx += fermi_T_dE[ie]*TDF[ie,spin,1]*(TDF_Energy[ie] - mu)
                sigmaS_Sumxy += fermi_T_dE[ie]*TDF[ie,spin,2]*(TDF_Energy[ie] - mu)
                sigmaS_Sumyy += fermi_T_dE[ie]*TDF[ie,spin,3]*(TDF_Energy[ie] - mu)
                sigmaS_Sumxz += fermi_T_dE[ie]*TDF[ie,spin,4]*(TDF_Energy[ie] - mu)
                sigmaS_Sumyz += fermi_T_dE[ie]*TDF[ie,spin,5]*(TDF_Energy[ie] - mu)
                sigmaS_Sumzz += fermi_T_dE[ie]*TDF[ie,spin,6]*(TDF_Energy[ie] - mu)
            end

            SigmaS[imu,spin,1] = sigmaS_Sumxx*dTDF/Temp[iTemp]
            SigmaS[imu,spin,2] = sigmaS_Sumxy*dTDF/Temp[iTemp]
            SigmaS[imu,spin,3] = sigmaS_Sumyy*dTDF/Temp[iTemp]
            SigmaS[imu,spin,4] = sigmaS_Sumxz*dTDF/Temp[iTemp]
            SigmaS[imu,spin,5] = sigmaS_Sumyz*dTDF/Temp[iTemp]
            SigmaS[imu,spin,6] = sigmaS_Sumzz*dTDF/Temp[iTemp]
        end

        Write_SigmaS_3D(filename, cal_type, SpinPol, muE, Temp[iTemp], SigmaS)
    end
end