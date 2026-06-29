function Calc_Sigma(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    plane_type = boltz_setup.plane_type
    if plane_type
        Calc_Sigma_2D(boltz_setup, TDF_Energy, TDF)
    else
        Calc_Sigma_3D(boltz_setup, TDF_Energy, TDF)
    end
end


function Calc_Sigma_2D(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    material = boltz_setup.material
    mat_type = boltz_setup.mat_type
    filename = boltz_setup.filename
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    TDF_Erange = boltz_setup.TDF_Erange
    Temp = boltz_setup.Temp
    NTemp = length(Temp)


    TDF_EneNum = length(TDF_Energy)
    dTDF = (TDF_Erange[2]-TDF_Erange[1])/(TDF_EneNum-1)
    Nmu = TDF_EneNum
    muE = TDF_Energy


    kBT = zeros(Float64, NTemp)
    for iTemp = 1:NTemp
        kBT[iTemp] = Temp[iTemp] * k_B_SI/elem_charge_SI
    end
    fermi_T_dE = zeros(Float64, TDF_EneNum)


    sigma = zeros(Float64, Nmu, spinsize, 3)

    for iTemp = 1:NTemp
        
        fill!(sigma, 0.0)
        for spin = 1:spinsize, imu = 1:Nmu
            mu = muE[imu]
            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            Sumxx = 0.0
            Sumxy = 0.0
            Sumyy = 0.0
            for ie = 1:TDF_EneNum
                Sumxx += fermi_T_dE[ie]*TDF[ie,spin,1]        # xx direction
                Sumxy += fermi_T_dE[ie]*TDF[ie,spin,2]        # xy direction
                Sumyy += fermi_T_dE[ie]*TDF[ie,spin,3]        # yy direction
            end

            sigma[imu,spin,1] = Sumxx*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,2] = Sumxy*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,3] = Sumyy*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
        end

        Write_Sigma_2D(filename, mat_type, SpinPol, muE, Temp[iTemp], sigma)
    end    
end


function Calc_Sigma_3D(boltz_setup::Boltz_Setup, TDF_Energy, TDF)

    material = boltz_setup.material
    mat_type = boltz_setup.mat_type
    filename = boltz_setup.filename
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    TDF_dE = boltz_setup.TDF_dE
    TDF_Erange = boltz_setup.TDF_Erange
    Temp = boltz_setup.Temp
    NTemp = length(Temp)
    

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


    sigma = zeros(Float64, Nmu, spinsize, 6)

    for iTemp = 1:NTemp
        
        fill!(sigma, 0.0)
        for spin = 1:spinsize, imu = 1:Nmu
            mu = muE[imu]
            Set_dfdE!(fermi_T_dE, mu, TDF_Energy, TDF_EneNum, kBT[iTemp])

            Sumxx = 0.0
            Sumxy = 0.0
            Sumyy = 0.0
            Sumxz = 0.0
            Sumyz = 0.0
            Sumzz = 0.0

            for ie = 1:TDF_EneNum
                Sumxx += fermi_T_dE[ie] * TDF[ie,spin,1]        # xx direction
                Sumxy += fermi_T_dE[ie] * TDF[ie,spin,2]        # xy direction
                Sumyy += fermi_T_dE[ie] * TDF[ie,spin,3]        # yy direction
                Sumxz += fermi_T_dE[ie] * TDF[ie,spin,4]        # xz direction
                Sumyz += fermi_T_dE[ie] * TDF[ie,spin,5]        # yz direction
                Sumzz += fermi_T_dE[ie] * TDF[ie,spin,6]        # zz direction
            end

            sigma[imu,spin,1] = Sumxx*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,2] = Sumxy*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,3] = Sumyy*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,4] = Sumxz*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,5] = Sumyz*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
            sigma[imu,spin,6] = Sumzz*dTDF*elem_charge_SI^3/(hbar_SI^2)*1.0e-5
        end

        Write_Sigma_3D(filename, mat_type, SpinPol, muE, Temp[iTemp], sigma)
    end
end