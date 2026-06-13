function Write_Sigma_2D(filename::String, cal_type::String, SpinPol::String, muE::Vector{Float64}, Temp, sigma_mu)

    Nmu = length(muE)
    data = open(filename*"."*cal_type*"_sigma"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# mu  sigmaxx(mu,up)  sigmaxy(mu,up)  sigmayy(mu,up)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            2*sigma_mu[imu,1,1], 2*sigma_mu[imu,1,2], 2*sigma_mu[imu,1,3])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  sigmaxx(mu,up)  sigmaxy(mu,up)  sigmayy(mu,up)   sigmaxx(mu,dn)  sigmaxy(mu,dn)  sigmayy(mu,dn)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigma_mu[imu,1,1], sigma_mu[imu,1,2], sigma_mu[imu,1,3],
                            sigma_mu[imu,2,1], sigma_mu[imu,2,2], sigma_mu[imu,2,3])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  sigmaxx(mu)  sigmaxy(mu)  sigmayy(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigma_mu[imu,1,1], sigma_mu[imu,1,2], sigma_mu[imu,1,3])
        end
    end
end


function Write_Sigma_3D(filename::String, cal_type::String, SpinPol::String, muE::Vector{Float64}, Temp, sigma_mu)

    Nmu = length(muE)
    data = open(filename*"."*cal_type*"_sigma"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# mu  sigmaxx(mu,up)  sigmaxy(mu,up)  sigmayy(mu,up)   sigmaxz(mu,up)  sigmayz(mu,up)  sigmazz(mu,up)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            2*sigma_mu[imu,1,1], 2*sigma_mu[imu,1,2], 2*sigma_mu[imu,1,3],
                            2*sigma_mu[imu,1,4], 2*sigma_mu[imu,1,5], 2*sigma_mu[imu,1,6])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  sigmaxx(mu,up)  sigmaxy(mu,up)  sigmayy(mu,up)   sigmaxz(mu,up)  sigmayz(mu,up)  sigmazz(mu,up)  sigmaxx(mu,dn)  sigmaxy(mu,dn)  sigmayy(mu,dn)   sigmaxz(mu,dn)  sigmayz(mu,dn)  sigmazz(mu,dn)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigma_mu[imu,1,1], sigma_mu[imu,1,2], sigma_mu[imu,1,3],
                            sigma_mu[imu,1,4], sigma_mu[imu,1,5], sigma_mu[imu,1,6],
                            sigma_mu[imu,2,1], sigma_mu[imu,2,2], sigma_mu[imu,2,3],
                            sigma_mu[imu,2,4], sigma_mu[imu,2,5], sigma_mu[imu,2,6])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  sigmaxx(mu)  sigmaxy(mu)  sigmayy(mu)   sigmaxz(mu)  sigmayz(mu)  sigmazz(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigma_mu[imu,1,1], sigma_mu[imu,1,2], sigma_mu[imu,1,3],
                            sigma_mu[imu,1,4], sigma_mu[imu,1,5], sigma_mu[imu,1,6])
        end
    end
end


function Write_Sigma_decomp_2D(filename::String, cal_type::String, SpinPol::String, Nwann, mu, Temp, sigma_mu)

    data = open(filename*"."*cal_type*"_sigma_decomp_"*string(mu)*"_"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Nwann = %d\n", Nwann)
    @printf(data, "# mu = %5.15f\n", mu)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# wann  sigmaxx(wann,up)  sigmaxy(wann,up)  sigmayy(wann,up)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            2*sigma_mu[ist,1,1], 2*sigma_mu[ist,1,2], 2*sigma_mu[ist,1,3])
        end
    elseif SpinPol == "on"
        @printf(data, "# wann  sigmaxx(wann,up)  sigmaxy(wann,up)  sigmayy(wann,up)  sigmaxx(wann,dn)  sigmaxy(wann,dn)  sigmayy(wann,dn)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            sigma_mu[ist,1,1], sigma_mu[ist,1,2], sigma_mu[ist,1,3],
                            sigma_mu[ist,2,1], sigma_mu[ist,2,2], sigma_mu[ist,2,3])
        end
    elseif SpinPol == "nc"
        @printf(data, "# wann  sigmaxx(wann)  sigmaxy(wann)  sigmayy(wann)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            sigma_mu[ist,1,1], sigma_mu[ist,1,2], sigma_mu[ist,1,3])
        end
    end
end


function Write_Sigma_decomp_3D(filename::String, cal_type::String, SpinPol::String, Nwann, mu, Temp, sigma_mu)

    data = open(filename*"."*cal_type*"_sigma_decomp_"*string(mu)*"_"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Nwann = %d\n", Nwann)
    @printf(data, "# mu = %5.15f\n", mu)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# wann  sigmaxx(wann,up)  sigmaxy(wann,up)  sigmayy(wann,up)   sigmaxz(wann,up)  sigmayz(wann,up)  sigmazz(wann,up)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            2*sigma_mu[ist,1,1], 2*sigma_mu[ist,1,2], 2*sigma_mu[ist,1,3],
                            2*sigma_mu[ist,1,4], 2*sigma_mu[ist,1,5], 2*sigma_mu[ist,1,6])
        end
    elseif SpinPol == "on"
        @printf(data, "# wann  sigmaxx(wann,up)  sigmaxy(wann,up)  sigmayy(wann,up)   sigmaxz(wann,up)  sigmayz(wann,up)  sigmazz(wann,up)  sigmaxx(wann,dn)  sigmaxy(wann,dn)  sigmayy(wann,dn)   sigmaxz(wann,dn)  sigmayz(wann,dn)  sigmazz(wann,dn)\n")
        for ist = 1:Nwann
            @printf(data, "%3d   %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            sigma_mu[ist,1,1], sigma_mu[ist,1,2], sigma_mu[ist,1,3],
                            sigma_mu[ist,1,4], sigma_mu[ist,1,5], sigma_mu[ist,1,6],
                            sigma_mu[ist,2,1], sigma_mu[ist,2,2], sigma_mu[ist,2,3],
                            sigma_mu[ist,2,4], sigma_mu[ist,2,5], sigma_mu[ist,2,6])
        end
    elseif SpinPol == "nc"
        @printf(data, "# wann  sigmaxx(wann)  sigmaxy(wann)  sigmayy(wann)   sigmaxz(wann)  sigmayz(wann)  sigmazz(wann)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            sigma_mu[ist,1,1], sigma_mu[ist,1,2], sigma_mu[ist,1,3],
                            sigma_mu[ist,1,4], sigma_mu[ist,1,5], sigma_mu[ist,1,6])
        end
    end
end