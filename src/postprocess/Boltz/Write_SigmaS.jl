function Write_SigmaS_2D(filename::String, cal_type::String, SpinPol::String, muE::Vector{Float64}, Temp, sigmas_mu)

    Nmu = length(muE)
    data = open(filename*"."*cal_type*"_sigmas"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# mu  sigmasxx(mu,up)  sigmasxy(mu,up)  sigmasyy(mu,up)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            2*sigmas_mu[imu,1,1], 2*sigmas_mu[imu,1,2], 2*sigmas_mu[imu,1,3])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  sigmasxx(mu,up)  sigmasxy(mu,up)  sigmasyy(mu,up)  sigmasxx(mu,dn)  sigmasxy(mu,dn)  sigmasyy(mu,dn)")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigmas_mu[imu,1,1], sigmas_mu[imu,1,2], sigmas_mu[imu,1,3],
                            sigmas_mu[imu,2,1], sigmas_mu[imu,2,2], sigmas_mu[imu,2,3])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  sigmasxx(mu)  sigmasxy(mu)  sigmasyy(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigmas_mu[imu,1,1], sigmas_mu[imu,1,2], sigmas_mu[imu,1,3])
        end
    end
end


function Write_SigmaS_3D(filename::String, cal_type::String, SpinPol::String, muE::Vector{Float64}, Temp, sigmas_mu)

    Nmu = length(muE)
    data = open(filename*"."*cal_type*"_sigmas"*string(Temp)*".dat", "w")

    @printf(data, "# %s", now())
    @printf(data, "# cal_type = %s\n", cal_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")
    
    if SpinPol == "off"
        @printf(data, "# mu  sigmasxx(mu,up)  sigmasxy(mu,up)  sigmasyy(mu,up)  sigmasxz(mu,up)  sigmasyz(mu,up)  sigmaszz(mu,up)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            2*sigmas_mu[imu,1,1], 2*sigmas_mu[imu,1,2], 2*sigmas_mu[imu,1,3],
                            2*sigmas_mu[imu,1,4], 2*sigmas_mu[imu,1,5], 2*sigmas_mu[imu,1,6])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  sigmasxx(mu,up)  sigmasxy(mu,up)  sigmasyy(mu,up)  sigmasxz(mu,up)  sigmasyz(mu,up)  sigmaszz(mu,up)  sigmasxx(mu,dn)  sigmasxy(mu,dn)  sigmasyy(mu,dn)  sigmasxz(mu,dn)  sigmasyz(mu,dn)  sigmaszz(mu,dn)")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigmas_mu[imu,1,1], sigmas_mu[imu,1,2], sigmas_mu[imu,1,3],
                            sigmas_mu[imu,1,4], sigmas_mu[imu,1,5], sigmas_mu[imu,1,6],
                            sigmas_mu[imu,2,1], sigmas_mu[imu,2,2], sigmas_mu[imu,2,3],
                            sigmas_mu[imu,2,4], sigmas_mu[imu,2,5], sigmas_mu[imu,2,6])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  sigmasxx(mu)  sigmasxy(mu)  sigmasyy(mu)  sigmasxz(mu)  sigmasyz(mu)  sigmaszz(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            sigmas_mu[imu,1,1], sigmas_mu[imu,1,2], sigmas_mu[imu,1,3],
                            sigmas_mu[imu,1,4], sigmas_mu[imu,1,5], sigmas_mu[imu,1,6])
        end
    end
end