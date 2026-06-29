function Write_Seebeck_2D(filename::String, mat_type::String, SpinPol::String, muE::Vector{Float64}, Temp, seebeck_mu)

    Nmu = length(muE)
    data = open(filename*"_"*mat_type*"_seebeck"*string(Temp)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# mu  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up) \n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3], seebeck_mu[imu,1,4])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckxx(mu,dn)  seebeckxy(mu,dn)  seebeckyx(mu,dn)  seebeckyy(mu,dn)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3], seebeck_mu[imu,1,4],
                            seebeck_mu[imu,2,1], seebeck_mu[imu,2,2], seebeck_mu[imu,2,3], seebeck_mu[imu,2,4])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  seebeckxx(mu)  seebeckxy(mu)  seebeckyx(mu)  seebeckyy(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3], seebeck_mu[imu,1,4])
        end
    end
end


function Write_Seebeck_3D(filename::String, mat_type::String, SpinPol::String, muE::Vector{Float64}, Temp, seebeck_mu)

    Nmu = length(muE)
    data = open(filename*"_"*mat_type*"_seebeck"*string(Temp)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "# Nmu = %d\n", Nmu)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# mu  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckxz(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckyz(mu,up)  seebeckzx(mu,up)  seebeckzy(mu,up)  seebeckzz(mu,up)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3],
                            seebeck_mu[imu,1,4], seebeck_mu[imu,1,5], seebeck_mu[imu,1,6],
                            seebeck_mu[imu,1,7], seebeck_mu[imu,1,8], seebeck_mu[imu,1,9])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckxz(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckyz(mu,up)  seebeckzx(mu,up)  seebeckzy(mu,up)  seebeckzz(mu,up)  seebeckxx(mu,dn)  seebeckxy(mu,dn)  seebeckxz(mu,dn)  seebeckyx(mu,dn)  seebeckyy(mu,dn)  seebeckyz(mu,dn)  seebeckzx(mu,dn)  seebeckzy(mu,dn)  seebeckzz(mu,dn)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3],
                            seebeck_mu[imu,1,4], seebeck_mu[imu,1,5], seebeck_mu[imu,1,6],
                            seebeck_mu[imu,1,7], seebeck_mu[imu,1,8], seebeck_mu[imu,1,9],
                            seebeck_mu[imu,2,1], seebeck_mu[imu,2,2], seebeck_mu[imu,2,3],
                            seebeck_mu[imu,2,4], seebeck_mu[imu,2,5], seebeck_mu[imu,2,6],
                            seebeck_mu[imu,2,7], seebeck_mu[imu,2,8], seebeck_mu[imu,2,9])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu  seebeckxx(mu)  seebeckxy(mu)  seebeckxz(mu)  seebeckyx(mu)  seebeckyy(mu)  seebeckyz(mu)  seebeckzx(mu)  seebeckzy(mu)  seebeckzz(mu)\n")
        for imu = 1:Nmu
            @printf(data, "%15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            muE[imu], 
                            seebeck_mu[imu,1,1], seebeck_mu[imu,1,2], seebeck_mu[imu,1,3],
                            seebeck_mu[imu,1,4], seebeck_mu[imu,1,5], seebeck_mu[imu,1,6],
                            seebeck_mu[imu,1,7], seebeck_mu[imu,1,8], seebeck_mu[imu,1,9])
        end
    end
end


function Write_Seebeck_decomp_2D(filename::String, mat_type::String, SpinPol::String, Nwann, mu, Temp, seebeck_mu)

    data = open(filename*"_"*mat_type*"_seebeck_decomp_"*string(mu)*"_"*string(Temp)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Nwann = %d\n", Nwann)
    @printf(data, "# mu = %5.15f\n", mu)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# wann  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up) \n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist,
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3], seebeck_mu[ist,1,4])
        end
    elseif SpinPol == "on"
        @printf(data, "# wann  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckxx(mu,dn)  seebeckxy(mu,dn)  seebeckyx(mu,dn)  seebeckyy(mu,dn)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3], seebeck_mu[ist,1,4],
                            seebeck_mu[ist,2,1], seebeck_mu[ist,2,2], seebeck_mu[ist,2,3], seebeck_mu[ist,2,4])
        end
    elseif SpinPol == "nc"
        @printf(data, "# wann  seebeckxx(mu)  seebeckxy(mu)  seebeckyx(mu)  seebeckyy(mu)\n")
        for ist = 1:Nwann
            @printf(data, "%3d  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3], seebeck_mu[ist,1,4])
        end
    end
end


function Write_Seebeck_decomp_3D(filename::String, mat_type::String, SpinPol::String, Nwann, mu, Temp, seebeck_mu)

    data = open(filename*"_"*mat_type*"_seebeck_decomp_"*string(mu)*"_"*string(Temp)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# Nwann = %d\n", Nwann)
    @printf(data, "# mu = %5.15f\n", mu)
    @printf(data, "# Temp = %5.15f\n", Temp)
    @printf(data, "\n")

    if SpinPol == "off"
        @printf(data, "# wann  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckxz(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckyz(mu,up)  seebeckzx(mu,up)  seebeckzy(mu,up)  seebeckzz(mu,up)\n")
        for ist = 1:Nwann
            @printf(data, "%d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3],
                            seebeck_mu[ist,1,4], seebeck_mu[ist,1,5], seebeck_mu[ist,1,6],
                            seebeck_mu[ist,1,7], seebeck_mu[ist,1,8], seebeck_mu[ist,1,9])
        end
    elseif SpinPol == "on"
        @printf(data, "# wann  seebeckxx(mu,up)  seebeckxy(mu,up)  seebeckxz(mu,up)  seebeckyx(mu,up)  seebeckyy(mu,up)  seebeckyz(mu,up)  seebeckzx(mu,up)  seebeckzy(mu,up)  seebeckzz(mu,up)  seebeckxx(mu,dn)  seebeckxy(mu,dn)  seebeckxz(mu,dn)  seebeckyx(mu,dn)  seebeckyy(mu,dn)  seebeckyz(mu,dn)  seebeckzx(mu,dn)  seebeckzy(mu,dn)  seebeckzz(mu,dn)\n")
        for ist = 1:Nwann
            @printf(data, "%d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3],
                            seebeck_mu[ist,1,4], seebeck_mu[ist,1,5], seebeck_mu[ist,1,6],
                            seebeck_mu[ist,1,7], seebeck_mu[ist,1,8], seebeck_mu[ist,1,9],
                            seebeck_mu[ist,2,1], seebeck_mu[ist,2,2], seebeck_mu[ist,2,3],
                            seebeck_mu[ist,2,4], seebeck_mu[ist,2,5], seebeck_mu[ist,2,6],
                            seebeck_mu[ist,2,7], seebeck_mu[ist,2,8], seebeck_mu[ist,2,9])
        end
    elseif SpinPol == "nc"
        @printf(data, "# wann  seebeckxx(mu)  seebeckxy(mu)  seebeckxz(mu)  seebeckyx(mu)  seebeckyy(mu)  seebeckyz(mu)  seebeckzx(mu)  seebeckzy(mu)  seebeckzz(mu)\n")
        for ist = 1:Nwann
            @printf(data, "%d  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f  %15.12f\n", 
                            ist, 
                            seebeck_mu[ist,1,1], seebeck_mu[ist,1,2], seebeck_mu[ist,1,3],
                            seebeck_mu[ist,1,4], seebeck_mu[ist,1,5], seebeck_mu[ist,1,6],
                            seebeck_mu[ist,1,7], seebeck_mu[ist,1,8], seebeck_mu[ist,1,9])
        end
    end
end