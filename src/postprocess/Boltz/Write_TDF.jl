function Write_TDF_3elements(filename::String, mat_type::String, SpinPol::String, muE::Vector{Float64}, tau, TDF)

    Nmu = length(muE)
    data = open(filename*"_"*mat_type*"_TDF"*string(tau)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# tau = %5.15f\n", tau)
    @printf(data, "# Nmu = %d\n", Nmu)

    if SpinPol == "off"
        @printf(data, "# mu TDF_xx(mu,up)  TDF_xy(mu,up)  TDF_yy(mu,up)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,1,2], TDF[ie,1,3])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu TDF_xx(mu,up)  TDF_xx(mu,dn)  TDF_xy(mu,up)  TDF_xy(mu,dn)  TDF_yy(mu,up)  TDF_yy(mu,dn)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,2,1], TDF[ie,1,2], 
                            TDF[ie,2,2], TDF[ie,1,3], TDF[ie,2,3])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu TDF_xx(mu) TDF_xy(mu) TDF_yy(mu)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,1,2], TDF[ie,1,3])
        end
    end

    close(data)
end


function Write_TDF_6elements(filename::String, mat_type::String, SpinPol::String, muE::Vector{Float64}, tau, TDF)
    
    Nmu = length(muE)
    data = open(filename*"_"*mat_type*"_TDF"*string(tau)*".dat", "w")

    @printf(data, "# %s\n", now())
    @printf(data, "# mat_type = %s\n", mat_type)
    @printf(data, "# SpinPol = %s\n", SpinPol)
    @printf(data, "# tau = %5.15f\n", tau)
    @printf(data, "# Nmu = %d\n", Nmu)

    if SpinPol == "off"
        @printf(data, "# mu TDF_xx(mu,up)  TDF_xy(mu,up)  TDF_yy(mu,up)  TDF_xz(mu,up)  TDF_yz(mu,up) TDF_zz(mu,up)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,1,2], TDF[ie,1,3], 
                            TDF[ie,1,4], TDF[ie,1,5], TDF[ie,1,6])
        end
    elseif SpinPol == "on"
        @printf(data, "# mu TDF_xx(mu,up)  TDF_xx(mu,dn)  TDF_xy(mu,up)  TDF_xy(mu,dn)  TDF_yy(mu,up)  TDF_yy(mu,dn)  TDF_xz(mu,up)  TDF_xz(mu,dn)  TDF_yz(mu,up)  TDF_yz(mu,dn)  TDF_zz(mu,up)  TDF_zz(mu,dn)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,2,1], TDF[ie,1,2], TDF[ie,2,2], 
                            TDF[ie,1,3], TDF[ie,2,3], TDF[ie,1,4], TDF[ie,2,4],
                            TDF[ie,1,6], TDF[ie,2,6], TDF[ie,1,6], TDF[ie,2,6])
        end
    elseif SpinPol == "nc"
        @printf(data, "# mu TDF_xx(mu)  TDF_xy(mu)  TDF_yy(mu)  TDF_xz(mu)  TDF_yz(mu) TDF_zz(mu)\n")
        for ie = 1:Nmu
            @printf(data, "%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n", 
                            muE[ie], 
                            TDF[ie,1,1], TDF[ie,1,2], TDF[ie,1,3],
                            TDF[ie,1,4], TDF[ie,1,5], TDF[ie,1,6])
        end
    end

    close(data)
end