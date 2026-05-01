function Calc_Ngrid(Ecut::Float64, Latvecs::Matrix{Float64})

    if Ecut < 0.0
        error("please check Ecut.")
    end

    CellVolume = abs(det(Latvecs)) # Bohr^3
    temp = CellVolume/pi * sqrt(Ecut)

    @views Ns1 = temp/norm(cross(Latvecs[2,:], Latvecs[3,:]))
    @views Ns2 = temp/norm(cross(Latvecs[3,:], Latvecs[1,:]))
    @views Ns3 = temp/norm(cross(Latvecs[1,:], Latvecs[2,:]))
    

    Ngrid1 = _Calc_Ngrid(Ns1)
    Ngrid2 = _Calc_Ngrid(Ns2)
    Ngrid3 = _Calc_Ngrid(Ns3)

    return (Ngrid1, Ngrid2, Ngrid3)
end


function _Calc_Ngrid(Ns::Float64)

    LgN = log(Ns)
    MinD = 1e+10

    popt = 1
    qopt = 1
    ropt = 1
    sopt = 1

    log2 = log(2)
    log3 = log(3)
    log5 = log(5)
    log7 = log(7)

    pmax = ceil(Int64, LgN/log2)
    for p = 0:pmax
        qmax = ceil(Int64, (LgN-p*log2)/log3)
        for q = 0:qmax
            rmax = ceil(Int64, (LgN-p*log2-q*log3)/log5)
            for r = 0:rmax
                smax = ceil(Int64, (LgN-p*log2-q*log3-r*log5)/log7)
                for s = 0:smax
                    
                    LgTN = p*log2+q*log3+r*log5+s*log7
                    
                    if abs(LgTN-LgN)<MinD
                        MinD = abs(LgTN-LgN)
                        popt = p
                        qopt = q
                        ropt = r
                        sopt = s
                    end
                end
            end
        end
    end

    k = 1
    for _ = 1:popt
        k *= 2
    end
    for _ = 1:qopt
        k *= 3
    end
    for _ = 1:ropt
        k *= 5
    end
    for _ = 1:sopt
        k *= 7
    end

    return k
end

