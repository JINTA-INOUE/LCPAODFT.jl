function Ylm_table(l::Int, m::Int, theta, phi)

    if l < abs(m)
        error("please check l, m.")
    end

    if l == 0

        return 0.5*sqrt(1.0/pi) + im*0.0

    elseif l == 1

        if m == 0
            return 0.5*sqrt(3.0/pi)*cos(theta) + im*0.0
        else    # if m == -1 || m == 1
            return sign(-m)*0.5*sqrt(3.0/2.0/pi)*exp(im*m*phi)*sin(theta)
        end

    elseif l == 2

        if m == 0
            return 0.25*sqrt(5.0/pi)*(3*cos(theta)^2 - 1) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*0.5*sqrt(15.0/2.0/pi)*exp(im*m*phi)*sin(theta)*cos(theta)
        else    # if m == -2 || m == 2
            return 0.25*sqrt(15.0/2.0/pi)*exp(im*m*phi)*sin(theta)^2
        end
        
    elseif l == 3

        if m == 0
            return 1/4*sqrt(7/pi)*(5*cos(theta)^3 - 3*cos(theta)) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*1/8*sqrt(21/pi)*exp(im*m*phi)*sin(theta)*(5*cos(theta)^2 - 1)
        elseif m == -2 || m == 2
            return 1/4*sqrt(105/2/pi)*exp(im*m*phi)*sin(theta)^2*cos(theta)
        else    # if m == -3 || m == 3
            return sign(-m)*1/8*sqrt(35/pi)*exp(im*m*phi)*sin(theta)^3
        end

    elseif l == 4

        if m == 0
            return 3/16*sqrt(1/pi)*(35*cos(theta)^4 - 30*cos(theta)^2 + 3) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*3/8*sqrt(5/pi)*exp(im*m*phi)*sin(theta)*(7*cos(theta)^3 - 3*cos(theta))
        elseif m == -2 || m == 2
            return 3/8*sqrt(5/pi/2)*exp(im*m*phi)*sin(theta)^2*(7*cos(theta)^2 - 1)
        elseif m == -3 || m == 3
            return sign(-m)*3/8*sqrt(35/pi)*exp(im*m*phi)*sin(theta)^3*cos(theta)
        else    # if m == -4 || m == 4
            return 3/16*sqrt(35/pi/2)*exp(im*m*phi)*sin(theta)^4
        end

    elseif l == 5

        if m == 0
            return 1/16*sqrt(11/pi)*(63*cos(theta)^5 - 70*cos(theta)^3 + 15*cos(theta)) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*1/16*sqrt(165/pi/2)*exp(im*m*phi)*sin(theta)*(21*cos(theta)^4 - 14*cos(theta)^2 + 1)
        elseif m == -2 || m == 2
            return 1/8*sqrt(1155/pi/2)*exp(im*m*phi)*sin(theta)^2*(3*cos(theta)^3 - cos(theta))
        elseif m == -3 || m == 3
            return sign(-m)*1/32*sqrt(385/pi)*exp(im*m*phi)*sin(theta)^3*(9*cos(theta)^2 - 1)
        elseif m == -4 || m == 4
            return 3/16*sqrt(385/pi/2)*exp(im*m*phi)*sin(theta)^4*cos(theta)
        else    # if m == -5 || m == 5
            return sign(-m)*3/32*sqrt(77/pi)*exp(im*m*phi)*sin(theta)^5
        end

    elseif l == 6

        if m == 0
            return 1/32*sqrt(13/pi)*(231*cos(theta)^6 - 315*cos(theta)^4 + 105*cos(theta)^2 - 5) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*1/16*sqrt(273/pi/2)*exp(im*m*phi)*sin(theta)*(33*cos(theta)^5 - 30*cos(theta)^3 + 5*cos(theta))
        elseif m == -2 || m == 2
            return 1/64*sqrt(1365/pi)*exp(im*m*phi)*sin(theta)^2*(33*cos(theta)^4 - 18*cos(theta)^2 + 1)
        elseif m == -3 || m == 3
            return sign(-m)*1/32*sqrt(1365/pi)*exp(im*m*phi)*sin(theta)^3*(11*cos(theta)^3 - 3*cos(theta))
        elseif m == -4 || m == 4
            return 3/32*sqrt(91/pi/2)*exp(im*m*phi)*sin(theta)^4*(11*cos(theta)^2 - 1)
        elseif m == -5 || m == 5
            return sign(-m)*3/32*sqrt(1001/pi)*exp(im*m*phi)*sin(theta)^5*cos(theta)
        else    # if m == -6 || m == 6
            return 1/64*sqrt(3003/pi)*exp(im*m*phi)*sin(theta)^6
        end

    elseif l == 7

        if m == 0
            return 1/32*sqrt(15/pi)*(429*cos(theta)^7- 693*cos(theta)^5 + 315*cos(theta)^3 - 35*cos(theta)) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*1/64*sqrt(105/pi/2)*exp(im*m*phi)*sin(theta)*(429*cos(theta)^6 - 495*cos(theta)^4 + 135*cos(theta)^2 - 5)
        elseif m == -2 || m == 2
            return 3/64*sqrt(35/pi)*exp(im*m*phi)*sin(theta)^2*(143*cos(theta)^5 - 110*cos(theta)^3 + 15*cos(theta))
        elseif m == -3 || m == 3
            return sign(-m)*3/64*sqrt(35/pi/2)*exp(im*m*phi)*sin(theta)^3*(143*cos(theta)^4 - 66*cos(theta)^2 + 3)
        elseif m == -4 || m == 4
            return 3/32*sqrt(385/pi/2)*exp(im*m*phi)*sin(theta)^4*(13*cos(theta)^3 - 3*cos(theta))
        elseif m == -5 || m == 5
            return sign(-m)*3/64*sqrt(385/pi/2)*exp(im*m*phi)*sin(theta)^5*(13*cos(theta)^2 - 1)
        elseif m == -6 || m == 6
            return 3/64*sqrt(5005/pi)*exp(im*m*phi)*sin(theta)^6*cos(theta)
        else    # if m == -7 || m == 7
            return sign(-m)*3/64*sqrt(715/pi/2)*exp(im*m*phi)*sin(theta)^7
        end

    elseif l == 8

        if m == 0
            return 1/256*sqrt(17/pi)*(6435*cos(theta)^8 - 12012*cos(theta)^6 + 6930*cos(theta)^4 - 1260*cos(theta)^2 + 35) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*3/64*sqrt(17/pi/2)*exp(im*m*phi)*sin(theta)*(715*cos(theta)^7 - 1001*cos(theta)^5 + 385*cos(theta)^3 - 35*cos(theta))
        elseif m == -2 || m == 2
            return 3/128*sqrt(595/pi)*exp(im*m*phi)*sin(theta)^2*(143*cos(theta)^6 - 143*cos(theta)^4 + 33*cos(theta)^2 - 1)
        elseif m == -3 || m == 3
            return sign(-m)*1/64*sqrt(19635/pi/2)*exp(im*m*phi)*sin(theta)^3*(39*cos(theta)^5 - 26*cos(theta)^3 + 3*cos(theta))
        elseif m == -4 || m == 4
            return 3/128*sqrt(1309/pi/2)*exp(im*m*phi)*sin(theta)^4*(65*cos(theta)^4 - 26*cos(theta)^2 + 1)
        elseif m == -5 || m == 5
            return sign(-m)*3/64*sqrt(17017/pi/2)*exp(im*m*phi)*sin(theta)^5*(5*cos(theta)^3 - cos(theta))
        elseif m == -6 || m == 6
            return 1/128*sqrt(7293/pi)*exp(im*m*phi)*sin(theta)^6*(15*cos(theta)^2 - 1)
        elseif m == -7 || m == 7
            return sign(-m)*3/64*sqrt(12155/pi/2)*exp(im*m*phi)*sin(theta)^7*cos(theta)
        else    # if m == -8 || m == 8
            return 3/256*sqrt(12155/pi/2)*exp(im*m*phi)*sin(theta)^8
        end

    elseif l == 9

        if m == 0
            return 1/256*sqrt(19/pi)*(12155*cos(theta)^9 - 25740*cos(theta)^7 + 18018*cos(theta)^5 - 4620*cos(theta)^3 + 315*cos(theta)) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*3/256*sqrt(95/pi/2)*exp(im*m*phi)*sin(theta)*(2431*cos(theta)^8 - 4004*cos(theta)^6 + 2002*cos(theta)^4 - 308*cos(theta)^2 + 7)
        elseif m == -2 || m == 2
            return 3/128*sqrt(1045/pi)*exp(im*m*phi)*sin(theta)^2*(221*cos(theta)^7 - 273*cos(theta)^5 + 91*cos(theta)^3 - 7*cos(theta))
        elseif m == -3 || m == 3
            return sign(-m)*1/256*sqrt(21945/pi)*exp(im*m*phi)*sin(theta)^3*(221*cos(theta)^6 - 195*cos(theta)^4 + 39*cos(theta)^2 - 1)
        elseif m == -4 || m == 4
            return 3/128*sqrt(95095/pi/2)*exp(im*m*phi)*sin(theta)^4*(17*cos(theta)^5 - 10*cos(theta)^3 + cos(theta))
        elseif m == -5 || m == 5
            return sign(-m)*3/256*sqrt(2717/pi)*exp(im*m*phi)*sin(theta)^5*(85*cos(theta)^4 - 30*cos(theta)^2 + 1)
        elseif m == -6 || m == 6
            return 1/128*sqrt(40755/pi)*exp(im*m*phi)*sin(theta)^6*(17*cos(theta)^3 - 3*cos(theta))
        elseif m == -7 || m == 7
            return sign(-m)*3/512*sqrt(13585/pi)*exp(im*m*phi)*sin(theta)^7*(17*cos(theta)^2 - 1)
        elseif m == -8 || m == 8
            return 3/256*sqrt(230945/pi/2)*exp(im*m*phi)*sin(theta)^8*cos(theta)    
        else    # if m == -9 || m == 9
            return sign(-m)*1/512*sqrt(230945/pi)*exp(im*m*phi)*sin(theta)^9
        end

    elseif l == 10

        if m == 0
            return 1/512*sqrt(21/pi)*(46189*cos(theta)^10 - 109395*cos(theta)^8 + 90090*cos(theta)^6 - 30030*cos(theta)^4 + 3465*cos(theta)^2 - 63) + im*0.0
        elseif m == -1 || m == 1
            return sign(-m)*1/256*sqrt(1155/pi/2)*exp(im*m*phi)*sin(theta)*(4199*cos(theta)^9 - 7956*cos(theta)^7 + 4914*cos(theta)^5 - 1092*cos(theta)^3 + 63*cos(theta))
        elseif m == -2 || m == 2
            return 3/512*sqrt(385/pi/2)*exp(im*m*phi)*sin(theta)^2*(4199*cos(theta)^8 - 6188*cos(theta)^6 + 2730*cos(theta)^4 - 364*cos(theta)^2 + 7)
        elseif m == -3 || m == 3
            return sign(-m)*3/256*sqrt(5005/pi)*exp(im*m*phi)*sin(theta)^3*(323*cos(theta)^7 - 357*cos(theta)^5 + 105*cos(theta)^3 - 7*cos(theta))
        elseif m == -4 || m == 4
            return 3/256*sqrt(5005/pi/2)*exp(im*m*phi)*sin(theta)^4*(323*cos(theta)^6 - 255*cos(theta)^4 + 45*cos(theta)^2 - 1)
        elseif m == -5 || m == 5
            return sign(-m)*3/256*sqrt(1001/pi)*exp(im*m*phi)*sin(theta)^5*(323*cos(theta)^5 - 170*cos(theta)^3 + 15*cos(theta))
        elseif m == -6 || m == 6
            return 3/1024*sqrt(5005/pi)*exp(im*m*phi)*sin(theta)^6*(323*cos(theta)^4 - 102*cos(theta)^2 + 3)
        elseif m == -7 || m == 7
            return sign(-m)*3/512*sqrt(85085/pi)*exp(im*m*phi)*sin(theta)^7*(19*cos(theta)^3 - 3*cos(theta))
        elseif m == -8 || m == 8
            return 1/512*sqrt(255255/pi/2)*exp(im*m*phi)*sin(theta)^8*(19*cos(theta)^2 - 1)
        elseif m == -9 || m == 9
            return sign(-m)*1/512*sqrt(4849845/pi)*exp(im*m*phi)*sin(theta)^9*cos(theta)
        else    # if m == -10 || m == 10
            return 1/1024*sqrt(969969/pi)*exp(im*m*phi)*sin(theta)^10
        end
    else
        error("please check l.")
    end
end


# calculate complex spherical harmonics
function calc_Ylm!(l, m, fact, theta, phi, SH)
    
    tmp = sqrt((2*l+1)/(4*pi))*fact
    Plm = Associated_Legendre(l, abs(m), cos(theta))

    co = cos(m*phi)
    si = sin(m*phi)

    if m >= 0
        SH[1] = tmp*Plm*co
        SH[2] = tmp*Plm*si
    else
        if iszero(abs(m)%2)
            SH[1] = tmp*Plm*co
            SH[2] = tmp*Plm*si
        else
            SH[1] = -tmp*Plm*co
            SH[2] = -tmp*Plm*si
        end
    end
end


function calc_Ylm!(l, m, fact, theta, phi, SH, dSHt, dSHp)
    
    tmp = sqrt((2*l+1)/(4*pi))*fact
    Plm, dPlm = Associated_Legendre2(l, abs(m), cos(theta))

    co = cos(m*phi)
    si = sin(m*phi)

    if m >= 0
        SH[1]   = tmp*Plm*co
        SH[2]   = tmp*Plm*si
        dSHt[1] = tmp*dPlm*co
        dSHt[2] = tmp*dPlm*si
        dSHp[1] = -m*tmp*Plm*si
        dSHp[2] =  m*tmp*Plm*co
    else
        if iszero(abs(m)%2)
            SH[1]   = tmp*Plm*co
            SH[2]   = tmp*Plm*si
            dSHt[1] = tmp*dPlm*co
            dSHt[2] = tmp*dPlm*si
            dSHp[1] = -m*tmp*Plm*si
            dSHp[2] =  m*tmp*Plm*co
        else
            SH[1]   = -tmp*Plm*co
            SH[2]   = -tmp*Plm*si
            dSHt[1] = -tmp*dPlm*co
            dSHt[2] = -tmp*dPlm*si
            dSHp[1] =  m*tmp*Plm*si
            dSHp[2] = -m*tmp*Plm*co
        end
    end
end


function Ylm_complex(l::Integer, m::Integer, x, y, z)
    
    _, theta, phi = xyz_to_spherical(x, y, z)

    if 0 <= l <= 10
        return Ylm_table(l,m,theta,phi)
    else
        error("Invalid l, m, please check Ylm_complex.jl.")
    end
end

# For Set_ProExpn_VNA.jl
function Ylm_complex!(l::Integer, m::Integer, fact, x, y, z, SH)
    
    _, theta, phi = xyz_to_spherical(x, y, z)

    if 0 <= l <= 10
        ylm = Ylm_table(l,m,theta,phi)
        SH[1] = real(ylm)
        SH[2] = imag(ylm)
    else
        calc_Ylm!(l,m,fact,theta,phi,SH)
    end
end


# For Force.jl
function Ylm_complex!(l::Integer, m::Integer, fact, theta, phi, SH, dSHt, dSHp)
    calc_Ylm!(l,m,fact,theta,phi,SH,dSHt,dSHp)
end