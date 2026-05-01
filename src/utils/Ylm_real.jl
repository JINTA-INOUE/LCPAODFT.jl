function Ylm_real(l::Integer, m::Integer, theta, phi)

    if l < abs(m)
        error("please check l, m.")
    end

    siQ = sin(theta)
    coQ = cos(theta)
    siP = sin(phi)
    coP = cos(phi)

    if l == 0
        return 0.5*sqrt(1/pi)
    elseif l == 1
        if m == -1         # py
            return sqrt(3/4/pi)*siQ*siP
        elseif m == 0      # pz
            return sqrt(3/4/pi)*coQ
        else               # px
            return sqrt(3/4/pi)*siQ*coP
        end
    elseif l == 2
        if m == -2         # dxy
            return 0.5*sqrt(15/pi) * siQ^2 * siP * coP
        elseif m == -1     # dyz
            return 0.5*sqrt(15/pi) * siQ * coQ * siP
        elseif m == 0      # dz^2
            return 0.25*sqrt(5/pi) * (3*coQ^2 - 1)
        elseif m == 1      # dxz
            return 0.5*sqrt(15/pi) * siQ * coQ * coP
        else               # dx^2-y^2
            return 0.25*sqrt(15/pi) * siQ^2 * (1 - 2*siP^2)
        end   
    elseif l == 3
        if m == -3          # y(3x^2-y^2)
            return 0.25*sqrt(35/2/pi) * siQ^3 * (3*siP - 4*siP^3)
        elseif m == -2      # xyz
            return 0.5*sqrt(105/pi) * siQ^2 * coQ * siP * coP
        elseif m == -1      # yz^2
            return 0.25*sqrt(21/2/pi)*siQ*( 5*coQ^2 - 1 )*siP
        elseif m == 0       # z^3
            return 0.25*sqrt(7/pi)*( 5*coQ^3 - 3*coQ )
        elseif m == 1       # xz^2
            return 0.25*sqrt(21/2/pi)*siQ*( 5*coQ^2 - 1 )*coP
        elseif m == 2       # z(x^2-y^2)
            return 0.25*sqrt(105/pi) * siQ^2 * coQ * (coP^2 - siP^2)
        else                # x(x^2-3y^2)
            return 0.25*sqrt(35/2/pi) * siQ^3 * (4*coP^3 - 3*coP)   
        end
    else
        error("please check Ylm_real.jl")
    end
end


function dYlmdtheta_real(l::Int64, m::Int64, theta, phi)

    if l < abs(m)
        error("please check l, m.")
    end

    siQ = sin(theta)
    coQ = cos(theta)
    siP = sin(phi)
    coP = cos(phi)

    if l == 0
        return 0.0
    elseif l == 1
        if m == -1         # py
            return sqrt(3/4/pi)*coQ*siP
        elseif m == 0      # pz
            return -sqrt(3/4/pi)*siQ
        else               # px
            return sqrt(3/4/pi)*coQ*coP
        end
    elseif l == 2
        if m == -2         # dxy
            # return 0.5*sqrt(15/pi) * siQ * coQ * siP * coP
            return sqrt(15/pi) * siQ * coQ * siP * coP
        elseif m == -1     # dyz
            return 0.5*sqrt(15/pi) * (1 - 2*siQ^2) * siP
        elseif m == 0      # dz^2
            return -1.5*sqrt(5/pi) * coQ * siQ
        elseif m == 1      # dxz
            return 0.5*sqrt(15/pi) * (1 - 2*siQ^2) * coP
        else               # dx^2-y^2
            # return 0.25*sqrt(15/pi) * siQ*coQ*(1 - 2*siP^2)
            return 0.5*sqrt(15/pi) * siQ*coQ*(1 - 2*siP^2)
        end
    elseif l == 3
        if m == -3          # y(3x^2-y^2)
            return 0.75*sqrt(35/2/pi)*siQ^2*coQ*siP*(3 - 4*siP^2)
        elseif m == -2      # xyz
            return 0.5*sqrt(105/pi)*coP*siP*siQ*(2*coQ^2 - siQ^2)
        elseif m == -1      # yz^2
            return 0.25*sqrt(21/2/pi)*siP*coQ*(15*coQ^2 - 11)
        elseif m == 0       # z^3
            return 0.25*sqrt(7/pi)*siQ*(-15*coQ^2 + 3)
        elseif m == 1       # xz^2
            return 0.25*sqrt(21/2/pi)*coP*coQ*(15*coQ^2 - 11)
        elseif m == 2       # z(x^2-y^2)
            return 0.25*sqrt(105/pi)*(coP^2-siP^2)*siQ*(2*coQ^2-siQ^2)
        else                # x(x^2-3y^2)
            return 0.75*sqrt(35/2/pi)*coP*coQ*siQ^2*(-3 + 4*coP^2) 
        end
    else
        error("please check Ylm_real.jl")
    end
end


function dYlmdphi_real(l::Int64, m::Int64, theta, phi)

    if l < abs(m)
        error("please check l, m.")
    end

    siQ = sin(theta)
    coQ = cos(theta)
    siP = sin(phi)
    coP = cos(phi)

    if l == 0
        return 0.0
    elseif l == 1
        if m == -1         # py
            return -sqrt(3/4/pi)*siP
        elseif m == 0      # pz
            return 0.0
        else               # px
            return sqrt(3/4/pi)*coP
        end
    elseif l == 2
        if m == -2         # dxy
            return 0.5*sqrt(15/pi)*siQ*(1 - 2*siP^2)
        elseif m == -1     # dyz
            return 0.5*sqrt(15/pi)*coQ*coP
        elseif m == 0      # dz^2
            return 0.0
        elseif m == 1      # dxz
            return -0.5*sqrt(15/pi)*coQ*siP
        else               # dx^2-y^2
            return -sqrt(15/pi)*siQ*siP*coP
        end
    elseif l == 3
        if m == -3          # y(3x^2-y^2)
            return 0.75*sqrt(35/2/pi)*coP*siQ^2*(1 - 4*siP^2)
        elseif m == -2      # xyz
            return 0.5*sqrt(105/pi)*coQ*siQ*(coP^2 - siP^2)
        elseif m == -1      # yz^2
            return 0.25*sqrt(21/2/pi)*coP*( 5*coQ^2 - 1)
        elseif m == 0      # z^3
            return 0.0
        elseif m == 1       # xz^2
            return 0.25*sqrt(21/2/pi)*siP*(-5*coQ^2 + 1)
        elseif m == 2       # z(x^2-y^2)
            return -sqrt(105/pi) * coP*coQ*siP*siQ
        else                # x(x^2-3y^2)
            return 0.75*sqrt(35/2/pi) * siP*siQ^2*(1 - 4*coP^2)      
        end
    else
        error("please check Ylm_real.jl")
    end
end