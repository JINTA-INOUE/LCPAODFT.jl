function xyz_to_spherical(x, y, z)

    Min_r = 1e-15

    dum = x^2 + y^2
    r = sqrt(dum + z^2)
    r1 = sqrt(dum)

    if Min_r <= r

        dum1 = ifelse(r < abs(z), sign(z), z/r)
        theta = acos(dum1)

        if Min_r <= r1
            if x >= 0.0
                dum1 = ifelse(r1 < abs(y), sign(y), y/r1)
                phi = asin(dum1)
            else
                dum1 = ifelse(r1 < abs(y), sign(y), y/r1)
                phi = pi - asin(dum1)
            end
        else
            phi = 0.0
        end
    else
        theta = 0.5*pi
        phi = 0.0
    end

    return r, theta, phi
end