function Calc_Grid_Origin(Natom, Gxyz_AU, gLatvecs, Ngrid, unit)

    Grid_Origin = zeros(Float64, 3)
    if lowercase(unit) ≠ "frac"
        xc = yc = zc = 0.0
        for atom = 1:Natom
            xc += Gxyz_AU[atom][1]
            yc += Gxyz_AU[atom][2]
            zc += Gxyz_AU[atom][3]
        end
        xc = xc/Natom
        yc = yc/Natom
        zc = zc/Natom
        sn1 = 0.5*mod(Ngrid[1]+1, 2)
        sn2 = 0.5*mod(Ngrid[2]+1, 2)
        sn3 = 0.5*mod(Ngrid[3]+1, 2)

        xm = (div(Ngrid[1],2)-sn1)*gLatvecs[1,1] + (div(Ngrid[2],2)-sn2)*gLatvecs[2,1] + (div(Ngrid[3],2)-sn3)*gLatvecs[3,1]
        ym = (div(Ngrid[1],2)-sn1)*gLatvecs[1,2] + (div(Ngrid[2],2)-sn2)*gLatvecs[2,2] + (div(Ngrid[3],2)-sn3)*gLatvecs[3,2]
        zm = (div(Ngrid[1],2)-sn1)*gLatvecs[1,3] + (div(Ngrid[2],2)-sn2)*gLatvecs[2,3] + (div(Ngrid[3],2)-sn3)*gLatvecs[3,3]
        @. Grid_Origin = [xc-xm, yc-ym, zc-zm]
    end

    return Grid_Origin
end