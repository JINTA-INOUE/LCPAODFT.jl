function Calc_Gxyz_frac(Natom, Gxyz_AU, Recvecs)

    Gxyz_frac = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Gxyz_frac[atom] = zeros(Float64, 3)
    end

    for atom = 1:Natom    
        Gxyz_frac[atom][1] = dot(Gxyz_AU[atom], Recvecs[1,:])*0.5/pi
        Gxyz_frac[atom][2] = dot(Gxyz_AU[atom], Recvecs[2,:])*0.5/pi
        Gxyz_frac[atom][3] = dot(Gxyz_AU[atom], Recvecs[3,:])*0.5/pi
        for i = 1:3
            tmp = floor(Int64, Gxyz_frac[atom][i])
            if Gxyz_frac[atom][i] > 1.0
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]-tmp)
            elseif Gxyz_frac[atom][i] < -1e-13
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]+abs(tmp)+1)
            end
        end
    end


    return Gxyz_frac
end