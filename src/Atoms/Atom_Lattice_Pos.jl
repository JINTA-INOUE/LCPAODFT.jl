struct Atompos
    Natom::Int32
    Gxyz::Vector{Vector{Float64}}
    Gxyz_frac::Vector{Vector{Float64}}
    unit::String
end


function Base.:*(Gxyz::Vector{Vector{T}}, unit::AbstractString) where {T <: Real}
    
    Natom = length(Gxyz)
    unit = lowercase(strip(unit))
    
    Gxyz_frac = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Gxyz_frac[atom] = [-99999.0, -99999.0, -99999.0]
    end

    if unit == "au"
        unit = "au"
	elseif unit == "ang"
		unit = "ang"
        for atom = 1:Natom
            Gxyz[atom] = Gxyz[atom]*Ang_to_bohr
        end
	elseif unit == "frac"
		unit = "frac"

        for atom = 1:Natom
            Gxyz_frac[atom][1] = Gxyz[atom][1]
            Gxyz_frac[atom][2] = Gxyz[atom][2]
            Gxyz_frac[atom][3] = Gxyz[atom][3]
        end

        for atom = 1:Natom
            for i = 1:3
                if !(zero(T) <= Gxyz[atom][i] <= one(T))
                    error("please check Atom positions")
                end

                if Gxyz[atom][i] > one(T)
                    Gxyz[atom][i] = Gxyz[atom][i] - trunc(Gxyz[atom][i])
                elseif Gxyz[atom][i] < zero(T)
                    Gxyz[atom][i] = Gxyz[atom][i] + trunc(Gxyz[atom][i]) + one(T)
                end

            end
        end
	else
        println("not support ", unit)
		error("please check Atom positions unit")
	end

    return Atompos(Natom, Gxyz, Gxyz_frac, unit)
end


struct Lattice
    Latvecs::Matrix{Float64}
    unit::String
    dim::Int32
end


function Base.:*(Latvecs::Matrix{T}, unit::String) where{T<:Real}
    unit = lowercase(unit)
	if unit == "au"
		unit = "au"
	elseif unit == "ang"
		unit = "ang"
        Latvecs = Latvecs*Ang_to_bohr
	else
        println("not support ", unit)
		error("please check Lattice Vector unit")
	end

    if size(Latvecs) == (1,1)
        dim = 1
    elseif size(Latvecs) == (2,2)
        dim = 2
    elseif size(Latvecs) == (3,3)
        dim = 3
    else
        println("now Latvecs = $Latvecs")
        error("please check lattice Vectors")
    end
        
    return Lattice(Latvecs, unit, dim)
end
