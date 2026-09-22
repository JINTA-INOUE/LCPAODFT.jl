function parse_vector(io::IO, T::Type, n_elements::Integer)
    vec = zeros(T, n_elements)

    counter = 0
    while counter < n_elements
        !eof(io) || error("unexpected end of file")
        line = strip(readline(io))
        splitted = split(line)
        n_splitted = length(splitted)
        vec[(counter + 1):(counter + n_splitted)] = parse.(T, splitted)
        counter += n_splitted
    end

    return vec
end

@inline function parse_vector(parts::AbstractVector{<:AbstractString}, T::Type = Float64)
    return map(x -> parse(T, x), parts)
end
@inline function parse_vector(line::AbstractString, T::Type = Float64)
    return parse_vector(split(line), T)
end


function chk2MLWF_model(mlwf_setup::MLWF_Setup, seedname::String)
    
    MLWF_Outer_Window_Bottom = mlwf_setup.MLWF_Outer_Window_Bottom
    MLWF_Outer_Window_Top = mlwf_setup.MLWF_Outer_Window_Top
    MLWF_Inner_Window_Bottom = mlwf_setup.MLWF_Inner_Window_Bottom
    MLWF_Inner_Window_Top = mlwf_setup.MLWF_Inner_Window_Top
    Dis_Energy = [MLWF_Outer_Window_Bottom, MLWF_Inner_Window_Bottom, MLWF_Inner_Window_Top, MLWF_Outer_Window_Top]
    kmesh = mlwf_setup.MLWF_kmesh
    material = mlwf_setup.material
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    SO_switch = material.SO_switch
    ChemP = material.ChemP
    scf_inputfile = material.scf_inputfile
    mlwf_inputfile = [""]

    data = open("$(seedname)_hr.dat", "r")
    header = strip(readline(data))
    Nwann = parse(Int32, strip(readline(data)))
    NCell = parse(Int32, strip(readline(data)))

    Rdegens = parse_vector(data, Int, NCell)
    cell_list_ijk = Vector{Vector{Int32}}(undef, NCell)
    for cell = 1:NCell
        cell_list_ijk[cell] = zeros(Int32, 3)
    end

    HmnR = zeros(ComplexF64, Nwann, Nwann, NCell, spinsize)
    for cell = 1:NCell, n = 1:Nwann, m = 1:Nwann
        line = split(strip(readline(data)))
        cell_list_ijk[cell][1] = parse(Int, line[1])
        cell_list_ijk[cell][2] = parse(Int, line[2])
        cell_list_ijk[cell][3] = parse(Int, line[3])
        m == parse(Int, line[4]) || error(line)
        n == parse(Int, line[5]) || error(line)
        HmnR[m,n,cell,1] = complex(parse(Float64, line[6]), parse(Float64, line[7]))
    end
    close(data)

    println("\tWrite $(seedname).MLWF.jld2")
    jldopen("$(seedname).MLWF.jld2", "w") do file
        file["Dates"] = now()
        file["spinsize"] = spinsize
        file["Nwann"] = Nwann
        file["Latvecs"] = Latvecs
        file["Recvecs"] = Recvecs
        file["SpinPol"] = SpinPol
        file["SO_switch"] = SO_switch
        file["kmesh"] = kmesh
        file["Dis_Energy"] = Dis_Energy
        file["NCell"] = NCell
        file["cell_list_ijk"] = cell_list_ijk
        file["Rdegens"] = Rdegens
        file["HmnR"] = HmnR
        file["ChemP"] = ChemP
        file["scf_inputfile"] = scf_inputfile
        file["mlwf_inputfile"] = mlwf_inputfile
    end
end