struct PAO
    Atom_symbol::String
    Atom_orbital::String
    Atom_extra::String
    Spe_MaxL_Basis::Int64
    Spe_Num_Basis::Vector{Int64}
    Spe_Total_NumOrbs::Int32
    Spe_Num_Mesh_PAO::Int64
    Spe_Atom_Cut1::Float64
    Spe_PAO_RV::Vector{Float64}
    Spe_Atomic_Den::Vector{Float64}
    Spe_PAO_Lmax::Int64
    Spe_PAO_Mul::Int64
    Spe_PAO_RWF::Vector{Vector{Vector{Float64}}}
    Spe_RF_Bessel::Vector{Vector{Vector{Float64}}}
    paofile::String
end


function Print_PAO(pao::PAO)

    Atom_symbol = pao.Atom_symbol
    paofile = pao.paofile
    Spe_PAO_Lmax = pao.Spe_PAO_Lmax
    Spe_PAO_Mul = pao.Spe_PAO_Mul
    Atom_orbital = pao.Atom_orbital
    Atom_extra = pao.Atom_extra
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Num_Mesh_PAO = pao.Spe_Num_Mesh_PAO

    println("<PAO>")
    println("\tAtom_symbol : $(Atom_symbol)")
    println("\tfilepath : $(paofile)")
    println("\tSpe_PAO_Lmax : $(Spe_PAO_Lmax)")
    println("\tSpe_PAO_Mul : $(Spe_PAO_Mul)")
    println("\tAtom_orbital : $(Atom_orbital)")
    if Atom_extra ≠ ""
        println("\tAtom_extra : $(Atom_extra)")
    end
    println("\tSpe_Atom_Cut1 : $(Spe_Atom_Cut1)")
    println("\tSpe_Num_Mesh_PAO : $(Spe_Num_Mesh_PAO)")
end


function get_ialpha_index(Spe_orbitals::String)

    if length(Spe_orbitals)%2 == 1
        println("Spe_orbitals = ", Spe_orbitals)
        error("please check input data")
    end


    Spe_MaxL_Basis = div(length(Spe_orbitals),2)-1
    Spe_Num_Basis = zeros(Int64, Spe_MaxL_Basis+1)

    for L = 0:Spe_MaxL_Basis
        l = L+1
        if L == 0
            Spe_Num_Basis[l] = parse(Int64, split(Spe_orbitals[1:2], "s")[2])
        elseif L == 1
            Spe_Num_Basis[l] = parse(Int64, split(Spe_orbitals[3:4], "p")[2])
        elseif L == 2
            Spe_Num_Basis[l] = parse(Int64, split(Spe_orbitals[5:6], "d")[2])
        elseif L == 3
            Spe_Num_Basis[l] = parse(Int64, split(Spe_orbitals[7:8], "f")[2])
        else
            error("please check basis")
        end
    end

    return Spe_MaxL_Basis, Spe_Num_Basis
end



function Read_Spe_PAO_RWF!(filename, L, Spe_Num_Mesh_PAO, Spe_PAO_Mul, Spe_PAO_RWF)
    mark_pseudo_atomic_orbitals = "pseudo.atomic.orbitals.L="
    po0 = 0
    mark_begin_pseudo_atomic_orbitals ="<"*mark_pseudo_atomic_orbitals*string(L)
    mark_end_pseudo_atomic_orbitals = mark_pseudo_atomic_orbitals*string(L)*">"
    open(filename) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()

            if occursin(mark_begin_pseudo_atomic_orbitals, line)
                for i = 1:Spe_Num_Mesh_PAO
                    line = srline()
                    for iMul = 1:Spe_PAO_Mul
                        Spe_PAO_RWF[iMul][i] = parse(Float64, split(line)[iMul+2])
                    end
                end
                
                line = srline()
                if !occursin(mark_end_pseudo_atomic_orbitals, line)
                    error("please check pseudo_atomic_orbitals")
                else
                    po0 = 1
                end
            end
        end
    end

    if po0 == 0
        error("please check pseudo_atomic_orbitals")
    end
end


function _Read_PAO(Atom_cutoff, paofile::String)

    mark_grid_num_output = "grid.num.output"
    mark_radial_cutoff_pao = "radial.cutoff.pao"
    mark_begin_valence_charge_density = "<valence.charge.density"
    mark_end_valence_charge_density = "valence.charge.density>"
    mark_PAO_Lmax = "PAO.Lmax"
    mark_PAO_Mul = "PAO.Mul"


    po = zeros(Int64, 4)
    Spe_Num_Mesh_PAO = 0
    Spe_Atom_Cut1 = -1.0
    Spe_PAO_Lmax = 0
    Spe_PAO_Mul = 0

    open(paofile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_grid_num_output, line)
                Spe_Num_Mesh_PAO = parse(Int64, split(line)[2])
                po[1] = 1
                continue
            end

            if occursin(mark_radial_cutoff_pao, line)
                Spe_Atom_Cut1 = parse(Float64, split(line)[2])
                po[2] = 1
                continue
            end

            if occursin(mark_PAO_Lmax, line)
                Spe_PAO_Lmax = parse(Int64, split(line)[2])
                po[3] = 1
                continue
            end

            if occursin(mark_PAO_Mul, line)
                Spe_PAO_Mul = parse(Int64, split(line)[2])
                po[4] = 1
                continue
            end
        end
    end

    if abs(Spe_Atom_Cut1 - Atom_cutoff) > 1e-10
        @show Spe_Atom_Cut1
        @show Atom_cutoff
        error("Please check Spe_Atom_Cut1")
    end

    if Spe_Num_Mesh_PAO == 0 || Spe_Atom_Cut1 < 0.0
        error("Please check Spe_Num_Mesh_PAO")
    end

    if iszero(Spe_PAO_Lmax) || iszero(Spe_PAO_Mul)
        error("Please check spe_PAO_Lmax and spe_PAO_Mul")
    end

    if prod(po) == 0
        error("please check po")
    end



    po = 0
    Spe_PAO_RV = zeros(Float64, Spe_Num_Mesh_PAO)
    Spe_Atomic_Den = zeros(Float64, Spe_Num_Mesh_PAO+4)

    open(paofile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_begin_valence_charge_density, line)
                for i = 1:Spe_Num_Mesh_PAO
                    line = srline()
                    Spe_PAO_RV[i] = parse(Float64, split(line)[2])
                    Spe_Atomic_Den[i+1] = parse(Float64, split(line)[3])
                end
                
                line = srline()
                if !occursin(mark_end_valence_charge_density, line)
                    error("please check valence_charge_density")
                else
                    po = 1
                end
            end
        end
    end

    if Spe_Atomic_Den[2] < 1e-10 || po == 0
        error("please check Spe_Atomic_Den")
    end


    Spe_PAO_RWF = Vector{Vector{Vector{Float64}}}(undef, Spe_PAO_Lmax+1)
    for L = 0:Spe_PAO_Lmax
        Spe_PAO_RWF[L+1] = Vector{Vector{Float64}}(undef, Spe_PAO_Mul)
        for mul = 1:Spe_PAO_Mul
            Spe_PAO_RWF[L+1][mul] = zeros(Float64, Spe_Num_Mesh_PAO)
        end
    end


    for L = 0:Spe_PAO_Lmax
        Read_Spe_PAO_RWF!(paofile, L, Spe_Num_Mesh_PAO, Spe_PAO_Mul, Spe_PAO_RWF[L+1])
    end

    
    return Spe_Num_Mesh_PAO, Spe_Atom_Cut1, Spe_PAO_RV, Spe_Atomic_Den, Spe_PAO_Lmax, Spe_PAO_Mul, Spe_PAO_RWF
end


"""
    pao = read_PAO(...)

Create an instance of `PAO`.

Mandatory arguments:

- `Spe_Core_Charge` : Core charge
- `Atom_symbol`: Atom symbol (`H`, `He`, `Li` ...)
- `Atom_cutoff`: cut off radius ( Bohr unit )
- `Atom_orb`: Atom orbitals (`s2p1`, `s2p2d1` ...)
- `Atom_extra` : Atom Extra symbol for `Fe`, `Co`, `Ni`, `Cu`, `Zn`

# Examples
```
julia> Read_PAO(4.0, "C", 5.0, "s2p2d1", "")
```
"""
@timeit timer "Read_PAO" function Read_PAO( 
    Spe_Core_Charge, 
    Atom_symbol::String, 
    Atom_cutoff::Float64,
    Atom_orb::String, 
    Atom_extra::String;
    verbosity::Int=1 )

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    

    PAO_Name = Atom_symbol*string(Atom_cutoff)*Atom_extra
    if PAO_Name ∉ PAO_Files
        println("Cutoff radius $Atom_cutoff")
        error("$PAO_Name cannot find. No such file or directory")
    end

    filename = PAO_Name*".pao" 


    Spe_Num_Mesh_PAO, Spe_Atom_Cut1, Spe_PAO_RV, Spe_Atomic_Den, Spe_PAO_Lmax, Spe_PAO_Mul, Spe_PAO_RWF = _Read_PAO(Atom_cutoff, PAO_File_path*filename)


    Spe_MaxL_Basis, Spe_Num_Basis = get_ialpha_index( Atom_orb )
    Spe_Total_NumOrbs = 0
    for l = 0:Spe_MaxL_Basis
        Spe_Total_NumOrbs += Spe_Num_Basis[l+1]*(2*l+1)
    end
    


    # re-normalization of atomic charge density
    Sum = 0.0
    dx = log(Spe_PAO_RV[2]/Spe_PAO_RV[1])
    for i = 1:Spe_Num_Mesh_PAO
        Sum += Spe_Atomic_Den[i+1]*Spe_PAO_RV[i]^3
    end
    Sum *= 4*pi*dx

    @. Spe_Atomic_Den = Spe_Atomic_Den * Spe_Core_Charge/Sum

    Spe_Atomic_Den[1] = 2*Spe_Atomic_Den[2] - Spe_Atomic_Den[3]
    Spe_Atomic_Den[Spe_Num_Mesh_PAO+2] = 2*Spe_Atomic_Den[Spe_Num_Mesh_PAO+1] - Spe_Atomic_Den[Spe_Num_Mesh_PAO]
    

    
    Spe_RF_Bessel = Vector{Vector{Vector{Float64}}}(undef, Spe_MaxL_Basis+1)
    for l = 0:Spe_MaxL_Basis
        Spe_RF_Bessel[l+1] = Vector{Vector{Float64}}(undef, Spe_Num_Basis[l+1])
        for p = 1:Spe_Num_Basis[l+1]
            Spe_RF_Bessel[l+1][p] = zeros(Float64, NkGrid+1)
        end
    end
    FT_PAO!( Spe_Atom_Cut1, Spe_MaxL_Basis, Spe_Num_Basis, Spe_PAO_RV, Spe_PAO_RWF, Spe_RF_Bessel)


    pao = PAO( Atom_symbol, Atom_orb, Atom_extra,
               Spe_MaxL_Basis, Spe_Num_Basis, Spe_Total_NumOrbs,
               Spe_Num_Mesh_PAO, Spe_Atom_Cut1,
               Spe_PAO_RV, Spe_Atomic_Den,
               Spe_PAO_Lmax, Spe_PAO_Mul,
               Spe_PAO_RWF, Spe_RF_Bessel,
               PAO_File_path*filename)

    if myrank == 0 && verbosity>=1
        Print_PAO(pao)
    end

    return pao
end