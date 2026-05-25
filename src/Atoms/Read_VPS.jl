struct Pspot
    Atom_symbol::String
    Atom_extra::String
    xc_type::String
    Spe_Core_Charge::Float64
    VPS_j_dependency::Int32
    SO_switch::Bool
    Spe_Num_Mesh_VPS::Int64
    Spe_Num_RVPS::Int32
    Spe_VPS_List::Vector{Int64}
    Spe_VNLE::Array{Float64,2}
    Spe_VNL::Vector{Vector{Vector{Float64}}}
    Spe_VPS_XV::Vector{Float64}
    Spe_VPS_RV::Vector{Float64}
    Spe_Vcore::Vector{Float64}
    is_pcc::Bool
    Spe_Atomic_PCC::Vector{Float64}
    psfile::String
end


function Print_Pspot(pspot::Pspot)

    Atom_symbol = pspot.Atom_symbol
    Atom_extra = pspot.Atom_extra
    xc_type = pspot.xc_type
    VPS_j_dependency = pspot.VPS_j_dependency
    SO_switch = pspot.SO_switch
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    is_pcc = pspot.is_pcc
    psfile = pspot.psfile

    println("<Pspot>")
    println("\tAtom_symbol : $(Atom_symbol)")
    if Atom_extra ≠ ""
        println("\tAtom_extra : $(Atom_extra)")
    end
    println("\tfilepath : $(psfile)")
    println("\txc_type : $(xc_type)")
    println("\tVPS_j_dependency : $(VPS_j_dependency)")
    println("\tSO_switch : $(SO_switch)")
    println("\tSpe_Num_Mesh_VPS : $(Spe_Num_Mesh_VPS)")
    println("\tis_pcc : $(is_pcc)")
    println("")
end


function _Read_VPS(psfile)
    
    mark_AtomSpecies = "AtomSpecies"
    mark_total_electron = "total.electron"
    mark_valence_electron = "valence.electron"
    mark_grid_num_output = "grid.num.output"
    mark_charge_pcc_calc = "charge.pcc.calc"
    mark_begin_project_energies = "<project.energies"
    mark_end_project_energies = "project.energies>"
    mark_begin_Pseudo_Potentials = "<Pseudo.Potentials"
    mark_end_Pseudo_Potentials = "Pseudo.Potentials>"
    mark_begin_density_PCC = "<density.PCC"
    mark_end_density_PCC = "density.PCC>"

    po = zeros(Int64, 5)
    AtomSpecies = 100.0
    tot_ele = 100.0
    val_ele = 100.0
    Spe_Num_Mesh_VPS = -100
    is_pcc = false

    open(psfile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_AtomSpecies, line)
                AtomSpecies = parse(Float64, split(line)[2])
                po[1] = 1
                continue
            end

            if occursin(mark_total_electron, line)
                tot_ele = parse(Float64, split(line)[2])
                po[2] = 1
                continue
            end

            if occursin(mark_valence_electron, line)
                val_ele = parse(Float64, split(line)[2])
                po[3] = 1
                continue
            end

            if occursin(mark_grid_num_output, line)
                Spe_Num_Mesh_VPS = parse(Int64, split(line)[2])
                po[4] = 1
                continue
            end

            if occursin(mark_charge_pcc_calc, line)
                temp = lowercase(split(line)[2])
                if temp == "on"
                    is_pcc = true
                else
                    is_pcc = false
                end
                po[5] = 1
                continue
            end
        end
    end

    Spe_Core_Charge = val_ele - tot_ele + AtomSpecies
    if AtomSpecies == 100.0 || tot_ele == 100.0 || val_ele == 100.0
        @show AtomSpecies
        @show tot_ele
        @show val_ele
        @show po
        error("please check AtomSpecies total.electron valence.electron")
    end

    if Spe_Num_Mesh_VPS == -100
        error("check idata loop grid.num.output")
    end

    if prod(po) == 0
        error("please check")
    end



    # get Spe_Num_RVPS and Spe_VNLE
    po = 0
    Spe_Num_RVPS = -100
    Spe_VNLE_temp = zeros(Float64, 2, 50)
    Spe_VPS_List_temp = zeros(Int64, 50)

    open(psfile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_begin_project_energies, line)
                line = srline()
                Spe_Num_RVPS = parse(Int32, split(line)[1])
                for j = 1:Spe_Num_RVPS
                    line = srline()
                    Spe_VPS_List_temp[j] = parse(Int64, split(line)[1])
                    Spe_VNLE_temp[1,j] = parse(Float64, split(line)[2])
                    Spe_VNLE_temp[2,j] = parse(Float64, split(line)[3])
                end

                line = srline()
                if !occursin(mark_end_project_energies, line)
                    error("please check project_energies")
                else
                    po = 1
                end
                break
            end
        end
    end

    if Spe_Num_RVPS == -100 || po == 0
        error("please check project.energies")
    end

    @views Spe_VPS_List = Spe_VPS_List_temp[1:Spe_Num_RVPS]

    Spe_VNLE_wSOC = zeros(Float64, 2, Spe_Num_RVPS)
    Spe_VNLE_woSOC = zeros(Float64, 2, Spe_Num_RVPS)

    @. @views Spe_VNLE_wSOC[1,:] = Spe_VNLE_temp[1,1:Spe_Num_RVPS]
    @. @views Spe_VNLE_wSOC[2,:] = Spe_VNLE_temp[2,1:Spe_Num_RVPS]

    for l = 1:Spe_Num_RVPS
        LVPS = Spe_VPS_List[l]
        tmp = ((LVPS+1)*Spe_VNLE_wSOC[1,l]+LVPS*Spe_VNLE_wSOC[2,l])/(2*LVPS+1)
        Spe_VNLE_woSOC[1,l] = tmp
        Spe_VNLE_woSOC[2,l] = tmp
    end
    



    Spe_VNL_woSOC = Vector{Vector{Vector{Float64}}}(undef, 1)
    Spe_VNL_wSOC = Vector{Vector{Vector{Float64}}}(undef, 2)
    Spe_VNL_woSOC[1] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
    Spe_VNL_wSOC[1] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
    Spe_VNL_wSOC[2] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
    for l = 1:Spe_Num_RVPS
        Spe_VNL_woSOC[1][l] = zeros(Float64, Spe_Num_Mesh_VPS)
        Spe_VNL_wSOC[1][l] = zeros(Float64, Spe_Num_Mesh_VPS)
        Spe_VNL_wSOC[2][l] = zeros(Float64, Spe_Num_Mesh_VPS)
    end


    po = 0
    Spe_VPS_XV = zeros(Float64, Spe_Num_Mesh_VPS)
    Spe_VPS_RV = zeros(Float64, Spe_Num_Mesh_VPS)
    Spe_Vcore = zeros(Float64, Spe_Num_Mesh_VPS)

    open(psfile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_begin_Pseudo_Potentials, line)
                for i = 1:Spe_Num_Mesh_VPS
                    line = srline()
                    Spe_VPS_XV[i] = parse(Float64, split(line)[1])
                    Spe_VPS_RV[i] = parse(Float64, split(line)[2])
                    Spe_Vcore[i] = parse(Float64, split(line)[3])

                    for l = 1:Spe_Num_RVPS
                        Spe_VNL_wSOC[1][l][i] = parse(Float64, split(line)[2*l+2])
                        Spe_VNL_wSOC[2][l][i] = parse(Float64, split(line)[2*l+3])
                    end
                end

                line = srline()
                if !occursin(mark_end_Pseudo_Potentials, line)
                    error("please check Pseudo_Potentials")
                else
                    po = 1
                end
                break
            end
        end
    end

    if Spe_VPS_RV[1] < 1e-10 || po == 0
        error("please check Pseudo_Potentials")
    end


    for l = 1:Spe_Num_RVPS
        LVPS = Spe_VPS_List[l]
        for i = 1:Spe_Num_Mesh_VPS
            Spe_VNL_woSOC[1][l][i] = ((LVPS+1)*Spe_VNL_wSOC[1][l][i] + LVPS*Spe_VNL_wSOC[2][l][i])/(2*LVPS+1)
        end
    end



    po = 0
    Spe_Atomic_PCC = zeros(Float64, Spe_Num_Mesh_VPS+4)
    
    open(psfile) do io
        srline() = strip(readline(io))     
        while !eof(io)
            line = srline()
            if occursin(mark_begin_density_PCC, line)
                for i = 1:Spe_Num_Mesh_VPS
                    line = srline()
                    Spe_Atomic_PCC[i+1] = parse(Float64, split(line)[3])
                end

                line = srline()
                if !occursin(mark_end_density_PCC, line)
                    error("please check density_PCC")
                else
                    po = 1
                end
                break
            end
        end
    end
    Spe_Atomic_PCC[1] = 2*Spe_Atomic_PCC[2] - Spe_Atomic_PCC[3]
    Spe_Atomic_PCC[Spe_Num_Mesh_VPS+2] = 2*Spe_Atomic_PCC[Spe_Num_Mesh_VPS+1] - Spe_Atomic_PCC[Spe_Num_Mesh_VPS]

    if po == 0 && is_pcc
        error("please check density_PCC")
    end    

    
    return (Spe_Core_Charge, Spe_Num_Mesh_VPS, is_pcc, Spe_Num_RVPS, Spe_VPS_List, 
            Spe_VNLE_woSOC, Spe_VNLE_wSOC, Spe_VNL_woSOC, Spe_VNL_wSOC, Spe_VPS_XV, Spe_VPS_RV, Spe_Vcore, Spe_Atomic_PCC)
end


"""
    pspot = read_VPS(...)

Create an instance of `Pspot`.

Mandatory arguments:

- `Atom_symbol`: Atom symbol (`H`, `He`, `Li` ...)
- `Atom_extra`: Atom symbol (`H`, `S`)
- `xc_type`: which use xc type (`LDA`, `LSDA`, `GGA_PBE`)  
             default `LDA`  
             if use `LDA` or `LSDA` when CA type calculation

# Examples
```
julia> Read_VPS("C", "", "LDA", false)
julia> Read_VPS("Fe", "S", "LDA", false)
```
"""
@timeit timer "Read_VPS" function Read_VPS( 
    Atom_symbol::String, 
    Atom_extra::String, 
    xc_type::String, 
    SO_switch::Bool;
    verbosity = 1)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    filename = Atom_symbol*"_"
    if xc_type ∈ ["LDA", "LSDA"]
        filename = filename*"CA19"*Atom_extra*".vps"
    elseif xc_type ∈ ["GGA_PBE"]
        filename = filename*"PBE19"*Atom_extra*".vps"
    else
        println("xc_type is $xc_type")
        error("please check")
    end


    VPS_Name = split(filename, ".vps")[1]
    if VPS_Name ∉ VPS_Files
        error("$VPS_Name cannot find. No such file or directory")
    end


    (Spe_Core_Charge, Spe_Num_Mesh_VPS, is_pcc, Spe_Num_RVPS, Spe_VPS_List, 
    Spe_VNLE_woSOC, Spe_VNLE_wSOC, Spe_VNL_woSOC, Spe_VNL_wSOC, 
    Spe_VPS_XV, Spe_VPS_RV, Spe_Vcore, Spe_Atomic_PCC) = _Read_VPS(VPS_File_path*VPS_Name*".vps")


    if SO_switch
        Spe_VNLE = Spe_VNLE_wSOC
        Spe_VNL = Spe_VNL_wSOC
        VPS_j_dependency = 1
    else
        Spe_VNLE = Spe_VNLE_woSOC
        Spe_VNL = Spe_VNL_woSOC
        VPS_j_dependency = 0
    end


    pspot = Pspot( Atom_symbol, Atom_extra,
                   xc_type, 
                   Spe_Core_Charge,
                   VPS_j_dependency,
                   SO_switch,
                   Spe_Num_Mesh_VPS,
                   Spe_Num_RVPS, Spe_VPS_List,
                   Spe_VNLE, Spe_VNL, 
                   Spe_VPS_XV, Spe_VPS_RV, Spe_Vcore,
                   is_pcc, Spe_Atomic_PCC,
                   VPS_File_path*filename )

    if myrank == 0 && verbosity>=1
        Print_Pspot(pspot)
    end
            
                
    return pspot
end