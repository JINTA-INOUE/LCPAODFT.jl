struct Symmetry
    spacegroup_number::Int64
    hall_number::Int64
    international_symbol::String
    hall_symbol::String
    choice::String
    transformation_matrix::Matrix{Float64}
    origin_shift::Vector{Float64}
    n_atoms::Int64
    wyckoffs::Vector{Char}
    site_symmetry_symbols::Vector{String}
    equivalent_atoms::Vector{Int32}
    crystallographic_orbits::Vector{Int32}
    primitive_lattice::Matrix{Float64}
    mapping_to_primitive::Vector{Int32}
    n_std_atoms::Int64
    std_lattice::Matrix{Float64}
    std_types::Vector{Int32}
    std_positions::Vector{SVector{3, Float64}}
    std_rotation_matrix::Matrix{Float64}
    std_mapping_to_primitive::Vector{Int32}
    pointgroup_symbol::String
    # Nrots::Int32
    # Rots_direct::Vector{Matrix{Int32}}
    # Rots_cartesian::Vector{Matrix{Float64}}
    # translations::Vector{Vector{Float64}}
    # Rots_is_proper::Vector{Bool}
    # Rots_euler::Vector{Vector{Float64}}
    # Rots_Slm::Vector{Vector{Matrix{ComplexF64}}}
    # Rots_name::Vector{String}
end


function Print_Symmetry(symmetry::Symmetry)
    println("<Print_Symmetry>")
    println("\tspacegroup_number : $(symmetry.spacegroup_number)")
    println("\thall_number : $(symmetry.hall_number)")
    println("\tinternational_symbol : $(symmetry.international_symbol)")
    println("\thall_symbol : $(symmetry.hall_symbol)")
    # println("\tchoice : $(symmetry.choice)")
    println("\torigin_shift : $(symmetry.origin_shift)")
    println("\tn_atoms : $(symmetry.n_atoms)")
    println("\twyckoffs : $(symmetry.wyckoffs)")
    println("\tsite_symmetry_symbols : $(symmetry.site_symmetry_symbols)")
    println("\tequivalent_atoms : $(symmetry.equivalent_atoms)")
    # println("\tcrystallographic_orbits : $(symmetry.crystallographic_orbits)")
    # println("\tprimitive_lattice : $(transpose(symmetry.primitive_lattice))")
    # println("\tmapping_to_primitive : $(symmetry.mapping_to_primitive)")
    # println("\tn_std_atoms : $(symmetry.n_std_atoms)")
    # println("\tstd_lattice : $(symmetry.std_lattice)")
    # println("\tstd_types : $(symmetry.std_types)")
    # println("\tstd_positions : $(symmetry.std_positions)")
    # println("\tstd_rotation_matrix : $(symmetry.std_rotation_matrix)")
    # println("\tstd_mapping_to_primitive : $(symmetry.std_mapping_to_primitive)")
    println("\tpointgroup_symbol : $(symmetry.pointgroup_symbol)")
end


function Get_Symmetry_Spglib(Latvecs::Matrix{Float64}, Gxyz_frac::Vector{Vector{Float64}}, atom2spe; tol_symmetry=1e-5)

    lattice = Spglib.Lattice(Latvecs[1,:],Latvecs[2,:],Latvecs[3,:])
    cell = Spglib.Cell(lattice, Gxyz_frac, atom2spe)

    dataset = Spglib.get_dataset(cell, tol_symmetry)
           
    
    return Symmetry(
            dataset.spacegroup_number, dataset.hall_number, dataset.international_symbol, dataset.hall_symbol,
            dataset.choice, dataset.transformation_matrix, dataset.origin_shift,
            dataset.n_atoms, dataset.wyckoffs, dataset.site_symmetry_symbols, dataset.equivalent_atoms,
            dataset.crystallographic_orbits,
            dataset.primitive_lattice, dataset.mapping_to_primitive,
            dataset.n_std_atoms, dataset.std_lattice,
            dataset.std_types, dataset.std_positions, dataset.std_rotation_matrix, dataset.std_mapping_to_primitive,
            dataset.pointgroup_symbol
    )
end