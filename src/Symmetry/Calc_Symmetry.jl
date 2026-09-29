struct Symmetry
    # lattice::Vector{Vector{Float64}}
    # positions::Vector{Vector{Float64}}
    # atoms::Vector{Int32}
    # magmoms::Vector{Float64}
    Spglib_type::Type
    SpinPol::String
    Latvecs::Matrix{Float64}
    Gxyz_frac::Vector{Vector{Float64}}
    atom2spe::Vector{Int32}
    Atoms_Nspin::Vector{Vector{Float64}}
    Atoms_Angle::Vector{Vector{Float64}}
    spacegroup::Int32
    transformation_matrix::Matrix{Float64}
    n_operations::Int32
    rotations::Vector{Matrix{Int32}}
    translations::Vector{Vector{Float64}}
    time_reversals::Vector{Int32}
    tol_symmetry::Float64
    # spacegroup_number::Int64
    # hall_number::Int64
    # international_symbol::String
    # hall_symbol::String
    # choice::String
    # transformation_matrix::Matrix{Float64}
    # origin_shift::Vector{Float64}
    # n_atoms::Int64
    # wyckoffs::Vector{Char}
    # site_symmetry_symbols::Vector{String}
    # equivalent_atoms::Vector{Int32}
    # crystallographic_orbits::Vector{Int32}
    # primitive_lattice::Matrix{Float64}
    # mapping_to_primitive::Vector{Int32}
    # n_std_atoms::Int64
    # std_lattice::Matrix{Float64}
    # std_types::Vector{Int32}
    # std_positions::Vector{SVector{3, Float64}}
    # std_rotation_matrix::Matrix{Float64}
    # std_mapping_to_primitive::Vector{Int32}
    # pointgroup_symbol::String
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
    println("\tn_operations : $(symmetry.n_operations)")
    # println("\tspacegroup_number : $(symmetry.spacegroup_number)")
    # println("\thall_number : $(symmetry.hall_number)")
    # println("\tinternational_symbol : $(symmetry.international_symbol)")
    # println("\thall_symbol : $(symmetry.hall_symbol)")
    # println("\tchoice : $(symmetry.choice)")
    # println("\torigin_shift : $(symmetry.origin_shift)")
    # println("\tn_atoms : $(symmetry.n_atoms)")
    # println("\twyckoffs : $(symmetry.wyckoffs)")
    # println("\tsite_symmetry_symbols : $(symmetry.site_symmetry_symbols)")
    # println("\tequivalent_atoms : $(symmetry.equivalent_atoms)")
    # println("\tcrystallographic_orbits : $(symmetry.crystallographic_orbits)")
    # println("\tprimitive_lattice : $(transpose(symmetry.primitive_lattice))")
    # println("\tmapping_to_primitive : $(symmetry.mapping_to_primitive)")
    # println("\tn_std_atoms : $(symmetry.n_std_atoms)")
    # println("\tstd_lattice : $(symmetry.std_lattice)")
    # println("\tstd_types : $(symmetry.std_types)")
    # println("\tstd_positions : $(symmetry.std_positions)")
    # println("\tstd_rotation_matrix : $(symmetry.std_rotation_matrix)")
    # println("\tstd_mapping_to_primitive : $(symmetry.std_mapping_to_primitive)")
    # println("\tpointgroup_symbol : $(symmetry.pointgroup_symbol)")
end


function Set_Spglib_magmoms(SpinPol::String, Natom, Atoms_Nspin::Vector{Vector{Float64}}, Atoms_Angle::Vector{Vector{Float64}})

    magmoms = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        magmoms[atom] = zeros(Float64, 3)
        M = Atoms_Nspin[atom][2] - Atoms_Nspin[atom][1]
        theta, phi = Atoms_Angle[atom]
        theta = deg2rad(theta)
        phi = deg2rad(phi)
        magmoms[atom][1] = M*sin(theta)*cos(phi)
        magmoms[atom][2] = M*sin(theta)*sin(phi)
        magmoms[atom][3] = M*cos(theta)
    end

    return magmoms
end


function Get_Symmetry_Spglib_Col(Latvecs, Gxyz_frac, atom2spe, tol_symmetry)

    Natom = length(Gxyz_frac)
    Atoms_Nspin = [[0.0, 0.0] for atom = 1:Natom]
    Atoms_Angle = [[0.0, 0.0] for atom = 1:Natom]
    cell = Spglib.SpglibCell(transpose(Latvecs), Gxyz_frac, atom2spe)
    dataset = Spglib.get_dataset(cell, tol_symmetry)
    spacegroup = dataset.spacegroup_number

    return Symmetry( 
            typeof(dataset), "off", Latvecs, Gxyz_frac, atom2spe, Atoms_Nspin, Atoms_Angle, 
            spacegroup, dataset.transformation_matrix,
            dataset.n_operations, dataset.rotations, dataset.translations, zeros(Int32, dataset.n_operations),
            tol_symmetry)
    #=
    return Symmetry(
            cell,
            dataset.spacegroup_number, dataset.hall_number, dataset.international_symbol, dataset.hall_symbol,
            dataset.choice, dataset.transformation_matrix, dataset.origin_shift,
            dataset.n_atoms, dataset.wyckoffs, dataset.site_symmetry_symbols, dataset.equivalent_atoms,
            dataset.crystallographic_orbits,
            dataset.primitive_lattice, dataset.mapping_to_primitive,
            dataset.n_std_atoms, dataset.std_lattice,
            dataset.std_types, dataset.std_positions, dataset.std_rotation_matrix, dataset.std_mapping_to_primitive,
            dataset.pointgroup_symbol)=#
end


function Get_Symmetry_Spglib_NonCol(SpinPol::String, Latvecs, Gxyz_frac, atom2spe, Atoms_Nspin, Atoms_Angle, tol_symmetry)

    Natom = length(Gxyz_frac)
    magmoms = Set_Spglib_magmoms(SpinPol, Natom, Atoms_Nspin, Atoms_Angle)
    cell = Spglib.SpglibCell(transpose(Latvecs), Gxyz_frac, atom2spe, magmoms)
    dataset = Spglib.get_magnetic_dataset(cell, tol_symmetry)
    spacegroup = Int(Spglib.get_spacegroup_type(dataset.hall_number).number)

    return Symmetry(
            typeof(dataset), SpinPol, Latvecs, Gxyz_frac, atom2spe, Atoms_Nspin, Atoms_Angle, 
            spacegroup, dataset.transformation_matrix,
            dataset.n_operations, dataset.rotations, dataset.translations, dataset.time_reversals,
            tol_symmetry)
    #=
    return Symmetry(
            cell,
            dataset.spacegroup_number, dataset.hall_number, dataset.international_symbol, dataset.hall_symbol,
            dataset.choice, dataset.transformation_matrix, dataset.origin_shift,
            dataset.n_atoms, dataset.wyckoffs, dataset.site_symmetry_symbols, dataset.equivalent_atoms,
            dataset.crystallographic_orbits,
            dataset.primitive_lattice, dataset.mapping_to_primitive,
            dataset.n_std_atoms, dataset.std_lattice,
            dataset.std_types, dataset.std_positions, dataset.std_rotation_matrix, dataset.std_mapping_to_primitive,
            dataset.pointgroup_symbol)=#
end


@timeit timer "Get_Symmetry" function Get_Symmetry_Spglib(
    SpinPol::String, 
    Latvecs::Matrix{Float64}, 
    Gxyz_frac::Vector{Vector{Float64}}, 
    atom2spe::Vector{<:Integer}, 
    Atoms_Nspin::Vector{Vector{Float64}}, 
    Atoms_Angle::Vector{Vector{Float64}}; 
    tol_symmetry=1e-5)
    
    if SpinPol == "nc"
        return Get_Symmetry_Spglib_NonCol(SpinPol, Latvecs, Gxyz_frac, atom2spe, Atoms_Nspin, Atoms_Angle, tol_symmetry)
    else
        return Get_Symmetry_Spglib_Col(Latvecs, Gxyz_frac, atom2spe, tol_symmetry)
    end
end