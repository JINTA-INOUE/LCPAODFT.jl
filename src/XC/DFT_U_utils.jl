"""
    Set_Hub_U_Basis(...)
    
    # Examples
    ```
        Natom = 4
        Atom_MaxL_Basis = Int32[2, 2, 2, 2]
        Atom_Num_Basis = [[2, 1, 1], [2, 1, 1], [1, 1, 1], [1, 1, 1]]
        Hub_U_Atom = [[0.0,0.0,0.0,4.0], [0.0,0.0,0.0,4.0], [0.0,0.0,0.0], [0.0,0.0,0.0]]
        Hub_U_Basis = Set_Hub_U_Basis(Natom, Atom_MaxL_Basis, Atom_Num_Basis, Hub_U_Atom)
    ```
"""
function Set_Hub_U_Basis(Natom, Atom_MaxL_Basis, Atom_Num_Basis, Hub_U_Atom)

    if Natom ≠ length(Hub_U_Atom)
        error("please check Hubbard U values.")
    end

    for atom = 1:Natom
        if length(Hub_U_Atom[atom]) ≠ sum(Atom_Num_Basis[atom])
            @show Hub_U_Atom
            @show Atom_MaxL_Basis
            @show Atom_Num_Basis
            error("please check Set_Hub_U_Basis.")
        end
    end

    Hub_U_Basis = Vector{Vector{Vector{Float64}}}(undef, Natom)
    for atom = 1:Natom
        Hub_U_Basis[atom] = Vector{Vector{Float64}}(undef, Atom_MaxL_Basis[atom]+1)
        for l = 0:Atom_MaxL_Basis[atom]
            Hub_U_Basis[atom][l+1] = zeros(Float64, Atom_Num_Basis[atom][l+1])
        end
    end


    for atom = 1:Natom
        ist = 0
        for l = 0:Atom_MaxL_Basis[atom], p = 1:Atom_Num_Basis[atom][l+1]
            ist += 1
            Hub_U_Basis[atom][l+1][p] = Hub_U_Atom[atom][ist]/eV2Hartree
        end
    end


    return Hub_U_Basis
end