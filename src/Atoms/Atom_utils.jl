function Get_Atoms_data(orbitals::Vector{String})
    
    Natom = length(orbitals)
    atoms_Symbol = Vector{String}(undef, Natom)
    atoms_cutoff = Vector{Float64}(undef, Natom)
    atoms_orb = Vector{String}(undef, Natom)
    atoms_Extrasymbol = Vector{String}(undef, Natom)
    
    for atom = 1:Natom
        basefile, atoms_orb[atom] = split(orbitals[atom], "-")
        
        if length(basefile) == 4    # temp = "C5.0" case
            atoms_Symbol[atom] = basefile[1:1]
            atoms_cutoff[atom] = parse(Float64, basefile[2:end])
            atoms_Extrasymbol[atom] = ""
        elseif length(basefile) == 5    # temp = "Si5.0" or "S10.0" case
            cutoff = parse(Float64, basefile[3:end])
            if cutoff ∈ (0.0, 1.0, 2.0, 3.0, 4.0)
                atoms_Symbol[atom] = basefile[1:1]
                atoms_cutoff[atom] = parse(Float64, basefile[2:end])
                atoms_Extrasymbol[atom] = ""
            else
                atoms_Symbol[atom] = basefile[1:2]
                atoms_cutoff[atom] = cutoff
                atoms_Extrasymbol[atom] = ""
            end
        elseif length(basefile) == 6    # temp = "Cu6.0H" or "Ag11.0" case
            Atom_Symbol = basefile[1:2]
            cutoff = parse(Float64, basefile[3:5])
            if basefile[end] == 'H'
                atoms_Symbol[atom] = Atom_Symbol
                atoms_cutoff[atom] = cutoff
                atoms_Extrasymbol[atom] = "H"
            elseif basefile[end] == 'S'
                atoms_Symbol[atom] = Atom_Symbol
                atoms_cutoff[atom] = cutoff
                atoms_Extrasymbol[atom] = "S"
            else
                atoms_Symbol[atom] = Atom_Symbol
                atoms_cutoff[atom] = parse(Float64, basefile[3:end])
                atoms_Extrasymbol[atom] = ""
            end
        elseif length(basefile) == 7    # temp = "Ni10.0H" case
            Atom_Symbol = basefile[1:2]
            cutoff = parse(Float64, basefile[3:5])
            if basefile[end] == 'H'
                atoms_Symbol[atom] = Atom_Symbol
                atoms_cutoff[atom] = cutoff
                atoms_Extrasymbol[atom] = "H"
            elseif basefile[end] == 'S'
                atoms_Symbol[atom] = Atom_Symbol
                atoms_cutoff[atom] = cutoff
                atoms_Extrasymbol[atom] = "S"
            else
                error("please check")
            end
        else
            println("input = ", basefile)
            error("please check input Atoms orbital")
        end
    end
    

    return atoms_Symbol, atoms_cutoff, atoms_orb, atoms_Extrasymbol
end
