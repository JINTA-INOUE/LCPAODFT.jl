struct PackedOrbitalsGrid
    atoms::Vector{Int32}
    data::Vector{Matrix{Float64}}
end

Base.length(grid::PackedOrbitalsGrid) = length(grid.data)

@inline function _local_atom_index(grid::PackedOrbitalsGrid, atom::Integer)
    global_atom = Int32(atom)
    local_atom = searchsortedfirst(grid.atoms, global_atom)
    if local_atom > length(grid.atoms) || grid.atoms[local_atom] != global_atom
        throw(BoundsError(grid, atom))
    end
    return local_atom
end
@inline atom_matrix(grid::PackedOrbitalsGrid, atom::Integer) = @inbounds grid.data[_local_atom_index(grid, atom)]

"""Sorted atoms whose orbital grids are needed by the current MPI rank."""
function _rank_orbital_atoms(system_grid::System_Grid)
    return sort!(unique(vcat(system_grid.MPI_atom, system_grid.MPI_natn)))
end

function _empty_orbitals_grid(atoms::Vector{Int32}, total_orbitals, grid_points)
    data = [zeros(Float64, total_orbitals[atom], grid_points[atom])
            for atom in atoms]
    return PackedOrbitalsGrid(atoms, data)
end

@inline function _set_real_harmonics!(values::Vector{Vector{Float64}},
                                      max_l::Integer,
                                      theta::Float64, phi::Float64)
    max_l <= 3 || error("PAO angular momentum l > 3 is not supported")
    for l = 0:max_l
        if l == 0
            values[1][1] = Ylm_real(0, 0, theta, phi)
        elseif l == 1
            values[2][1] = Ylm_real(1, 1, theta, phi)
            values[2][2] = Ylm_real(1, -1, theta, phi)
            values[2][3] = Ylm_real(1, 0, theta, phi)
        elseif l == 2
            values[3][1] = Ylm_real(2, 0, theta, phi)
            values[3][2] = Ylm_real(2, 2, theta, phi)
            values[3][3] = Ylm_real(2, -2, theta, phi)
            values[3][4] = Ylm_real(2, 1, theta, phi)
            values[3][5] = Ylm_real(2, -1, theta, phi)
        else
            values[4][1] = Ylm_real(3, 0, theta, phi)
            values[4][2] = Ylm_real(3, 1, theta, phi)
            values[4][3] = Ylm_real(3, -1, theta, phi)
            values[4][4] = Ylm_real(3, 2, theta, phi)
            values[4][5] = Ylm_real(3, -2, theta, phi)
            values[4][6] = Ylm_real(3, 3, theta, phi)
            values[4][7] = Ylm_real(3, -3, theta, phi)
        end
    end
    return nothing
end

@timeit timer "Set_Orbitals_Grid" function Set_Orbitals_Grid(
    pao::Vector{PAO}, ucell::UCell)
    atoms = _rank_orbital_atoms(ucell.system_grid)
    orbitals_grid = _empty_orbitals_grid(
        atoms, ucell.system_grid.Total_NumOrbs, ucell.GridN_Atom)
    Set_Orbitals_Grid!(orbitals_grid, pao, ucell)
    return orbitals_grid
end


function Set_Orbitals_Grid!(orbitals_grid::PackedOrbitalsGrid,
                            pao::Vector{PAO}, ucell::UCell)

    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    atv = system_grid.atv
    Grid_Origin = system_grid.Grid_Origin
    Gxyz = system_grid.Gxyz
    atom2spe = system_grid.atom2spe
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom

    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    maxSpe_MaxL_Basis = maximum(species.Spe_MaxL_Basis for species in pao)
    pmax = maximum(maximum(species.Spe_Num_Basis) for species in pao)

    RF = Vector{Vector{Float64}}(undef, maxSpe_MaxL_Basis+1)
    AF = Vector{Vector{Float64}}(undef, maxSpe_MaxL_Basis+1)
    for l = 0:maxSpe_MaxL_Basis
        RF[l+1] = zeros(Float64, pmax)
        AF[l+1] = zeros(Float64, 2*l+1)
    end


    for (local_atom, atom) in pairs(orbitals_grid.atoms)
        atom_orbitals = orbitals_grid.data[local_atom]

        spe = atom2spe[atom]
        Spe_MaxL_Basis = pao[spe].Spe_MaxL_Basis
        Spe_Num_Basis = pao[spe].Spe_Num_Basis
        Spe_Num_Mesh_PAO = pao[spe].Spe_Num_Mesh_PAO
        Spe_PAO_RV = pao[spe].Spe_PAO_RV
        Spe_PAO_RWF = pao[spe].Spe_PAO_RWF

        for xyz = 1:GridN_Atom[atom]

            GNc = GridListAtom[atom][xyz]
            GRc = CellListAtom[atom][xyz]+1

            # GNc = n1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3
            n1 = div(GNc, Ngrid2*Ngrid3)
            n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
            n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3
            
            Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
            Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
            Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

            x = Cx + atv[GRc][1] - Gxyz[atom][1]
            y = Cy + atv[GRc][2] - Gxyz[atom][2]
            z = Cz + atv[GRc][3] - Gxyz[atom][3]
            R, theta, phi = xyz_to_spherical(x,y,z)
            
            po = 0
            mp_min = 1
            mp_max = Spe_Num_Mesh_PAO

            for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
                RF[l+1][p] = 0.0
            end


            if Spe_PAO_RV[end] < R
                po = 1
            elseif R < Spe_PAO_RV[begin]

                m = 5
                rm = Spe_PAO_RV[m]
        
                h1 = Spe_PAO_RV[m-1] - Spe_PAO_RV[m-2]
                h2 = Spe_PAO_RV[m]   - Spe_PAO_RV[m-1]
                h3 = Spe_PAO_RV[m+1] - Spe_PAO_RV[m]
        
                x1 = rm - Spe_PAO_RV[m-1]
                x2 = rm - Spe_PAO_RV[m]
                y1 = x1/h2
                y2 = x2/h2
                y12 = y1*y1
                y22 = y2*y2
        
                dum = h1 + h2
                dum1 = h1/h2/dum
                dum2 = h2/h1/dum
                dum = h2 + h3
                dum3 = h2/h3/dum
                dum4 = h3/h2/dum
        
                for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
        
                    f1 = Spe_PAO_RWF[l+1][p][m-2]
                    f2 = Spe_PAO_RWF[l+1][p][m-1]
                    f3 = Spe_PAO_RWF[l+1][p][m]
                    f4 = Spe_PAO_RWF[l+1][p][m+1]
                
                    dum = f3 - f2
                    g1 = dum*dum1 + (f2-f1)*dum2
                    g2 = (f4-f3)*dum3 + dum*dum4

                    A = 3*f2 + h2*g1
                    B = 2*f2 + h2*g1
                    C = 3*f3 - h2*g2
                    D = 2*f3 - h2*g2
                    E = A + B*y2
                    F = C - D*y1
                    f = y22*E + y12*F
                    df = (2*y2/h2)*E + y22*B/h2 + (2*y1/h2)*F - y12*D/h2
                
                    if l == 0
                        a = 0.0
                        b = 0.5*df/rm
                        c = 0.0
                        d = f - b*rm^2
                    elseif l == 1
                        a = (rm*df - f)/(2*rm^3)
                        b = 0.0
                        c = df - 3*a*rm^2
                        d = 0.0
                    else
                        b = (3*f - rm*df)/rm^2
                        a = (f - b*rm^2)/rm^3
                        c = 0.0
                        d = 0.0
                    end
                
                    RF[l+1][p] = a*R^3 + b*R^2 + c*R + d
                end
            else
                while (mp_max-mp_min) ≠ 1
                    m = div(mp_min + mp_max, 2)
                    if (Spe_PAO_RV[m]<R)
                        mp_min = m
                    else 
                        mp_max = m
                    end
                end

                m = mp_max
            
                h1 = Spe_PAO_RV[m-1] - Spe_PAO_RV[m-2]
                h2 = Spe_PAO_RV[m]   - Spe_PAO_RV[m-1]
                h3 = Spe_PAO_RV[m+1] - Spe_PAO_RV[m]
            
                x1 = R - Spe_PAO_RV[m-1]
                x2 = R - Spe_PAO_RV[m]
                y1 = x1/h2
                y2 = x2/h2
                y12 = y1*y1
                y22 = y2*y2
            
                dum = h1 + h2
                dum1 = h1/h2/dum
                dum2 = h2/h1/dum
                dum = h2 + h3
                dum3 = h2/h3/dum
                dum4 = h3/h2/dum
          
                for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
            
                    f1 = Spe_PAO_RWF[l+1][p][m-2]
                    f2 = Spe_PAO_RWF[l+1][p][m-1]
                    f3 = Spe_PAO_RWF[l+1][p][m]
                    f4 = Spe_PAO_RWF[l+1][p][m+1]
                
                    if m == 1
                        h1 = -(h2+h3)
                        f1 = f4
                    elseif m == (Spe_Num_Mesh_PAO-1)
                        h3 = -(h1+h2)
                        f4 = f1
                    end
                
                    dum = f3 - f2
                    g1 = dum*dum1 + (f2-f1)*dum2
                    g2 = (f4-f3)*dum3 + dum*dum4
            
                    f = (y22*(3*f2 + h2*g1 + (2*f2 + h2*g1)*y2)+ y12*(3*f3 - h2*g2 - (2*f3 - h2*g2)*y1))
            
                    RF[l+1][p] = f
                end
            end
            
            # Real spherical harmonics in the PAO orbital ordering.
            if po == 0
                _set_real_harmonics!(AF, Spe_MaxL_Basis, theta, phi)
            end

            ist = 0
            @inbounds for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1], m = 1:2*l+1
                ist += 1
                atom_orbitals[ist, xyz] = RF[l+1][p]*AF[l+1][m]
            end
        end
    end
end
