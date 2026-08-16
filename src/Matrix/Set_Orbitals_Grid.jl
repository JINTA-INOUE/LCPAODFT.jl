struct PackedOrbitalsGrid
    data::Vector{Matrix{Float64}}
end

struct _PackedAtomOrbitals{M<:AbstractMatrix{Float64}}
    data::M
end

struct _PackedOrbitalColumn{M<:AbstractMatrix{Float64}}
    data::M
    column::Int
end

Base.length(grid::PackedOrbitalsGrid) = length(grid.data)
@inline Base.getindex(grid::PackedOrbitalsGrid, atom::Integer) = _PackedAtomOrbitals(@inbounds grid.data[Int(atom)])
Base.length(atom::_PackedAtomOrbitals) = size(atom.data, 2)
@inline Base.getindex(atom::_PackedAtomOrbitals, column::Integer) = _PackedOrbitalColumn(atom.data, Int(column))
Base.length(column::_PackedOrbitalColumn) = size(column.data, 1)
@inline Base.getindex(column::_PackedOrbitalColumn, orbital::Integer) = @inbounds column.data[Int(orbital), column.column]
@inline function Base.setindex!(column::_PackedOrbitalColumn, value, orbital::Integer)
    @inbounds column.data[Int(orbital), column.column] = value
    return value
end


@timeit timer "Set_Orbitals_Grid" function Set_Orbitals_Grid(pao::Vector{PAO}, ucell::UCell)
    
    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    GridN_Atom = ucell.GridN_Atom
    
    packed = Vector{Matrix{Float64}}(undef, Natom)
    for atom = 1:Natom
        packed[atom] = zeros(Float64, Total_NumOrbs[atom], GridN_Atom[atom])
    end
    Orbs_Grid = PackedOrbitalsGrid(packed)
    Set_Orbitals_Grid!(Orbs_Grid, pao, ucell)

    
    return Orbs_Grid
end


function Set_Orbitals_Grid!(Orbs_Grid, pao::Vector{PAO}, ucell::UCell)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Nspecies = length(pao)
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


    pmax = 4
    maxSpe_MaxL_Basis = 0
    for spe = 1:Nspecies
        maxSpe_MaxL_Basis = max(maxSpe_MaxL_Basis, pao[spe].Spe_MaxL_Basis)
    end

    RF = Vector{Vector{Float64}}(undef, maxSpe_MaxL_Basis+1)
    AF = Vector{Vector{Float64}}(undef, maxSpe_MaxL_Basis+1)
    for l = 0:maxSpe_MaxL_Basis
        RF[l+1] = zeros(Float64, pmax)
        AF[l+1] = zeros(Float64, 2*l+1)
    end


    for atom = 1:Natom

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
            
            # multiple real spherical harmics function Ylm
            if po == 0
                for l = 0:Spe_MaxL_Basis
                    if l == 0
                        AF[1][1] = Ylm_real(0,0,theta,phi)      # s
                    elseif l == 1
                        AF[2][1] = Ylm_real(1,1,theta,phi)      # px
                        AF[2][2] = Ylm_real(1,-1,theta,phi)     # py
                        AF[2][3] = Ylm_real(1,0,theta,phi)      # pz
                    elseif l == 2
                        AF[3][1] = Ylm_real(2,0,theta,phi)      # dz^2
                        AF[3][2] = Ylm_real(2,2,theta,phi)      # dx^2-y^2
                        AF[3][3] = Ylm_real(2,-2,theta,phi)     # dxy
                        AF[3][4] = Ylm_real(2,1,theta,phi)      # dxz
                        AF[3][5] = Ylm_real(2,-1,theta,phi)     # dyz
                    elseif l == 3
                        AF[4][1] = Ylm_real(3,0,theta,phi)      # z^3
                        AF[4][2] = Ylm_real(3,1,theta,phi)      # xz^2
                        AF[4][3] = Ylm_real(3,-1,theta,phi)     # yz^2
                        AF[4][4] = Ylm_real(3,2,theta,phi)      # z(x^2-y^2)
                        AF[4][5] = Ylm_real(3,-2,theta,phi)     # xyz
                        AF[4][6] = Ylm_real(3,3,theta,phi)      # x(x^2-3y^2)
                        AF[4][7] = Ylm_real(3,-3,theta,phi)     # y(3x^2-y^2)
                    else
                        error("not support PAO l > 4")
                    end
                end
            end

            ist = 0
            @inbounds for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1], m = 1:2*l+1
                ist += 1
                Orbs_Grid[atom][xyz][ist] = RF[l+1][p]*AF[l+1][m]
            end
        end
    end
end
