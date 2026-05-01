struct System_Grid
    CpyCell::Int32
    Natom::Int32
    atom2spe::Vector{Int32}
    Latvecs::Matrix{Float64}
    atv::Vector{Vector{Float64}}
    atv_ijk::Vector{Vector{Int32}}
    Gxyz::Vector{Vector{Float64}}
    GridVol::Float64
    Grid_Origin::Vector{Float64}
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    ncn::Vector{Vector{Int32}}
    Dis::Vector{Vector{Float64}}
    RMI::Vector{Vector{Vector{Int32}}}
    Total_Hsize::Int32
    MPI_Hsize::Vector{Int32}
    MPHks::Vector{Int32}
    Nloop::Int32
    MPI_size::Int32
    MPI_atom::Vector{Int32}
    MPI_FNAN::Vector{Int32}
    MPI_natn::Vector{Int32}
    MPI_ncn::Vector{Int32}
    MPI_Dis::Vector{Float64}
    MPI_RMI::Vector{Vector{Int32}}
    Atom_Cut1::Vector{Float64}
    Total_NumOrbs::Vector{Int32}
    MP::Vector{Int32}
    Ngrid::Tuple{Int32,Int32,Int32}
end


struct UCell
    system_grid::System_Grid
    GridN_Atom::Vector{Int32}
    GridListAtom::Vector{Vector{Int32}}
    CellListAtom::Vector{Vector{Int32}} 
    MPI_NumOLG::Vector{Int32}
    MPI_GListTAtoms1::Vector{Vector{Int32}}
    MPI_GListTAtoms2::Vector{Vector{Int32}}
    Ngrid::Tuple{Int32,Int32,Int32}
end



function split_system_grid(Natom, FNAN, natn, ncn, Dis, RMI, Total_NumOrbs)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Nloop = sum(FNAN.+1)
    OneD2atom = zeros(Int32, Nloop)
    OneD2FNAN = zeros(Int32, Nloop)
    OneD2natn = zeros(Int32, Nloop)
    OneD2ncn = zeros(Int32, Nloop)
    OneD2Dis = zeros(Float64, Nloop)
    OneD2RMI = Vector{Vector{Int32}}(undef, Nloop)

    count = 1
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        OneD2atom[count] = atom
        OneD2FNAN[count] = Rn
        OneD2natn[count] = natn[atom][Rn]
        OneD2ncn[count] = ncn[atom][Rn]
        OneD2Dis[count] = Dis[atom][Rn]


        OneD2RMI[count] = zeros(Int32, FNAN[atom]+1)
        for Rm = 1:FNAN[atom]+1
            OneD2RMI[count][Rm] = RMI[atom][Rn][Rm]
        end

        count += 1
    end


    

    myrange = split_evenly(1:Nloop, nprocs)
    MPI_size = length(myrange[myrank+1])

    MPI_atom = OneD2atom[myrange[myrank+1]]
    MPI_FNAN = OneD2FNAN[myrange[myrank+1]]
    MPI_natn = OneD2natn[myrange[myrank+1]]
    MPI_ncn = OneD2ncn[myrange[myrank+1]]
    MPI_Dis = OneD2Dis[myrange[myrank+1]]
    MPI_RMI = OneD2RMI[myrange[myrank+1]]

    
    myHsize = 0
    for loop = 1:MPI_size, _ = 1:Total_NumOrbs[MPI_atom[loop]], _ = 1:Total_NumOrbs[MPI_natn[loop]]
        myHsize += 1
    end
    Total_Hsize = MPI.Allreduce(myHsize, MPI.SUM, comm)
    MPI.Barrier(comm)

    
	_counts1 = zeros(Int32, nprocs)
    for id = 1:nprocs
        if id-1 == myrank
            _counts1[id] = myHsize
        end
        MPI.Barrier(comm)
    end
    MPI.Barrier(comm)

    MPI_Hsize = zeros(Int32, nprocs)
    MPI.Allreduce!(_counts1,MPI_Hsize,nprocs,MPI.SUM,comm) 
    MPI.Barrier(comm)


    MPHks = zeros(Int32, nprocs)
    Sum = 0
    for id = 1:nprocs
        MPHks[id] = Sum
        Sum += MPI_Hsize[id]
    end


    return Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Dis, MPI_RMI 
end


@timeit timer "UCell" function UCell(Latvecs, Natom, atom2spe, Gxyz, Atom_Cut1, Ngrid, Grid_Origin; Total_NumOrbs=nothing, verbosity=1)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    CpyCell, FNAN, natn, ncn, Dis = Get_FNAN(Latvecs, Natom, Gxyz, Atom_Cut1)


    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)
    Generation_ATV_ijk!(CpyCell, atv_ijk)

    RMI = Get_RMI(Natom, CpyCell, FNAN, natn, ncn)
    if myrank == 0
        system = Check_system(FNAN, ncn, atv_ijk)
        if verbosity >= 1
            println("<Check_System> The system is $system.")
        end
    end
    

    if !isnothing(Total_NumOrbs)
        MP = zeros(Int32, Natom+1)
        Sum = 0
        for atom = 1:Natom
            MP[atom+1] = Sum + Total_NumOrbs[atom]
            Sum += Total_NumOrbs[atom]
        end
    else
        Total_NumOrbs = zeros(Int32, Natom)
        MP = zeros(Int32, Natom+1)
    end

    GridVol = abs(det(Latvecs))/prod(Ngrid)



    # split element for MPI
    # one dimensionalization Natom, FNAN, natn, ncn, Dis, RMI
    Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Dis, MPI_RMI = split_system_grid(Natom, FNAN, natn, ncn, Dis, RMI, Total_NumOrbs)
    


    system_grid = System_Grid(CpyCell, Natom, atom2spe, Latvecs,
                              atv, atv_ijk, Gxyz, GridVol, Grid_Origin,
                              FNAN, natn, ncn, Dis, RMI, 
                              Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, 
                              MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Dis, MPI_RMI,
                              Atom_Cut1, Total_NumOrbs, MP, Ngrid)
    

    GridN_Atom, GridListAtom, CellListAtom = Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atom_Cut1, Ngrid, Grid_Origin)
    
    
    MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2 = Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)


    return UCell(system_grid, 
                 GridN_Atom, GridListAtom, CellListAtom,
                 MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2,
                 Ngrid)
end


function Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atom_Cut1, Ngrid, Grid_Origin)
    
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int32, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    Generation_ATV!(CpyCell, Latvecs, atv, atv_ijk, ratv)


    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid[1]
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid[2]
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid[3]


    # reciprocal grid 
    gRecvecs = 2*pi*inv(gLatvecs')



    # if system is non-periodic when CellListAtom nouse
    GridN_Atom = zeros(Int32, Natom)
    GridListAtom = Vector{Vector{Int32}}(undef, Natom)
    CellListAtom = Vector{Vector{Int32}}(undef, Natom)


    nmin = zeros(Int32, 3)
    nmax = zeros(Int32, 3)
    Cxyz = zeros(Float64, 3)
    NOC = zeros(Int32, 4)

    for atom = 1:Natom

        rcut = Atom_Cut1[atom] + 0.5

        for k = 1:3
            if k == 1
                i = 2
                j = 3
            elseif k == 2
                i = 3
                j = 1
            elseif k == 3
                i = 1
                j = 2
            end

            b = Latvecs[i,:]
            c = Latvecs[j,:]

            v = cross(b, c)
            coef = 1/norm(v)
            v = coef*v

            Cx = Gxyz[atom][1] + rcut*v[1] - Grid_Origin[1]
            Cy = Gxyz[atom][2] + rcut*v[2] - Grid_Origin[2]
            Cz = Gxyz[atom][3] + rcut*v[3] - Grid_Origin[3]
            nmax[k] = trunc(Int32, (Cx*gRecvecs[k,1] + Cy*gRecvecs[k,2] + Cz*gRecvecs[k,3])*0.5/pi)

            Cx = Gxyz[atom][1] - rcut*v[1] - Grid_Origin[1]
            Cy = Gxyz[atom][2] - rcut*v[2] - Grid_Origin[2]
            Cz = Gxyz[atom][3] - rcut*v[3] - Grid_Origin[3]
            nmin[k] = trunc(Int32, (Cx*gRecvecs[k,1] + Cy*gRecvecs[k,2] + Cz*gRecvecs[k,3])*0.5/pi)

            if nmax[k] < nmin[k]
                nmin[k], nmax[k] = nmax[k], nmin[k]
            end
        end


        Np = floor(Int32, prod(nmax.-nmin.+1)*3/2)

        Nct = 0
        rcut = Atom_Cut1[atom]

        tmp_GridListAtom = zeros(Int32, Np)
        tmp_CellListAtom = zeros(Int32, Np)

        for n1 = nmin[1]:nmax[1], n2 = nmin[2]:nmax[2], n3 = nmin[3]:nmax[3]

            Find_CGrids!(NOC, Cxyz, CpyCell, Ngrid, n1, n2, n3, atv, ratv, gLatvecs, Grid_Origin)

            Rn = NOC[1]
            l1 = NOC[2]
            l2 = NOC[3]
            l3 = NOC[4]
            N = l1*Ngrid[2]*Ngrid[3] + l2*Ngrid[3] + l3

            dx = Cxyz[1] - Gxyz[atom][1]
            dy = Cxyz[2] - Gxyz[atom][2]
            dz = Cxyz[3] - Gxyz[atom][3]

            R = sqrt(dx^2 + dy^2 + dz^2)

            if R <= rcut
                Nct = Nct + 1
                tmp_GridListAtom[Nct] = N
                tmp_CellListAtom[Nct] = Rn
            end
        end

        GridN_Atom[atom] = Nct
        GridListAtom[atom] = Vector{Int32}(undef, Nct)
        CellListAtom[atom] = Vector{Int32}(undef, Nct)
        for Nc = 1:Nct
            GridListAtom[atom][Nc] = tmp_GridListAtom[Nc]
            CellListAtom[atom][Nc] = tmp_CellListAtom[Nc]
        end
    end


    for atom = 1:Natom
        Grid_sort = sortperm(GridListAtom[atom])
        @. GridListAtom[atom] = GridListAtom[atom][Grid_sort]
        @. CellListAtom[atom] = CellListAtom[atom][Grid_sort]
    end

    

    return GridN_Atom, GridListAtom, CellListAtom
end



"""
    atomoverlap_grid = AtomOverlap_Grid(...)

Create an instance of `AtomOverlap_Grid`.

Mandatory arguments:

- `system_grid`: an instance of `System_Grid`
- `atoms_grid`: an instance of `Atoms_Grid`
"""
function Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)
    
    ratv = zeros(Int64, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)
    Generation_RATV!(CpyCell, ratv)


    # Find overlap grids between two orbitals
    # GListTAtoms0, GListTAtoms1, GListTAtoms2
    MPI_NumOLG = Vector{Int32}(undef, MPI_size)
    MPI_GListTAtoms1 = Vector{Vector{Int32}}(undef, MPI_size)
    MPI_GListTAtoms2 = Vector{Vector{Int32}}(undef, MPI_size)
    
    
    TAtoms1 = zeros(Int32, maximum(GridN_Atom))
    TAtoms2 = zeros(Int32, maximum(GridN_Atom))


    for loop = 1:MPI_size

        fill!(TAtoms1, 0.0)
        fill!(TAtoms2, 0.0)

        atom = MPI_atom[loop]
        jatom = MPI_natn[loop]
        cell = MPI_ncn[loop]

        l1 = atv_ijk[cell+1][1]
        l2 = atv_ijk[cell+1][2]
        l3 = atv_ijk[cell+1][3]

        Nog = -1
        Nc = 0

        for Nh = 1:GridN_Atom[jatom]

            GNh = GridListAtom[jatom][Nh]
            GRh = CellListAtom[jatom][Nh]

            ll1 = atv_ijk[GRh+1][1]
            ll2 = atv_ijk[GRh+1][2]
            ll3 = atv_ijk[GRh+1][3]
                    
            lll1 = l1 + ll1
            lll2 = l2 + ll2
            lll3 = l3 + ll3
                    
            if GridListAtom[atom][1] <= GNh

                if GNh == 0
                    Nc = 0
                else
                    while GNh <= GridListAtom[atom][Nc+1] && Nc ≠ 0
                        Nc = Nc - 10
                        if Nc < 0
                            Nc = 0
                        end
                    end
                end

                # find whether there is the overlapping or not
                if abs(lll1)<=CpyCell && abs(lll2)<=CpyCell && abs(lll3)<=CpyCell
                            
                    GRh1 = ratv[lll1+CpyCell+1,lll2+CpyCell+1,lll3+CpyCell+1]

                    po = 0

                    while po == 0 && Nc<GridN_Atom[atom]

                        GNc = GridListAtom[atom][Nc+1]
                        GRc = CellListAtom[atom][Nc+1]

                        if GNc==GNh && GRc==GRh1
                            Nog += 1

                            TAtoms1[Nog+1] = Nc
                            TAtoms2[Nog+1] = Nh-1

                            po = 1
                        elseif GNh < GNc
                            po = 1
                        end

                        Nc += 1
                    end

                    Nc -= 1

                    if Nc < 0
                         Nc = 0
                    end
                end
            end

            MPI_NumOLG[loop] = Nog + 1
        end

        MPI_GListTAtoms1[loop] = Vector{Int32}(undef, MPI_NumOLG[loop])
        MPI_GListTAtoms2[loop] = Vector{Int32}(undef, MPI_NumOLG[loop])

        for Nog = 1:MPI_NumOLG[loop]
            MPI_GListTAtoms1[loop][Nog] = TAtoms1[Nog]
            MPI_GListTAtoms2[loop][Nog] = TAtoms2[Nog]
        end
    end
    

    return MPI_NumOLG, MPI_GListTAtoms1, MPI_GListTAtoms2
end


# get the FNAN, natn, ncn, Dis
function Get_FNAN(LatVecs, Natom, Gxyz, Atom_Cut1 )
    
    po = 0
    CpyCell = 0
    TFNAN = 0
    count = 1

    countmax = 7

    # calc CpyCell
    while po == 0 && count < countmax

        CpyCell = CpyCell + 1
        atv, TCpyCell = Set_Periodic(LatVecs, CpyCell)

        TFNAN_temp = TFNAN

        FNAN = Estimate_Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)
        TFNAN = sum(FNAN)

        if TFNAN == TFNAN_temp
            po = 1
        end

        count += 1
    end

    if count == countmax || po == 0
        println("count = $(count), po = $(po)")
        error("please check Get_FNAN")
    end

    atv, TCpyCell = Set_Periodic(LatVecs, CpyCell)
    FNAN, natn, ncn, Dis = Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)
    
    return CpyCell, FNAN, natn, ncn, Dis
end


function Set_Periodic(Latvecs, CpyCell)
    
    TN = (2*CpyCell + 1)^3 - 1

    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)

    return atv, TN
end


function Estimate_Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)
    
    FNAN = zeros(Int32, Natom)

    for atom = 1:Natom

        rcutA = Atom_Cut1[atom]
        FNAN[atom] = 0

        for jatom = 1:Natom

            rcutB = Atom_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell

                if atom == jatom && iszero(Rn)
                    continue
                else
                    dx = abs(Gxyz[atom][1] - Gxyz[jatom][1] - atv[Rn+1][1])
                    dy = abs(Gxyz[atom][2] - Gxyz[jatom][2] - atv[Rn+1][2])
                    dz = abs(Gxyz[atom][3] - Gxyz[jatom][3] - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut

                        r = sqrt(dx^2 + dy^2 + dz^2)

                        if r <= rcut
                            FNAN[atom] = FNAN[atom] + 1
                        end
                    end
                end
            end
        end
    end

    return FNAN
end



function Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)

    Max_FNAN = maximum(Estimate_Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell))
    
    FNAN = zeros(Int32, Natom)

    natn_atom = zeros(Int32, Max_FNAN+1)
    ncn_atom = zeros(Int32, Max_FNAN+1)
    Dis_atom = zeros(Float64, Max_FNAN+1)

    natn = Vector{Vector{Int32}}(undef, Natom)
    ncn = Vector{Vector{Int32}}(undef, Natom)
    Dis = Vector{Vector{Float64}}(undef, Natom)


    for atom = 1:Natom
        
        FNAN[atom] = 0
        rcutA = Atom_Cut1[atom]

        for jatom = 1:Natom
            rcutB = Atom_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell

                if atom == jatom && iszero(Rn)
                    natn_atom[1] = atom
                    ncn_atom[1] = 0
                    Dis_atom[1] = 0.0
                else

                    dx = abs(Gxyz[atom][1] - Gxyz[jatom][1] - atv[Rn+1][1])
                    dy = abs(Gxyz[atom][2] - Gxyz[jatom][2] - atv[Rn+1][2])
                    dz = abs(Gxyz[atom][3] - Gxyz[jatom][3] - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut
                        r = sqrt(dx^2 + dy^2 + dz^2)
                        if r <= rcut
                            FNAN[atom] = FNAN[atom] + 1
                            natn_atom[FNAN[atom]+1] = jatom
                            ncn_atom[FNAN[atom]+1] = Rn
                            Dis_atom[FNAN[atom]+1] = r
                        end
                    end
                end
            end
        end

        natn[atom] = zeros(Int32, FNAN[atom]+1)
        ncn[atom] = zeros(Int32, FNAN[atom]+1)
        Dis[atom] = zeros(Float64, FNAN[atom]+1)
        for k = 1:FNAN[atom]+1
            natn[atom][k] = natn_atom[k]
            ncn[atom][k] = ncn_atom[k]
            Dis[atom][k] = Dis_atom[k]
        end
    end


    return FNAN, natn, ncn, Dis
end



function Get_RMI(Natom, CpyCell, FNAN, natn, ncn)

    RMI = Vector{Vector{Vector{Int32}}}(undef, Natom)
    for atom = 1:Natom
        RMI[atom] = Vector{Vector{Int32}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            RMI[atom][Rn] = zeros(Int32, FNAN[atom]+1)
        end
    end


    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv_ijk[cell] = zeros(Int32, 3)
    end
    ratv = zeros(Int64, 2*CpyCell+4, 2*CpyCell+4, 2*CpyCell+4)

    Generation_ATV_ijk!(CpyCell, atv_ijk)
    Generation_RATV!(CpyCell, ratv)



    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        ig = natn[atom][Rn]
        Rni = ncn[atom][Rn]+1
        for Rm = 1:FNAN[atom]+1

            jg = natn[atom][Rm]
            Rnj = ncn[atom][Rm]+1
            l1 = atv_ijk[Rnj][1] - atv_ijk[Rni][1]
            l2 = atv_ijk[Rnj][2] - atv_ijk[Rni][2]
            l3 = atv_ijk[Rnj][3] - atv_ijk[Rni][3]
            m1 = ifelse(l1<0, -l1, l1)
            m2 = ifelse(l2<0, -l2, l2)
            m3 = ifelse(l3<0, -l3, l3)  

            if m1 <= CpyCell && m2 <= CpyCell && m3 <= CpyCell  
                cell = ratv[l1+CpyCell+1,l2+CpyCell+1,l3+CpyCell+1] 
                k = 0
                po = 0
                RMI[atom][Rn][Rm] = -1  

                while po == 0 && k <= FNAN[ig]
                    if natn[ig][k+1]==jg && ncn[ig][k+1] == cell
                        RMI[atom][Rn][Rm] = k
                        po = 1
                    end
                    k = k + 1
                end
            else
                RMI[atom][Rn][Rm] = -1
            end
        end
    end
                
    return RMI
end


function Calc_Dis(Natom, FNAN, Gxyz, atv, Atom_Cut1, CpyCell)

    TCpyCell = (2*CpyCell+1)^3-1
    Max_FNAN = maximum(FNAN)
    Dis_atom = zeros(Float64, Max_FNAN+1)


    Dis = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Dis[atom] = zeros(Float64, FNAN[atom]+1)
    end

    for atom = 1:Natom    
        tmp = 0
        rcutA = Atom_Cut1[atom]
        for jatom = 1:Natom
            rcutB = Atom_Cut1[jatom]
            rcut = rcutA + rcutB

            for Rn = 0:TCpyCell

                if atom == jatom && iszero(Rn)
                    Dis_atom[1] = 0.0
                else
                    dx = abs(Gxyz[atom][1] - Gxyz[jatom][1] - atv[Rn+1][1])
                    dy = abs(Gxyz[atom][2] - Gxyz[jatom][2] - atv[Rn+1][2])
                    dz = abs(Gxyz[atom][3] - Gxyz[jatom][3] - atv[Rn+1][3])
                    
                    if dx <= rcut && dy <= rcut && dz <= rcut
                        r = sqrt(dx^2 + dy^2 + dz^2)
                        if r <= rcut
                            tmp = tmp + 1
                            Dis_atom[tmp+1] = r
                        end
                    end
                end
            end
        end

        for Rn = 1:FNAN[atom]+1
            Dis[atom][Rn] = Dis_atom[Rn]
        end
    end


    return Dis
end


function Check_system(FNAN, ncn, atv_ijk)

    Natom = length(FNAN)
    po = zeros(Int64, 3)

    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        cell = ncn[atom][Rn]
        if cell ≠ 0
            if atv_ijk[cell+1][1] ≠ 0
                po[1] = 1
            elseif atv_ijk[cell+1][2] ≠ 0
                po[2] = 1
            elseif atv_ijk[cell+1][3] ≠ 0
                po[3] = 1
            end
        end
    end

    num = sum(po)

    if num == 0
        system = "molecule"
    elseif num == 1
        system = "chain"
    elseif num == 2
        system = "slab"
    elseif num == 3
        system = "bluk"
    else
        error("please check Check_System")
    end

    return system
end