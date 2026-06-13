mutable struct Energy
    Eele::Float64
    Ekin::Float64
    EH0::Float64
    EH1::Float64
    Ena::Float64
    Enl::Float64
    Exc0::Float64
    Exc1::Float64
    Ecore::Float64
    EHub::Float64
    Eef::Float64
    Etot::Float64
    ChemP::Float64
end


"""
    energy = Energy()

Create an instance of `Energy` (Hartree unit).

- `Eele`:   band energy  
- `Ekin`:   kinetic energy  
- `EH0`:    electric part of screened Coulomb energy  
- `EH1`:    difference electron-electron Coulomb energy  
- `Ena`:    neutral atom potential energy  
- `Enl`:    non-local potential energy  
- `Exc0`:   exchange-correlation energy for alpha spin  
- `Exc1`:   exchange-correlation energy for beta spin  
- `Ecore`:  core-core Coulomb energy  
- `EHub`:   LDA+U energy
- `Etot`:   Total energy (= Ekin + EH0 + EH1 + Ena + Enl + Exc0 + Exc1 + Ecore)
- `ChemP`:  Chemical potential energy  
"""
function Init_Energy(; Eele=0.0, ChemP=0.0)
    return Energy(Eele, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ChemP)    
end



function Print_Energy(energy::Energy)
    @printf("  Eele  =     %5.12f\n", energy.Eele)
    @printf("  Ekin  =     %5.12f\n", energy.Ekin)
    @printf("  EH0   =     %5.12f\n", energy.EH0)
    @printf("  EH1   =     %5.12f\n", energy.EH1)
    @printf("  Ena   =     %5.12f\n", energy.Ena)
    @printf("  Enl   =     %5.12f\n", energy.Enl)
    @printf("  Exc0  =     %5.12f\n", energy.Exc0)
    @printf("  Exc1  =     %5.12f\n", energy.Exc1)
    @printf("  Ecore =     %5.12f\n", energy.Ecore)
    @printf("  EHub  =     %5.12f\n", energy.EHub)
    @printf("  Eef   =     %5.12f\n", energy.Eef)
    @printf("  Etot  =     %5.12f\n", energy.Etot)
    @printf("  ChemP =     %5.12f\n", energy.ChemP)
end


"""
    Total_Energy

Mandatory arguments:

- `DM`: Density matrix real part
- `iDM`: Density matrix imaginary part
- `ADensity_Grid`: atomic electron Density in real space
- `PCCDensity_Grid`: pcc electron Density in real space
- `Density_Grid`: electron Density in real space
- `dVHart_Grid`: Hartree potetial in real space
- `Ham`: an instance of `Hamiltonian`
- `system_grid`: an instance of `System_Grid`
- `pao`: an instance of `PAO`
- `pspot`: an instance of `Pspot`
"""
@timeit timer "Total_Energy" function Total_Energy!(
    energy::Energy, force::Force,
    DM, iDM, 
    ADensity_Grid, PCCDensity_Grid, Density_Grid, dVHart_Grid, 
    Ham::Hamiltonian, system_grid::System_Grid, 
    pao::Vector{PAO}, pspot::Vector{Pspot})

    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    GridVol = system_grid.GridVol
    SpinPol = Ham.SpinPol
    Hkin = Ham.Hkin
	HVNA = Ham.HVNA
	HNL = Ham.HNL
    iHNL = Ham.iHNL

    Core_Charge = zeros(Float64, Natom)
	for atom = 1:Natom
		spe = atom2spe[atom]
        Core_Charge[atom] = pspot[spe].Spe_Core_Charge
	end


    Enl = Calc_Enl(SpinPol, DM, iDM, HNL, iHNL)
    Ekin = Calc_Ekin(SpinPol, DM, Hkin)
    Ena = Calc_Ena(SpinPol, DM, HVNA)
    EH0, EH0Force = Calc_EH0( pao, pspot, system_grid)
    EH1 = Calc_EH1(SpinPol, GridVol, ADensity_Grid, Density_Grid, dVHart_Grid)
    Exc, ExcForce = Calc_EXC(SpinPol, pao, pspot, system_grid, ADensity_Grid, PCCDensity_Grid, Density_Grid)
    Ecore = Calc_Ecore(system_grid, Core_Charge)

    Etot = Ecore + Ekin + Ena + Enl + EH0 + EH1 + Exc[1] + Exc[2]

    energy.Ekin = Ekin
    energy.EH0 = EH0
    energy.EH1 = EH1
    energy.Ena = Ena
    energy.Enl = Enl
    energy.Exc0 = Exc[1]
    energy.Exc1 = Exc[2]
    energy.Ecore = Ecore
    energy.Etot = Etot

    force.EH0Force = EH0Force
    force.ExcForce = ExcForce
end


function Calc_Ecore(system_grid::System_Grid, Core_Charge)
    
    comm = MPI.COMM_WORLD

    MPI_size = system_grid.MPI_size
    MPI_atom = system_grid.MPI_atom
    MPI_natn = system_grid.MPI_natn
    MPI_FNAN = system_grid.MPI_FNAN
    Dis = system_grid.Dis

    Ecore = 0.0
    for loop = 1:MPI_size
        Rn = MPI_FNAN[loop]
        if Rn ≠ 1
            atom = MPI_atom[loop]
            Rn = MPI_FNAN[loop]
            jatom = MPI_natn[loop]
            dis = Dis[atom][Rn]
            Zc = Core_Charge[atom]
            Zh = Core_Charge[jatom]
            
            r = ifelse(dis<1e-10, 1e-10, dis) 

            Ecore += Zc*Zh/r
        end
    end
    Ecore = MPI.Allreduce(Ecore, MPI.SUM, comm)

    
    return 0.5*Ecore
end


function Calc_EH0(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)


    Gxyz = system_grid.Gxyz
    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv = system_grid.atv
    Dis = system_grid.Dis
    RMI = system_grid.RMI
    MPI_size = system_grid.MPI_size
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_ncn = system_grid.MPI_ncn
    

    Nspecies = length(pao)
    Scale_Grid_Ecut = 16.0 * 600.0
    Max_Nd = 0


    Spe_Atom_Cut1 = zeros(Float64, Nspecies)
    Spe_Core_Charge = zeros(Float64, Nspecies)
    for spe = 1:Nspecies
        Spe_Atom_Cut1[spe] = pao[spe].Spe_Atom_Cut1
        Spe_Core_Charge[spe] = pspot[spe].Spe_Core_Charge
    end

    
    for spe = 1:Nspecies
        bc = Spe_Atom_Cut1[spe]
        dx = pi/sqrt(Scale_Grid_Ecut)
        Nd = 2*Int64(div(bc,dx)) + 1
        if Max_Nd<Nd
            Max_Nd = Nd
        end
    end


    Max_TGN_EH0 = 0
    for spe = 1:Nspecies
        bc = Spe_Atom_Cut1[spe]
        dx = pi/sqrt(Scale_Grid_Ecut)

        Nd = 2*floor(Int64, bc/dx) + 1
        dx = 2.0*bc/(Nd-1)
        gnum = Nd*CoarseGL_Mesh
        if Max_TGN_EH0 < gnum
            Max_TGN_EH0 = gnum
        end
    end
    Max_TGN_EH0 = Max_TGN_EH0 + 10


    g0 = zeros(Float64, Max_Nd)
    TGN_EH0 = zeros(Int64, Nspecies)
    dv_EH0 = zeros(Float64, Nspecies)
    CoarseGL_x, CoarseGL_Weight = Gauss_Legendre(CoarseGL_Mesh)

    GridX_EH0 = Vector{Vector{Float64}}(undef, Nspecies)
    GridY_EH0 = Vector{Vector{Float64}}(undef, Nspecies)
    GridZ_EH0 = Vector{Vector{Float64}}(undef, Nspecies)
    Arho_EH0 = Vector{Vector{Float64}}(undef, Nspecies)
    Wt_EH0 = Vector{Vector{Float64}}(undef, Nspecies)

    for spe = 1:Nspecies
        GridX_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        GridY_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        GridZ_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        Arho_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        Wt_EH0[spe] = zeros(Float64, Max_TGN_EH0)
    end
    

    for spe = 1:Nspecies

        Spe_Num_Mesh_PAO = pao[spe].Spe_Num_Mesh_PAO
        Spe_PAO_XV = pao[spe].Spe_PAO_XV
        Spe_PAO_RV = pao[spe].Spe_PAO_RV
        Spe_Atomic_Den = pao[spe].Spe_Atomic_Den

        bc = Spe_Atom_Cut1[spe]
        dx = pi/sqrt(Scale_Grid_Ecut)
        Nd = 2*Int64(div(bc,dx)) + 1
        dx = 2.0*bc/(Nd-1)
        dv_EH0[spe] = dx

        for n1 = 1:Nd
            g0[n1] = dx*(n1-1) - bc
        end

        gnum = 0
        y = 0.0
        Sx = Spe_Atom_Cut1[spe]
        Dx = Spe_Atom_Cut1[spe]

        for n3 = 1:Nd
            z = g0[n3]
            tmp = z^2
            for n1 = 1:CoarseGL_Mesh
                x = 0.5*(Dx*CoarseGL_x[n1] + Sx)
                xx = 0.5*log(x^2 + tmp)
                
                gnum += 1

                GridX_EH0[spe][gnum] = x
                GridY_EH0[spe][gnum] = y
                GridZ_EH0[spe][gnum] = z

                Arho_EH0[spe][gnum] = KumoF(Spe_Num_Mesh_PAO, xx, Spe_PAO_XV, Spe_PAO_RV, Spe_Atomic_Den)
                Wt_EH0[spe][gnum] = pi*x*CoarseGL_Weight[n1]*Dx
            end
        end

        TGN_EH0[spe] = gnum
    end



    VH_Atom = Vector{Vector{Float64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_Num_Mesh_VPS = pspot[spe].Spe_Num_Mesh_VPS

        VH_Atom[spe] = zeros(Float64, Spe_Num_Mesh_VPS+2)
        Calc_Spe_VH_Atom!(pao[spe], pspot[spe], VH_Atom[spe])
    end



    # calculation of scaling factors
    EH0_scaling = zeros(Float64, Nspecies, Nspecies)
    for spe = 1:Nspecies

        r1cut = Spe_Atom_Cut1[spe]
        Z1 = Spe_Core_Charge[spe]

        for jspe = 1:Nspecies
            r2cut = Spe_Atom_Cut1[jspe]
            Z2 = Spe_Core_Charge[jspe]

            Spe_Num_Mesh_VPS = pspot[jspe].Spe_Num_Mesh_VPS
            Spe_VPS_XV = pspot[jspe].Spe_VPS_XV
            Spe_VPS_RV = pspot[jspe].Spe_VPS_RV
            Spe_VH_Atom = VH_Atom[jspe]

            Sum = 0.0
            rcut = r1cut + r2cut
            for n1 = 1:TGN_EH0[spe]
                x = GridX_EH0[spe][n1]
                y = GridY_EH0[spe][n1]
                z = GridZ_EH0[spe][n1]
                rho0 = Arho_EH0[spe][n1]
                wt = Wt_EH0[spe][n1]
                z2 = z - rcut
                r2 = x^2 + y^2 + z2^2
                r = sqrt(r2)
                xx = 0.5*log(r2)
                r = ifelse(r<1e-10, 1e-10, r)

                va0 = VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
                Sum += wt*va0*rho0
            end

            EH0ij = Sum*dv_EH0[spe]

            EH0_scaling[spe,jspe] = ifelse(abs(EH0ij)>1e-20, Z1*Z2/rcut/EH0ij, 0.0)
        end
    end




    EH0 = 0.0
    EH0Force = zeros(Float64, Natom, 3)

    for loop = 1:MPI_size
        atom = MPI_atom[loop]
        spe = atom2spe[atom]
        Rn = MPI_FNAN[loop]
        cell = MPI_ncn[loop]+1
        jatom = MPI_natn[loop]
        jspe = atom2spe[jatom]
        Z2 = Spe_Core_Charge[jspe]

        Spe_Num_Mesh_VPS = pspot[jspe].Spe_Num_Mesh_VPS
        Spe_VPS_XV = pspot[jspe].Spe_VPS_XV
        Spe_VPS_RV = pspot[jspe].Spe_VPS_RV
        Spe_VH_Atom = VH_Atom[jspe]

        TmpEh0 = 0.0
        Sum = 0.0
        Sumr = 0.0

        for n1 = 1:TGN_EH0[spe]
            x = GridX_EH0[spe][n1]
            y = GridY_EH0[spe][n1]
            z = GridZ_EH0[spe][n1]
            rho0 = Arho_EH0[spe][n1]
            wt = Wt_EH0[spe][n1]
            z2 = z - Dis[atom][Rn]
            r2 = x^2 + y^2 + z2^2
            r = sqrt(r2)
            xx = 0.5*log(r2)
            r = ifelse(r<1e-10, 1e-10, r)

            va0 = VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
            Sum += wt*va0*rho0

            if Rn ≠ 1 && r > 1.0e-14
                dr_va0 = Dr_VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
                Sumr -= wt*dr_va0*rho0*z2/r
            end
        end
        factor = ifelse(Rn==1, 1.0, EH0_scaling[spe,jspe])
        TmpEh0 = TmpEh0 - factor*Sum*dv_EH0[spe]*0.25

        if Rn ≠ 1
            r = Dis[atom][Rn]
            r = ifelse(r<1e-10, 1e-10, r)

            x = Gxyz[atom][1] - Gxyz[jatom][1] - atv[cell][1]
            y = Gxyz[atom][2] - Gxyz[jatom][2] - atv[cell][2]
            z = Gxyz[atom][3] - Gxyz[jatom][3] - atv[cell][3]

            Sumr = Sumr*dv_EH0[spe]
            Sumx = Sumr*x/r
            Sumy = Sumr*y/r
            Sumz = Sumr*z/r
        else
            Sumx = 0.0
            Sumy = 0.0
            Sumz = 0.0
        end

        EH0Force[atom,1] = EH0Force[atom,1] - 0.5*factor*Sumx
        EH0Force[atom,2] = EH0Force[atom,2] - 0.5*factor*Sumy
        EH0Force[atom,3] = EH0Force[atom,3] - 0.5*factor*Sumz


        Rm = RMI[atom][Rn][1]+1
        katom = natn[jatom][Rm]
        cell2 = ncn[jatom][Rm]+1
        kspe = atom2spe[katom]
        Z2 = Spe_Core_Charge[kspe]

        Spe_Num_Mesh_VPS = pspot[kspe].Spe_Num_Mesh_VPS
        Spe_VPS_XV = pspot[kspe].Spe_VPS_XV
        Spe_VPS_RV = pspot[kspe].Spe_VPS_RV
        Spe_VH_Atom = VH_Atom[kspe]


        Sum = 0.0
        Sumr = 0.0

        for n1 = 1:TGN_EH0[jspe]
            x = GridX_EH0[jspe][n1]
            y = GridY_EH0[jspe][n1]
            z = GridZ_EH0[jspe][n1]
            rho0 = Arho_EH0[jspe][n1]
            wt = Wt_EH0[jspe][n1]
            z2 = z - Dis[jatom][Rm]
            r2 = x^2 + y^2 + z2^2
            r = sqrt(r2)
            xx = 0.5*log(r2)
            r = ifelse(r<1e-10, 1e-10, r)

            va0 = VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
            Sum += wt*va0*rho0

            if Rm ≠ 1 && r > 1.0e-14
                dr_va0 = Dr_VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
                Sumr -= wt*dr_va0*rho0*z2/r
            end
        end
        factor = ifelse(Rn==1, 1.0, EH0_scaling[jspe,spe])
        TmpEh0 = TmpEh0 - factor*Sum*dv_EH0[spe]*0.25
        EH0 += TmpEh0

        if Rm ≠ 1
            r = Dis[jatom][Rm]
            r = ifelse(r<1e-10, 1e-10, r)

            x = Gxyz[jatom][1] - Gxyz[katom][1] - atv[cell2][1]
            y = Gxyz[jatom][2] - Gxyz[katom][2] - atv[cell2][2]
            z = Gxyz[jatom][3] - Gxyz[katom][3] - atv[cell2][3]

            Sumr = Sumr*dv_EH0[jspe]
            Sumx = Sumr*x/r
            Sumy = Sumr*y/r
            Sumz = Sumr*z/r
        else
            Sumx = 0.0
            Sumy = 0.0
            Sumz = 0.0
        end

        EH0Force[atom,1] = EH0Force[atom,1] + 0.5*factor*Sumx
        EH0Force[atom,2] = EH0Force[atom,2] + 0.5*factor*Sumy
        EH0Force[atom,3] = EH0Force[atom,3] + 0.5*factor*Sumz
    end


    EH0 = MPI.Allreduce(EH0, MPI.SUM, comm)
    MPI.Allreduce!(EH0Force, MPI.SUM, comm)


    return EH0, EH0Force
end


function Calc_Ekin(SpinPol::String, DM, Hkin)

    Ekin = 0.0
    if SpinPol == "off"
        Ekin += dot(DM[1], Hkin)
    elseif SpinPol ∈ ("on", "nc")
        Ekin += dot(DM[1], Hkin)
        Ekin += dot(DM[2], Hkin)
    else
        error("please check SpinPol")
    end


    if SpinPol == "off"
        Ekin = 2*Ekin
    end

    return Ekin
end


function Calc_Ena(SpinPol::String, DM, HVNA)
    
    Ena = 0.0
    if SpinPol == "off"
        Ena += dot(DM[1], HVNA)
    elseif SpinPol ∈ ("on", "nc")
        Ena += dot(DM[1], HVNA)
        Ena += dot(DM[2], HVNA)
    else
        error("please check SpinPol")
    end

    if SpinPol == "off"
        Ena = 2*Ena
    end

    return Ena
end


function Calc_Enl(SpinPol::String, DM, iDM, HNL, iHNL)
       
    Enl = 0.0
    if SpinPol == "off"
        Enl += dot(DM[1], HNL[1])
    elseif SpinPol == "on"
        Enl += dot(DM[1], HNL[1])
        Enl += dot(DM[2], HNL[1])
    elseif SpinPol == "nc"
        Enl +=   dot( DM[1],  HNL[1])
        Enl -=   dot(iDM[1], iHNL[1])
        Enl +=   dot( DM[2],  HNL[2])
        Enl -=   dot(iDM[2], iHNL[2])
        Enl += 2*dot( DM[3],  HNL[3])
        Enl -= 2*dot( DM[4], iHNL[3])
    else
        error("please check SpinPol")
    end

    if SpinPol == "off"
        Enl = 2*Enl
    end
    

    return Enl
end


function Calc_EXC(
    SpinPol::AbstractString,
    pao::Vector{PAO}, pspot::Vector{Pspot}, 
    system_grid::System_Grid, 
    ADensity_Grid, PCCDensity_Grid, Density_Grid)

    Nspecies = length(pao)
    if Nspecies ≠ length(pspot)
        error("please check Nspecies")
    end

    xc_type = Vector{String}(undef, Nspecies)
    for spe = 1:Nspecies
        xc_type[spe] = pspot[spe].xc_type
    end
    if !allequal(xc_type)
        error("please check xc_type")
    end

    
    EXC = Calc_EXC1(SpinPol, xc_type[1], system_grid, ADensity_Grid, PCCDensity_Grid, Density_Grid)
    Exc2, ExcForce = Calc_EXC2(pao, pspot, system_grid)

    EXC[1] += 0.5*Exc2
    EXC[2] += 0.5*Exc2


    if SpinPol == "off"
        EXC[2] = EXC[1]
    end


    return EXC, ExcForce
end


function Calc_EH1(SpinPol::AbstractString, GridVol, ADensity_Grid, Density_Grid, dVHart_Grid)

    EH1 = 0.0
    if SpinPol == "off"
        EH1 = dot(2*(Density_Grid[1] .- ADensity_Grid), dVHart_Grid)
    elseif SpinPol ∈ ("on", "nc")
        EH1 = dot((Density_Grid[1] + Density_Grid[2] .- 2*ADensity_Grid), dVHart_Grid)
    end
    EH1 = EH1*GridVol*0.5

    return EH1
end


function Calc_EXC1(SpinPol::String, xc_type::String, system_grid::System_Grid, ADensity_Grid, PCCDensity_Grid, Density_Grid)

    Ngrid = system_grid.Ngrid
    Latvecs = system_grid.Latvecs
    GridVol = system_grid.GridVol
    NN = prod(Ngrid)
    gtv = zeros(Float64, 3, 3)
    gtv[1,:] = Latvecs[1,:]/Ngrid[1]
	gtv[2,:] = Latvecs[2,:]/Ngrid[2]
	gtv[3,:] = Latvecs[3,:]/Ngrid[3]

    if SpinPol == "off"
        spinmax = 1
    elseif SpinPol ∈ ("on", "nc")
        spinmax = 2
    end

    Vxc_Grid = Calc_Vxc_Grid( xc_type, SpinPol, gtv, Ngrid, PCCDensity_Grid, Density_Grid )


    EXC = zeros(Float64, 2)
    for spin = 1:spinmax
        Sum = 0.0
        for i = 1:NN
            Sum += (Density_Grid[spin][i] + PCCDensity_Grid[i])*Vxc_Grid[spin][i]
            Sum -= (ADensity_Grid[i] + PCCDensity_Grid[i])*LDA_CA(2*(ADensity_Grid[i]+PCCDensity_Grid[i]), 0) 
        end
        EXC[spin] = Sum
    end

    EXC[1] = EXC[1]*GridVol
    EXC[2] = EXC[2]*GridVol


    return EXC
end


function Calc_EXC2(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    Gxyz = system_grid.Gxyz
    atv = system_grid.atv
    natn = system_grid.natn
    ncn = system_grid.ncn
    FNAN = system_grid.FNAN
    Atom_Cut1 = system_grid.Atom_Cut1
    Nspecies = length(pao)


    # calculation of Exc^(0) and its contribution to forces on the fine mesh
    Spe_Atomic_Den2 = Vector{Vector{Float64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_Atomic_Den2[spe] = Calc_Atomic_Den2(pao[spe], pspot[spe])
    end

    CoarseGL_Abscissae, CoarseGL_Weight = Gauss_Legendre(CoarseGL_Mesh)
    Leb_Grid_XYZW = Set_Lebedev_Grid()

    myrange = split_evenly(1:CoarseGL_Mesh, nprocs)
    MPI_size = length(myrange[myrank+1])
    MPI_CoarseGL_Weight = CoarseGL_Weight[myrange[myrank+1]]
    MPI_CoarseGL_Abscissae = CoarseGL_Abscissae[myrange[myrank+1]]


    maxFNAN = maximum(FNAN)
    gx = zeros(Float64, maxFNAN+1)
    gy = zeros(Float64, maxFNAN+1)
    gz = zeros(Float64, maxFNAN+1)
    sum_gx = zeros(Float64, maxFNAN+1)
    sum_gy = zeros(Float64, maxFNAN+1)
    sum_gz = zeros(Float64, maxFNAN+1)
    sum_rx = zeros(Float64, maxFNAN+1)
    sum_ry = zeros(Float64, maxFNAN+1)
    sum_rz = zeros(Float64, maxFNAN+1)


    # calculation of Exc^(0) and its contribution to forces on the fine mesh
    Exc2 = 0.0
    den0 = 0.0


    Force = zeros(Float64, Natom, 3)
    for atom = 1:Natom

        fill!(sum_rx, 0.0)
        fill!(sum_ry, 0.0)
        fill!(sum_rz, 0.0)

        Rcut = Atom_Cut1[atom]
        Sumr = 0.0
        
        for loop = 1:MPI_size

            fill!(sum_gx, 0.0)
            fill!(sum_gy, 0.0)
            fill!(sum_gz, 0.0)

            r = 0.5*(Rcut*MPI_CoarseGL_Abscissae[loop] + Rcut)
            sumt = 0.0

            for ia = 1:Num_Leb_Grid

                x0 = r*Leb_Grid_XYZW[ia,1] + Gxyz[atom][1]
                y0 = r*Leb_Grid_XYZW[ia,2] + Gxyz[atom][2]
                z0 = r*Leb_Grid_XYZW[ia,3] + Gxyz[atom][3]
                den = 0.0

                for Rn = 1:FNAN[atom]+1

                    jatom = natn[atom][Rn]
                    cell = ncn[atom][Rn]+1
                    jspe = atom2spe[jatom]
                    r2cut = Atom_Cut1[jatom]^2
                    
                    Spe_Num_Mesh_PAO = pao[jspe].Spe_Num_Mesh_PAO
                    Spe_PAO_XV = pao[jspe].Spe_PAO_XV
                    Spe_PAO_RV = pao[jspe].Spe_PAO_RV

                    x1 = Gxyz[jatom][1] + atv[cell][1]
                    y1 = Gxyz[jatom][2] + atv[cell][2]
                    z1 = Gxyz[jatom][3] + atv[cell][3]

                    dx = x1 - x0
                    dy = y1 - y0
                    dz = z1 - z0
                    
                    r2 = dx*dx + dy*dy + dz*dz
                
                    gx[Rn] = 0.0
                    gy[Rn] = 0.0
                    gz[Rn] = 0.0

                    if r2 < r2cut
                        x = 0.5*log(r2)
                        den += KumoF(Spe_Num_Mesh_PAO, x, Spe_PAO_XV, Spe_PAO_RV, Spe_Atomic_Den2[jspe])

                        if Rn == 1
                            den0 = den
                        end

                        if Rn ≠ 1
                            r1 = sqrt(r2)
                            gden = Dr_KumoF(Spe_Num_Mesh_PAO, x, r1, Spe_PAO_XV, Spe_PAO_RV, Spe_Atomic_Den2[jspe])

                            gx[Rn] = gden/r1*dx
                            gy[Rn] = gden/r1*dy
                            gz[Rn] = gden/r1*dz
                        end
                    end
                end

                exc0 = LDA_CA(den, 0)
                dexc0 = LDA_CA(den, 3)

                weight = Leb_Grid_XYZW[ia,4]
                sumt += weight*den0*exc0

                for Rn = 2:FNAN[atom]+1
                    sum_gx[Rn] += weight*den0*dexc0*gx[Rn]
                    sum_gy[Rn] += weight*den0*dexc0*gy[Rn]
                    sum_gz[Rn] += weight*den0*dexc0*gz[Rn]
                end
            end

            weight = r^2 * MPI_CoarseGL_Weight[loop]
            Sumr += weight*sumt
            for Rn = 2:FNAN[atom]+1
                sum_rx[Rn] += weight*sum_gx[Rn]
                sum_ry[Rn] += weight*sum_gy[Rn]
                sum_rz[Rn] += weight*sum_gz[Rn]
            end
        end

        Sum = MPI.Allreduce(Sumr, MPI.SUM, comm)
        Exc2 += 2*pi*Rcut*Sum

        MPI.Allreduce!(sum_rx, MPI.SUM, comm)
        MPI.Allreduce!(sum_ry, MPI.SUM, comm)
        MPI.Allreduce!(sum_rz, MPI.SUM, comm)

        for Rn = 2:FNAN[atom]+1
            jatom = natn[atom][Rn]
            
            Force[jatom,1] += 2*pi*Rcut*sum_rx[Rn]
            Force[jatom,2] += 2*pi*Rcut*sum_ry[Rn]
            Force[jatom,3] += 2*pi*Rcut*sum_rz[Rn]

            Force[atom,1] -= 2*pi*Rcut*sum_rx[Rn]
            Force[atom,2] -= 2*pi*Rcut*sum_ry[Rn]
            Force[atom,3] -= 2*pi*Rcut*sum_rz[Rn]
        end
    end
    MPI.Barrier(comm)


    return Exc2, Force
end


function Calc_Atomic_Den2(pao::PAO, pspot::Pspot)

    Spe_Num_Mesh_PAO = pao.Spe_Num_Mesh_PAO
    Spe_PAO_XV = pao.Spe_PAO_XV
    Spe_Atomic_Den = pao.Spe_Atomic_Den

    Spe_Core_Charge = pspot.Spe_Core_Charge
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_XV = pspot.Spe_VPS_XV
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Atomic_PCC = pspot.Spe_Atomic_PCC


    # calculate Spe_Atomic_Den2
    Spe_Atomic_Den2 = deepcopy(Spe_Atomic_Den)

    Sum = 0.0
    dx = Spe_PAO_XV[2] - Spe_PAO_XV[1]
    for i = 1:Spe_Num_Mesh_PAO
        Sum += Spe_Atomic_Den[i+1]*exp(3.0*Spe_PAO_XV[i])
    end
    Sum *= 4*pi*dx
    

    for i = 1:Spe_Num_Mesh_PAO
        Spe_Atomic_Den2[i+1] = Spe_Atomic_Den2[i+1]*Spe_Core_Charge/Sum + KumoF(Spe_Num_Mesh_VPS, Spe_PAO_XV[i], Spe_VPS_XV, Spe_VPS_RV, Spe_Atomic_PCC)
    end
    Spe_Atomic_Den2[1] = 2*Spe_Atomic_Den2[2] - Spe_Atomic_Den2[3]
    Spe_Atomic_Den2[Spe_Num_Mesh_PAO+2] = 2*Spe_Atomic_Den2[Spe_Num_Mesh_PAO+1] - Spe_Atomic_Den2[Spe_Num_Mesh_PAO]


    return Spe_Atomic_Den2
end
