mutable struct Force
    Natom::Int32
    Gxyz::Vector{Vector{Float64}}
    Atoms_Symbol::Vector{String}
    PCCForce::Matrix{Float64}
    HkinForce::Matrix{Float64}
    OLPForce::Matrix{Float64}
    HVNAForce::Matrix{Float64}
    VpotForce::Matrix{Float64}
    HNLForce::Matrix{Float64}
    CoreForce::Matrix{Float64}
    EH0Force::Matrix{Float64}
    ExcForce::Matrix{Float64}
    ForceAll::Matrix{Float64}
end


function Init_Force(Natom, Gxyz, Atoms_Symbol)

    PCCForce = zeros(Float64, Natom, 3)
    HkinForce = zeros(Float64, Natom, 3)
    OLPForce = zeros(Float64, Natom, 3)
    HVNAForce = zeros(Float64, Natom, 3)
    VpotForce = zeros(Float64, Natom, 3)
    HNLForce = zeros(Float64, Natom, 3)
    CoreForce = zeros(Float64, Natom, 3)
    EH0Force = zeros(Float64, Natom, 3)
    ExcForce = zeros(Float64, Natom, 3)
    ForceAll = zeros(Float64, Natom, 3)

    return Force(Natom, Gxyz, Atoms_Symbol,
                 PCCForce, HkinForce, OLPForce, HVNAForce, 
                 VpotForce, HNLForce, CoreForce, EH0Force,
                 ExcForce, ForceAll)
end


function Print_Force(title::String, Natom, force::Matrix{Float64})
    for atom = 1:Natom
        @printf("  %s  atom = %d  %15.12f  %15.12f  %15.12f\n", title, atom, force[atom,1], force[atom,2], force[atom,3])
    end
end


@timeit timer "Force" function Force!(
    force::Force, electron::CrystalBloch,
    DM, iDM, Orbs_Grid,
    ADensity_Grid, PCCDensity_Grid, 
    dVHart_Grid, Vxc_Grid, Vpot_Grid,
    Ham::Hamiltonian, ucell::UCell, pao::Vector{PAO}, pspot::Vector{Pspot}; 
    EH0_flag::Bool=false, Exc_flag::Bool=false)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    system_grid = ucell.system_grid
    SpinPol = Ham.SpinPol
    HNL = Ham.HNL
    iHNL = Ham.iHNL
    HVNA = Ham.HVNA
    Core_Charge = electron.Core_Charge
    Natom = system_grid.Natom


    HNL_Vec = Set_HNL2HNL_Vec(HNL, system_grid)
    if !isnothing(iHNL)
        iHNL_Vec = Set_HNL2HNL_Vec(iHNL, system_grid)
    end
    HVNA_Vec = Set_HVNA2HVNA_Vec(HVNA, system_grid)



    myrank == 0 && println("<PCC_Force>")
    PCCForce = PCC_Force(SpinPol, ADensity_Grid, PCCDensity_Grid, Vxc_Grid, dVHart_Grid, pao, pspot, ucell)
    myrank == 0 && Print_Force("PCC_Force", Natom, PCCForce)

    
    myrank == 0 && println("<Kinetic_Force>")
    OLP_force, Hkin_force = Set_OLP_Kinforce(pao, system_grid)
    HkinForce = Kinetic_Force(SpinPol, Hkin_force, DM, system_grid)
    myrank == 0 && Print_Force("Kinetic_Force", Natom, HkinForce)


    myrank == 0 && println("<Force3>")
    if SpinPol == "off"
        VpotForce = Force3_nospin(pao, Orbs_Grid, Vpot_Grid, DM, ucell)
    elseif SpinPol == "on"
        VpotForce = Force3_spin(pao, Orbs_Grid, Vpot_Grid, DM, ucell)
    else
        VpotForce = Force3_nc(pao, Orbs_Grid, dVHart_Grid, Vxc_Grid, DM, ucell)
    end
    myrank == 0 && Print_Force("Force3", Natom, VpotForce)



    
    myrank == 0 && println("<HVNA_Force>")
    DS_VNAforce = Set_DS_VNAforce(pao, pspot, system_grid)
    HVNA2force, HVNA3force = Set_HVNA2_3force(pao, pspot, system_grid)    
    HVNAForce = HVNA_Force(SpinPol, DS_VNAforce, HVNA2force, HVNA3force, HVNA_Vec, DM, system_grid)
    myrank == 0 && Print_Force("HVNA_Force", Natom, HVNAForce)


    myrank == 0 && println("<OLP_Force>")
    OLPForce = OLP_Force(OLP_force, electron, system_grid, false)
    myrank == 0 && Print_Force("OLP_Force", Natom, OLPForce)

    
    myrank == 0 && println("<HNL_Force>")
    NLPforce = Set_NLPforce(pao, pspot, system_grid)
    if SpinPol ∈ ("off", "on")
        HNLForce = HNL_Force(SpinPol, pspot, NLPforce, HNL_Vec, DM, system_grid)
    else 
        HNLForce = HNL_Force_NC(pspot, NLPforce, HNL_Vec, iHNL_Vec, DM, iDM, system_grid)
    end
    myrank == 0 && Print_Force("HNL_Force", Natom, HNLForce)

    
    myrank == 0 && println("<Core_Force>")
    CoreForce = Core_Force(Core_Charge, system_grid)
    myrank == 0 && Print_Force("CoreForce", Natom, CoreForce)


    myrank == 0 && println("<EH0_Force>")
    if EH0_flag
        EH0Force = EH0_Force( pao, pspot, system_grid )
    else
        EH0Force = force.EH0Force
    end
    myrank == 0 && Print_Force("EH0Force", Natom, EH0Force)
    

    myrank == 0 && println("<Exc_Force>")
    if Exc_flag
        ExcForce = Exc_Force(pao, pspot, system_grid)
    else
        ExcForce = force.ExcForce
    end
    myrank == 0 && Print_Force("ExcForce", Natom, ExcForce)
    

    myrank == 0 && println("<ForceAll>")
    ForceAll = PCCForce + HkinForce + VpotForce + HVNAForce + OLPForce + HNLForce + CoreForce + EH0Force + ExcForce
    myrank == 0 && Print_Force("ForceAll", Natom, ForceAll)


    force.PCCForce = PCCForce
    force.HkinForce = HkinForce
    force.OLPForce = OLPForce
    force.HVNAForce = HVNAForce
    force.VpotForce = VpotForce
    force.HNLForce = HNLForce
    force.CoreForce = CoreForce
    force.EH0Force = EH0Force
    force.ExcForce = ExcForce
    force.ForceAll = ForceAll
end


@timeit timer "PCC_Force" function PCC_Force(SpinPol, ADensity_Grid, PCCDensity_Grid, Vxc_Grid, dVHart_Grid, pao::Vector{PAO}, pspot::Vector{Pspot}, ucell::UCell)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    
    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Latvecs = system_grid.Latvecs
    atom2spe = system_grid.atom2spe
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom
    atv = system_grid.atv
    Gxyz = system_grid.Gxyz
    Grid_Origin = system_grid.Grid_Origin
    GridVol = system_grid.GridVol
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid


    gLatvecs = zeros(Float64,3,3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
    gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
    gLatvecs[3,:] = Latvecs[3,:]/Ngrid3


    Nloop = Natom
    OneD2atom = zeros(Int64, Nloop)

    counts = 1
    for atom = 1:Natom
        OneD2atom[counts] = atom
        counts += 1
    end

    myrange = split_evenly(1:Nloop, nprocs)
    MPI_size = length(myrange[myrank+1])
    MPI_atom = OneD2atom[myrange[myrank+1]]



    RefVxc_Grid = zeros(Float64, prod(Ngrid))
    @. RefVxc_Grid = LDA_CA(2*(ADensity_Grid+PCCDensity_Grid), 1)



    PCCForce = zeros(Float64, Natom, 3)


    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        spe = atom2spe[atom]
        Spe_Num_Mesh_PAO = pao[spe].Spe_Num_Mesh_PAO
        Spe_PAO_XV = pao[spe].Spe_PAO_XV
        Spe_PAO_RV = pao[spe].Spe_PAO_RV
        Spe_Atomic_Den = pao[spe].Spe_Atomic_Den
        Spe_Num_Mesh_VPS = pspot[spe].Spe_Num_Mesh_VPS
        Spe_VPS_XV = pspot[spe].Spe_VPS_XV
        Spe_VPS_RV = pspot[spe].Spe_VPS_RV
        Spe_Atomic_PCC = pspot[spe].Spe_Atomic_PCC

        Sumx = 0.0
        Sumy = 0.0
        Sumz = 0.0

        for Nc = 1:GridN_Atom[atom]

            GNc = GridListAtom[atom][Nc]
            GRc = CellListAtom[atom][Nc]+1
            GN = GNc + 1

            n1 = div(GNc, Ngrid2*Ngrid3)
            n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
            n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3

            x = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + atv[GRc][1] + Grid_Origin[1]
            y = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + atv[GRc][2] + Grid_Origin[2]
            z = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + atv[GRc][3] + Grid_Origin[3]

            dx = Gxyz[atom][1] - x
            dy = Gxyz[atom][2] - y
            dz = Gxyz[atom][3] - z
            
            r = sqrt(dx^2 + dy^2 + dz^2)
            xx = log(r)

            r = ifelse(r < 1e-10, 1e-10, r)

            
            if r > 1e-14

                tmp0 = Dr_KumoF(Spe_Num_Mesh_PAO, xx, r, Spe_PAO_XV, Spe_PAO_RV, Spe_Atomic_Den)
                tmp1 = dVHart_Grid[GN]*tmp0/r
                Sumx += tmp1*dx
                Sumy += tmp1*dy
                Sumz += tmp1*dz
            
                tmp1 = RefVxc_Grid[GN]*tmp0/r
                Sumx += tmp1*dx
                Sumy += tmp1*dy
                Sumz += tmp1*dz

                tmp0 = 0.5*Dr_KumoF(Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_Atomic_PCC)

                if SpinPol == "off"
                    tmp2 = 2*Vxc_Grid[1][GN]
                else
                    tmp2 = Vxc_Grid[1][GN] + Vxc_Grid[2][GN]
                end
            
                tmp1 = tmp2*tmp0/r
                Sumx -= tmp1*dx
                Sumy -= tmp1*dy
                Sumz -= tmp1*dz

                tmp2 = 2*RefVxc_Grid[GN]
                tmp1 = tmp2*tmp0/r
                Sumx += tmp1*dx
                Sumy += tmp1*dy
                Sumz += tmp1*dz
            end
        end

        PCCForce[atom,1] = -Sumx * GridVol
        PCCForce[atom,2] = -Sumy * GridVol
        PCCForce[atom,3] = -Sumz * GridVol
    end
    MPI.Allreduce!(PCCForce, MPI.SUM, comm)


    return PCCForce
end


@timeit timer "Kinetic_Force" function Kinetic_Force(SpinPol, Hkin_force, DM, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    
    Natom = system_grid.Natom
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    RMI = system_grid.RMI
    Total_NumOrbs = system_grid.Total_NumOrbs

    maxTotal_NumOrbs = maximum(Total_NumOrbs)
    Hx = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hy = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hz = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)

    
    HkinForce = zeros(Float64, Natom, 3)
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]

        start_Rm = ifelse(SpinPol == "nc", 1, Rn)

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for Rm = start_Rm:FNAN[atom]+1
            kl = RMI[atom][Rn][Rm]
            katom = natn[atom][Rm]

            if kl >= 0

                fill!(Hx, 0.0)
                fill!(Hy, 0.0)
                fill!(Hz, 0.0)

                if Rn == 1
                    for ist = 1:Total_NumOrbs[jatom], jst = 1:Total_NumOrbs[katom]
                        Hx[ist,jst] += Hkin_force[1][atom][Rm][ist][jst]
                        Hy[ist,jst] += Hkin_force[2][atom][Rm][ist][jst]
                        Hz[ist,jst] += Hkin_force[3][atom][Rm][ist][jst]
                    end
                elseif Rn ≠ 1 && Rm == 1
                    for ist = 1:Total_NumOrbs[jatom], jst = 1:Total_NumOrbs[katom]
                        Hx[ist,jst] += Hkin_force[1][atom][Rn][jst][ist]
                        Hy[ist,jst] += Hkin_force[2][atom][Rn][jst][ist]
                        Hz[ist,jst] += Hkin_force[3][atom][Rn][jst][ist]
                    end
                end

                if SpinPol == "off"
                    pref = ifelse(isequal(Rn,Rm), 2.0, 4.0)
                    for ist = 1:Total_NumOrbs[jatom], jst = 1:Total_NumOrbs[katom]
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                elseif SpinPol == "on"
                    pref = ifelse(isequal(Rn,Rm), 1.0, 2.0)
                    for ist = 1:Total_NumOrbs[jatom], jst = 1:Total_NumOrbs[katom]
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                        dEx += pref*DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                else 
                    for ist = 1:Total_NumOrbs[jatom], jst = 1:Total_NumOrbs[katom]
                        dEx += DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEx += DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEy += DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                        dEz += DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                end
            end
        end

        HkinForce[atom,1] += dEx
        HkinForce[atom,2] += dEy
        HkinForce[atom,3] += dEz
    end
    MPI.Allreduce!(HkinForce, MPI.SUM, comm)


    return HkinForce
end


@timeit timer "Force3" function Force3_nospin(pao::Vector{PAO}, Orbs_Grid, Vpot_Grid, DM, ucell::UCell)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_size = system_grid.MPI_size
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    GridVol = system_grid.GridVol
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2


    dOrbs_Grid = Set_dOrbitals_Grid(pao, ucell)
    VpotForce = zeros(Float64, Natom, 3)

    for loop = 1:MPI_size
        
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        sumx1, sumy1, sumz1 = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                             dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vpot_Grid[1], DM[1][atom][Rn])

        VpotForce[atom,1] += 4*sumx1*GridVol
        VpotForce[atom,2] += 4*sumy1*GridVol
        VpotForce[atom,3] += 4*sumz1*GridVol
    end

    MPI.Allreduce!(VpotForce, MPI.SUM, comm)


    return VpotForce
end


@timeit timer "Force3" function Force3_spin(pao::Vector{PAO}, Orbs_Grid, Vpot_Grid, DM, ucell::UCell)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_size = system_grid.MPI_size
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn

    GridVol = system_grid.GridVol
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2

    dOrbs_Grid = Set_dOrbitals_Grid(pao, ucell)


    VpotForce = zeros(Float64, Natom, 3)

    for loop = 1:MPI_size
        
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]

        sumx1_up, sumy1_up, sumz1_up = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                      dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vpot_Grid[1], DM[1][atom][Rn])

        sumx1_dn, sumy1_dn, sumz1_dn = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                      dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vpot_Grid[2], DM[2][atom][Rn])

        VpotForce[atom,1] += 2*(sumx1_up + sumx1_dn)*GridVol
        VpotForce[atom,2] += 2*(sumy1_up + sumy1_dn)*GridVol
        VpotForce[atom,3] += 2*(sumz1_up + sumz1_dn)*GridVol
    end

    MPI.Allreduce!(VpotForce, MPI.SUM, comm)


    return VpotForce
end


@timeit timer "Force3" function Force3_nc(pao::Vector{PAO}, Orbs_Grid, dVHart_Grid, Vxc_Grid, DM, ucell::UCell)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    system_grid = ucell.system_grid
    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MPI_size = system_grid.MPI_size
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn

    GridVol = system_grid.GridVol
    GridListAtom = ucell.GridListAtom
    MPI_NumOLG = ucell.MPI_NumOLG
    MPI_GListTAtoms1 = ucell.MPI_GListTAtoms1
    MPI_GListTAtoms2 = ucell.MPI_GListTAtoms2


    Vpot1 = dVHart_Grid + Vxc_Grid[1]
    Vpot2 = dVHart_Grid + Vxc_Grid[2]
    dOrbs_Grid = Set_dOrbitals_Grid(pao, ucell)


    VpotForce = zeros(Float64, Natom, 3)

    for loop = 1:MPI_size
        
        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        
        sumx11, sumy11, sumz11 = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vpot1, DM[1][atom][Rn])

        sumx12, sumy12, sumz12 = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vpot2, DM[2][atom][Rn])

        sumx13, sumy13, sumz13 = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vxc_Grid[3], DM[3][atom][Rn])

        sumx14, sumy14, sumz14 = _Calc_Force3_8(NO0, NO1, MPI_NumOLG[loop], GridListAtom[atom], MPI_GListTAtoms1[loop], MPI_GListTAtoms2[loop], 
                                                dOrbs_Grid[1][atom], dOrbs_Grid[2][atom], dOrbs_Grid[3][atom], Orbs_Grid[jatom], Vxc_Grid[4], DM[4][atom][Rn])

        VpotForce[atom,1] += (2*sumx11 + 2*sumx12 + 4*sumx13 - 4*sumx14)*GridVol
        VpotForce[atom,2] += (2*sumy11 + 2*sumy12 + 4*sumy13 - 4*sumy14)*GridVol
        VpotForce[atom,3] += (2*sumz11 + 2*sumz12 + 4*sumz13 - 4*sumz14)*GridVol
    end

    MPI.Allreduce!(VpotForce, MPI.SUM, comm)


    return VpotForce
end


@timeit timer "HVNA_Force" function HVNA_Force(SpinPol, DS_VNA, HVNA2, HVNA3, HVNA, DM, system_grid::System_Grid)

    comm = MPI.COMM_WORLD

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    RMI = system_grid.RMI
    Dis = system_grid.Dis
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    Total_NumOrbs = system_grid.Total_NumOrbs
    maxTotal_NumOrbs = maximum(Total_NumOrbs)
    VNATotal_Num = length(DS_VNA[1][1][1][1])

    Hx = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hy = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hz = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)

    
    HVNAForce = zeros(Float64, Natom, 3)
    
    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = 1
        Rm = MPI_FNAN[loop]
        jatom = natn[atom][Rn]
        katom = natn[atom][Rm]
        kl = RMI[atom][Rn][Rm]

        NO0 = Total_NumOrbs[jatom]
        NO1 = Total_NumOrbs[katom]

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        fill!(Hx, 0.0)
        fill!(Hy, 0.0)
        fill!(Hz, 0.0)

        dHVNA!( 0, atom, Rn, Rm, 
                Hx, Hy, Hz, 
                HVNA2[1][atom], HVNA2[2][atom], HVNA2[3][atom], 
                HVNA3[1][atom], HVNA3[2][atom], HVNA3[3][atom], 
                DS_VNA, HVNA[atom], system_grid)

        if SpinPol == "off"
            pref = ifelse(isequal(Rn,Rm), 2.0, 4.0)
            for ist = 1:NO0, jst = 1:NO1
                dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
            end
        else
            pref = ifelse(isequal(Rn,Rm), 1.0, 2.0)
            for ist = 1:NO0, jst = 1:NO1
                dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEx += pref*DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEy += pref*DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                dEz += pref*DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
            end
        end

        HVNAForce[atom,1] += dEx
        HVNAForce[atom,2] += dEy
        HVNAForce[atom,3] += dEz
    end
    

    OneDatom, _, OneDFNAN, HVNA_FNAN, MP_FNAN = Set_H_FNAN(Natom, FNAN, natn, ncn, atv_ijk)
    
    loop = 1
    atom_old = MPI_atom[loop]
    Set_DS_VNA_Natom!(atom_old, OneDatom, OneDFNAN, HVNA_FNAN, MP_FNAN, FNAN, Total_NumOrbs, VNATotal_Num, DS_VNA)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[jatom]

        if atom ≠ atom_old
            Set_DS_VNA_Natom!(atom, OneDatom, OneDFNAN, HVNA_FNAN, MP_FNAN, FNAN, Total_NumOrbs, VNATotal_Num, DS_VNA)
        end

        atom_old = atom
        
        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for Rm = Rn:FNAN[atom]+1
        
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            NO1 = Total_NumOrbs[kg]

            if kl >= 0

                fill!(Hx, 0.0)
                fill!(Hy, 0.0)
                fill!(Hz, 0.0)

                dHVNA!( 1, atom, Rn, Rm,
                        Hx, Hy, Hz, 
                        HVNA2[1][atom], HVNA2[2][atom], HVNA2[3][atom], 
                        HVNA3[1][atom], HVNA3[2][atom], HVNA3[3][atom], 
                        DS_VNA, HVNA[atom], system_grid)

                if SpinPol == "off"
                    pref = ifelse(isequal(Rn,Rm), 2.0, 4.0)
                    for ist = 1:NO0, jst = 1:NO1
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                else
                    pref = ifelse(isequal(Rn,Rm), 1.0, 2.0)
                    for ist = 1:NO0, jst = 1:NO1
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEx += pref*DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEy += pref*DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                        dEz += pref*DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                end
            end
        end

        HVNAForce[atom,1] += dEx
        HVNAForce[atom,2] += dEy
        HVNAForce[atom,3] += dEz
    end
    MPI.Allreduce!(HVNAForce, MPI.SUM, comm)
    

    return HVNAForce
end


@timeit timer "OLP_Force" function OLP_Force(OLP_force, electron::CrystalBloch, system_grid::System_Grid, occ_flag::Bool)

    SpinPol = electron.SpinPol
    EDM = Calc_EDM(electron, system_grid, occ_flag)
    OLPForce = OLP_Force(SpinPol, EDM, OLP_force, system_grid)

    return OLPForce
end


function OLP_Force(SpinPol, EDM, OLP_force, system_grid::System_Grid)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    
    OLPForce = zeros(Float64, Natom, 3)

    if SpinPol == "off"
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            NO0 = Total_NumOrbs[atom]
            jatom = natn[atom][Rn]
            NO1 = Total_NumOrbs[jatom]
            for ist = 1:NO0, jst = 1:NO1
                hst += 1
                if Rn ≠ 1
                    dum = 2*EDM[1][hst]

                    dx = dum*OLP_force[1][hst]
                    dy = dum*OLP_force[2][hst]
                    dz = dum*OLP_force[3][hst]

                    OLPForce[atom,1] = OLPForce[atom,1] - 2*dx
                    OLPForce[atom,2] = OLPForce[atom,2] - 2*dy
                    OLPForce[atom,3] = OLPForce[atom,3] - 2*dz
                end
            end
        end
    else
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            NO0 = Total_NumOrbs[atom]
            jatom = natn[atom][Rn]
            NO1 = Total_NumOrbs[jatom]
            for ist = 1:NO0, jst = 1:NO1
                hst += 1
                if Rn ≠ 1
                    dum = EDM[1][hst] + EDM[2][hst]

                    dx = dum*OLP_force[1][hst]
                    dy = dum*OLP_force[2][hst]
                    dz = dum*OLP_force[3][hst]

                    OLPForce[atom,1] = OLPForce[atom,1] - 2*dx
                    OLPForce[atom,2] = OLPForce[atom,2] - 2*dy
                    OLPForce[atom,3] = OLPForce[atom,3] - 2*dz
                end
            end
        end
    end


    return OLPForce
end


@timeit timer "HNL_Force" function HNL_Force(SpinPol, pspot::Vector{Pspot}, NLP, HNL, DM, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD

    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    RMI = system_grid.RMI
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    Total_NumOrbs = system_grid.Total_NumOrbs
    maxTotal_NumOrbs = maximum(Total_NumOrbs)



    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot[spe].Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end

        NLTotal_Num[atom] = tot
    end

    VNLE = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        VNLE[atom] = zeros(Float64, NLTotal_Num[atom])

        count = 0
        for lst = 1:pspot[spe].Spe_Num_RVPS
            ene = pspot[spe].Spe_VNLE[1,lst]
            L = pspot[spe].Spe_VPS_List[lst]
            for _ = -L:L
                count += 1
                VNLE[atom][count] = ene
            end
        end
    end

    tmpL = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        tmpL[atom] = zeros(Float64, NLTotal_Num[atom])
    end


    Hx = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hy = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)
    Hz = zeros(Float64, maxTotal_NumOrbs, maxTotal_NumOrbs)

    HNLForce = zeros(Float64, Natom, 3)


    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = 1
        Rm = MPI_FNAN[loop]
        jatom = natn[atom][Rn]
        kl = RMI[atom][Rn][Rm]
        katom = natn[atom][Rm]

        NO0 = Total_NumOrbs[jatom]
        NO1 = Total_NumOrbs[katom]

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        fill!(Hx, 0.0)
        fill!(Hy, 0.0)
        fill!(Hz, 0.0)

        dHNL!( 0, atom, Rn, Rm, tmpL, NLTotal_Num, VNLE,
               Hx, Hy, Hz, 
               NLP, HNL[1][atom],
               system_grid)

        if SpinPol == "off"
            pref = ifelse(isequal(Rn,Rm), 2.0, 4.0)
            for ist = 1:NO0, jst = 1:NO1
                dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
            end
        elseif SpinPol == "on"
            pref = ifelse(isequal(Rn,Rm), 1.0, 2.0)
            for ist = 1:NO0, jst = 1:NO1
                dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                dEx += pref*DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                dEy += pref*DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                dEz += pref*DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
            end
        end

        HNLForce[atom,1] += dEx
        HNLForce[atom,2] += dEy
        HNLForce[atom,3] += dEz
    end



    OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN = Set_H_FNAN(Natom, FNAN, natn, ncn, atv_ijk)


    loop = 1
    atom_old = MPI_atom[loop]
    Set_NLP_Natom!(atom_old, OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN, FNAN, Total_NumOrbs, NLTotal_Num, 0, NLP)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[jatom]

        if atom ≠ atom_old
            Set_NLP_Natom!(atom, OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN, FNAN, Total_NumOrbs, NLTotal_Num, 0, NLP)
        end

        atom_old = atom
        
        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for Rm = Rn:FNAN[atom]+1
        
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            NO1 = Total_NumOrbs[kg]

            if kl >= 0
                fill!(Hx, 0.0)
                fill!(Hy, 0.0)
                fill!(Hz, 0.0)

                dHNL!( 1, atom, Rn, Rm, tmpL, NLTotal_Num, VNLE,
                       Hx, Hy, Hz, 
                       NLP, HNL[1][atom],
                       system_grid)

                if SpinPol == "off"
                    pref = ifelse(isequal(Rn,Rm), 2.0, 4.0)
                    for ist = 1:NO0, jst = 1:NO1
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                elseif SpinPol == "on"
                    pref = ifelse(isequal(Rn,Rm), 1.0, 2.0)
                    for ist = 1:NO0, jst = 1:NO1
                        dEx += pref*DM[1][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[1][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[1][jatom][kl+1][ist][jst]*Hz[ist,jst]
                        dEx += pref*DM[2][jatom][kl+1][ist][jst]*Hx[ist,jst]
                        dEy += pref*DM[2][jatom][kl+1][ist][jst]*Hy[ist,jst]
                        dEz += pref*DM[2][jatom][kl+1][ist][jst]*Hz[ist,jst]
                    end
                end
            end
        end

        HNLForce[atom,1] += dEx
        HNLForce[atom,2] += dEy
        HNLForce[atom,3] += dEz
    end
    MPI.Allreduce!(HNLForce, MPI.SUM, comm)


    return HNLForce
end


# Calc HNL_Force for NonCollinear case
@timeit timer "HNL_Force" function HNL_Force_NC(pspot::Vector{Pspot}, NLP, HNL, iHNL, DM, iDM, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    
    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    RMI = system_grid.RMI
    MPI_atom = system_grid.MPI_atom
    MPI_FNAN = system_grid.MPI_FNAN
    MPI_natn = system_grid.MPI_natn
    MPI_size = system_grid.MPI_size
    Total_NumOrbs = system_grid.Total_NumOrbs
    maxTotal_NumOrbs = maximum(Total_NumOrbs)


    Atoms_Num_RVPS = zeros(Int64, Natom)
    NLTotal_Num = zeros(Int64, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        tot = 0
        List = pspot[spe].Spe_VPS_List
        for list in List
            tot += 2*list + 1
        end
        Atoms_Num_RVPS[atom] = pspot[spe].Spe_Num_RVPS
        NLTotal_Num[atom] = tot
    end

    Atoms_VPS_List = Vector{Vector{Int64}}(undef, Natom)
    Atoms_VNLE = Vector{Array{Float64,2}}(undef, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        Atoms_VPS_List[atom] = zeros(Int64, Atoms_Num_RVPS[atom])
        Atoms_VNLE[atom] = zeros(Float64, 2, Atoms_Num_RVPS[atom])

        for lst = 1:Atoms_Num_RVPS[atom]
            Atoms_VNLE[atom][1,lst] = pspot[spe].Spe_VNLE[1,lst]
            Atoms_VNLE[atom][2,lst] = pspot[spe].Spe_VNLE[2,lst]
            Atoms_VPS_List[atom][lst] = pspot[spe].Spe_VPS_List[lst]
        end
    end



    Hx = zeros(ComplexF64, maxTotal_NumOrbs, maxTotal_NumOrbs, 3)
    Hy = zeros(ComplexF64, maxTotal_NumOrbs, maxTotal_NumOrbs, 3)
    Hz = zeros(ComplexF64, maxTotal_NumOrbs, maxTotal_NumOrbs, 3)


    HNLForce = zeros(Float64, Natom, 3)


    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = 1
        Rm = MPI_FNAN[loop]
        jatom = natn[atom][Rn]
        kl = RMI[atom][Rn][Rm]
        katom = natn[atom][Rm]

        NO0 = Total_NumOrbs[jatom]
        NO1 = Total_NumOrbs[katom]

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        fill!(Hx, 0.0)
        fill!(Hy, 0.0)
        fill!(Hz, 0.0)

        dHNL_NC!( 0, atom, Rn, Rm,
                   Atoms_Num_RVPS, Atoms_VNLE, Atoms_VPS_List, 
                   Hx, Hy, Hz, 
                   NLP, HNL, iHNL,
                   system_grid)

        if Rn == Rm
            for ist = 1:NO0, jst = 1:NO1
                dEx +=   DM[1][jatom][kl+1][ist][jst]*real(Hx[ist,jst,1])
                dEx -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,1])
                dEx +=   DM[2][jatom][kl+1][ist][jst]*real(Hx[ist,jst,2])
                dEx -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,2])
                dEx += 2*DM[3][jatom][kl+1][ist][jst]*real(Hx[ist,jst,3])
                dEx -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,3])
        
                dEy +=   DM[1][jatom][kl+1][ist][jst]*real(Hy[ist,jst,1])
                dEy -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,1])
                dEy +=   DM[2][jatom][kl+1][ist][jst]*real(Hy[ist,jst,2])
                dEy -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,2])
                dEy += 2*DM[3][jatom][kl+1][ist][jst]*real(Hy[ist,jst,3])
                dEy -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,3])
        
                dEz +=   DM[1][jatom][kl+1][ist][jst]*real(Hz[ist,jst,1])
                dEz -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,1])
                dEz +=   DM[2][jatom][kl+1][ist][jst]*real(Hz[ist,jst,2])
                dEz -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,2])
                dEz += 2*DM[3][jatom][kl+1][ist][jst]*real(Hz[ist,jst,3])
                dEz -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,3])
            end
        else
            for ist = 1:NO0, jst = 1:NO1
                dEx +=   DM[1][jatom][kl+1][ist][jst]*real(Hx[ist,jst,1])
                dEx -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,1])
                dEx +=   DM[2][jatom][kl+1][ist][jst]*real(Hx[ist,jst,2])
                dEx -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,2])
                dEx += 2*DM[3][jatom][kl+1][ist][jst]*real(Hx[ist,jst,3])
                dEx -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,3])
        
                dEy +=   DM[1][jatom][kl+1][ist][jst]*real(Hy[ist,jst,1])
                dEy -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,1])
                dEy +=   DM[2][jatom][kl+1][ist][jst]*real(Hy[ist,jst,2])
                dEy -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,2])
                dEy += 2*DM[3][jatom][kl+1][ist][jst]*real(Hy[ist,jst,3])
                dEy -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,3])
        
                dEz +=   DM[1][jatom][kl+1][ist][jst]*real(Hz[ist,jst,1])
                dEz -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,1])
                dEz +=   DM[2][jatom][kl+1][ist][jst]*real(Hz[ist,jst,2])
                dEz -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,2])
                dEz += 2*DM[3][jatom][kl+1][ist][jst]*real(Hz[ist,jst,3])
                dEz -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,3])
            end
    
            fill!(Hx, 0.0)
            fill!(Hy, 0.0)
            fill!(Hz, 0.0)
    
            dHNL_NC!( 0, atom, Rm, Rn,
                Atoms_Num_RVPS, Atoms_VNLE, Atoms_VPS_List, 
                Hx, Hy, Hz, 
                NLP, HNL, iHNL,
                system_grid)
    
            kl1 = RMI[atom][Rm][Rn]
            for ist = 1:NO1, jst = 1:NO0
                dEx +=   DM[1][katom][kl1+1][ist][jst]*real(Hx[ist,jst,1])
                dEx -=  iDM[1][katom][kl1+1][ist][jst]*imag(Hx[ist,jst,1])
                dEx +=   DM[2][katom][kl1+1][ist][jst]*real(Hx[ist,jst,2])
                dEx -=  iDM[2][katom][kl1+1][ist][jst]*imag(Hx[ist,jst,2])
                dEx += 2*DM[3][katom][kl1+1][ist][jst]*real(Hx[ist,jst,3])
                dEx -= 2*DM[4][katom][kl1+1][ist][jst]*imag(Hx[ist,jst,3])
    
                dEy +=   DM[1][katom][kl1+1][ist][jst]*real(Hy[ist,jst,1])
                dEy -=  iDM[1][katom][kl1+1][ist][jst]*imag(Hy[ist,jst,1])
                dEy +=   DM[2][katom][kl1+1][ist][jst]*real(Hy[ist,jst,2])
                dEy -=  iDM[2][katom][kl1+1][ist][jst]*imag(Hy[ist,jst,2])
                dEy += 2*DM[3][katom][kl1+1][ist][jst]*real(Hy[ist,jst,3])
                dEy -= 2*DM[4][katom][kl1+1][ist][jst]*imag(Hy[ist,jst,3])
    
                dEz +=   DM[1][katom][kl1+1][ist][jst]*real(Hz[ist,jst,1])
                dEz -=  iDM[1][katom][kl1+1][ist][jst]*imag(Hz[ist,jst,1])
                dEz +=   DM[2][katom][kl1+1][ist][jst]*real(Hz[ist,jst,2])
                dEz -=  iDM[2][katom][kl1+1][ist][jst]*imag(Hz[ist,jst,2])
                dEz += 2*DM[3][katom][kl1+1][ist][jst]*real(Hz[ist,jst,3])
                dEz -= 2*DM[4][katom][kl1+1][ist][jst]*imag(Hz[ist,jst,3])
            end
        end

        HNLForce[atom,1] += dEx
        HNLForce[atom,2] += dEy
        HNLForce[atom,3] += dEz
    end



    OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN = Set_H_FNAN(Natom, FNAN, natn, ncn, atv_ijk)

    loop = 1
    atom_old = MPI_atom[loop]
    Set_NLP_Natom!(atom_old, OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN, FNAN, Total_NumOrbs, NLTotal_Num, 1, NLP)

    for loop = 1:MPI_size

        atom = MPI_atom[loop]
        Rn = MPI_FNAN[loop]
        jatom = MPI_natn[loop]
        NO0 = Total_NumOrbs[jatom]

        if atom ≠ atom_old
            Set_NLP_Natom!(atom, OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN, FNAN, Total_NumOrbs, NLTotal_Num, 1, NLP)
        end

        atom_old = atom
        
        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for Rm = 1:FNAN[atom]+1
            kg = natn[atom][Rm]
            kl = RMI[atom][Rn][Rm]
            NO1 = Total_NumOrbs[kg]
            if kl >= 0
                fill!(Hx, 0.0)
                fill!(Hy, 0.0)
                fill!(Hz, 0.0)

                dHNL_NC!( 1, atom, Rn, Rm,
                        Atoms_Num_RVPS, Atoms_VNLE, Atoms_VPS_List, 
                        Hx, Hy, Hz, 
                        NLP, HNL, iHNL,
                        system_grid)

                for ist = 1:NO0, jst = 1:NO1
                    dEx +=   DM[1][jatom][kl+1][ist][jst]*real(Hx[ist,jst,1])
                    dEx -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,1])
                    dEx +=   DM[2][jatom][kl+1][ist][jst]*real(Hx[ist,jst,2])
                    dEx -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,2])
                    dEx += 2*DM[3][jatom][kl+1][ist][jst]*real(Hx[ist,jst,3])
                    dEx -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hx[ist,jst,3])

                    dEy +=   DM[1][jatom][kl+1][ist][jst]*real(Hy[ist,jst,1])
                    dEy -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,1])
                    dEy +=   DM[2][jatom][kl+1][ist][jst]*real(Hy[ist,jst,2])
                    dEy -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,2])
                    dEy += 2*DM[3][jatom][kl+1][ist][jst]*real(Hy[ist,jst,3])
                    dEy -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hy[ist,jst,3])

                    dEz +=   DM[1][jatom][kl+1][ist][jst]*real(Hz[ist,jst,1])
                    dEz -=  iDM[1][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,1])
                    dEz +=   DM[2][jatom][kl+1][ist][jst]*real(Hz[ist,jst,2])
                    dEz -=  iDM[2][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,2])
                    dEz += 2*DM[3][jatom][kl+1][ist][jst]*real(Hz[ist,jst,3])
                    dEz -= 2*DM[4][jatom][kl+1][ist][jst]*imag(Hz[ist,jst,3])
                end
            end
        end

        HNLForce[atom,1] += dEx
        HNLForce[atom,2] += dEy
        HNLForce[atom,3] += dEz
    end
    MPI.Allreduce!(HNLForce, MPI.SUM, comm)


    return HNLForce
end


@timeit timer "Core_Force" function Core_Force(Core_Charge::Vector{Float64}, system_grid::System_Grid)

    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    atv = system_grid.atv
    natn = system_grid.natn
    ncn = system_grid.ncn
    FNAN = system_grid.FNAN
    Dis = system_grid.Dis
    

    CoreForce = zeros(Float64, Natom, 3)
    for atom = 1:Natom
        Zc = Core_Charge[atom]

        dEx = 0.0
        dEy = 0.0
        dEz = 0.0

        for Rn = 2:FNAN[atom]+1
            jatom = natn[atom][Rn]
            cell = ncn[atom][Rn]+1
            Zh = Core_Charge[jatom]

            r = ifelse(Dis[atom][Rn]<1e-10, 1e-10, Dis[atom][Rn]) 

            lx = (Gxyz[atom][1] - Gxyz[jatom][1] - atv[cell][1])/r
            ly = (Gxyz[atom][2] - Gxyz[jatom][2] - atv[cell][2])/r
            lz = (Gxyz[atom][3] - Gxyz[jatom][3] - atv[cell][3])/r

            dEx = dEx - lx * Zc*Zh/r^2
            dEy = dEy - ly * Zc*Zh/r^2
            dEz = dEz - lz * Zc*Zh/r^2
        end

        CoreForce[atom,1] = dEx
        CoreForce[atom,2] = dEy
        CoreForce[atom,3] = dEz
    end

    return CoreForce
end


@timeit timer "EH0_Force" function EH0_Force(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    
    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    natn = system_grid.natn
    ncn = system_grid.ncn
    FNAN = system_grid.FNAN
    Gxyz = system_grid.Gxyz
    atv = system_grid.atv
    RMI = system_grid.RMI
    Dis = system_grid.Dis
    
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
        Nd = 2*floor(Int64, bc/dx) + 1
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

        Spe_Num_Mesh_PAO = pao[spe].Spe_Num_Mesh_PAO
        Spe_PAO_XV = pao[spe].Spe_PAO_XV
        Spe_PAO_RV = pao[spe].Spe_PAO_RV
        Spe_Atomic_Den = pao[spe].Spe_Atomic_Den


        bc = Spe_Atom_Cut1[spe]
        dx = pi/sqrt(Scale_Grid_Ecut)
        Nd = 2*floor(Int64, bc/dx) + 1
        dx = 2.0*bc/(Nd-1)
        dv_EH0[spe] = dx

        GridX_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        GridY_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        GridZ_EH0[spe] = zeros(Float64, Max_TGN_EH0)
    
        Arho_EH0[spe] = zeros(Float64, Max_TGN_EH0)
        Wt_EH0[spe] = zeros(Float64, Max_TGN_EH0)

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



    EH0Force = zeros(Float64, Natom, 3)


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

            sum = 0.0
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
                sum += wt*va0*rho0
            end

            EH0ij = sum*dv_EH0[spe]

            EH0_scaling[spe,jspe] = ifelse(abs(EH0ij)>1e-20, Z1*Z2/rcut/EH0ij, 0.0)
        end
    end



    # -1/2\int n^a(r) V^a_H dr
    for atom = 1:Natom
        spe = atom2spe[atom]

        for Rn = 1:FNAN[atom]+1

            jatom = natn[atom][Rn]
            cell = ncn[atom][Rn]+1
            jspe = atom2spe[jatom]
            Z2 = Spe_Core_Charge[jspe]

            Spe_Num_Mesh_VPS = pspot[jspe].Spe_Num_Mesh_VPS
            Spe_VPS_XV = pspot[jspe].Spe_VPS_XV
            Spe_VPS_RV = pspot[jspe].Spe_VPS_RV
            Spe_VH_Atom = VH_Atom[jspe]

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

                if Rn ≠ 1 && r > 1e-14
                    dr_va0 = Dr_VH_AtomF( Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom )
                    Sumr += wt*dr_va0*rho0*z2/r
                end
            end

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

            factor = ifelse(Rn==1, 1.0, EH0_scaling[spe,jspe])


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

                if Rn ≠ 0 && r > 1e-14
                    dr_va0 = Dr_VH_AtomF(Z2, Spe_Num_Mesh_VPS, xx, r, Spe_VPS_XV, Spe_VPS_RV, Spe_VH_Atom)
                    Sumr += wt*dr_va0*rho0*z2/r
                end
            end

            if Rm ≠ 1
                r = Dis[jatom][Rm]
                r = ifelse(r<1e-10, 1e-10, r)

                x = Gxyz[jatom][1] - Gxyz[katom][1] - atv[cell2][1]
                y = Gxyz[jatom][2] - Gxyz[katom][2] - atv[cell2][2]
                z = Gxyz[jatom][3] - Gxyz[katom][3] - atv[cell2][3]

                Sumr = Sumr*dv_EH0[kspe]
                Sumx = Sumr*x/r
                Sumy = Sumr*y/r
                Sumz = Sumr*z/r
            else
                Sumx = 0.0
                Sumy = 0.0
                Sumz = 0.0
            end

            factor = ifelse(Rn==1, 1.0, EH0_scaling[jspe,spe])
            
            EH0Force[atom,1] = EH0Force[atom,1] + 0.5*factor*Sumx
            EH0Force[atom,2] = EH0Force[atom,2] + 0.5*factor*Sumy
            EH0Force[atom,3] = EH0Force[atom,3] + 0.5*factor*Sumz
        end
    end
    MPI.Allreduce!(EH0Force, MPI.SUM, comm)
   

    return EH0Force
end


@timeit timer "Exc_Force" function Exc_Force(pao::Vector{PAO}, pspot::Vector{Pspot}, system_grid::System_Grid)

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
    den0 = 0.0


    Force = zeros(Float64, Natom, 3)
    for atom = 1:Natom

        fill!(sum_rx, 0.0)
        fill!(sum_ry, 0.0)
        fill!(sum_rz, 0.0)

        Rcut = Atom_Cut1[atom]
        
        for loop = 1:MPI_size

            fill!(sum_gx, 0.0)
            fill!(sum_gy, 0.0)
            fill!(sum_gz, 0.0)

            r = 0.5*(Rcut*MPI_CoarseGL_Abscissae[loop] + Rcut)

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

                dexc0 = LDA_CA(den, 3)

                weight = Leb_Grid_XYZW[ia,4]
                for Rn = 2:FNAN[atom]+1
                    sum_gx[Rn] += weight*den0*dexc0*gx[Rn]
                    sum_gy[Rn] += weight*den0*dexc0*gy[Rn]
                    sum_gz[Rn] += weight*den0*dexc0*gz[Rn]
                end
            end

            weight = r^2 * MPI_CoarseGL_Weight[loop]
            for Rn = 2:FNAN[atom]+1
                sum_rx[Rn] += weight*sum_gx[Rn]
                sum_ry[Rn] += weight*sum_gy[Rn]
                sum_rz[Rn] += weight*sum_gz[Rn]
            end
        end

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


    return Force
end


function Set_DM_Vec2DM(DM, system_grid::System_Grid)

    Nspin = length(DM)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    Total_Hsize = system_grid.Total_Hsize

    DM_1D = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
		DM_1D[spin] = zeros(Float64, Total_Hsize)
	end
    
    for spin = 1:Nspin
        counts = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            counts += 1
            DM_1D[spin][counts] = DM[spin][atom][Rn][ist][jst]
        end
	end

    return DM_1D
end


function Set_DM2DM_Vec(DM, system_grid::System_Grid)

    Nspin = length(DM)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    DM_Vec = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspin)
	for spin = 1:Nspin
		DM_Vec[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
		for atom = 1:Natom
			DM_Vec[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
			for Rn = 1:FNAN[atom]+1
				DM_Vec[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
				for ist = 1:Total_NumOrbs[atom]
					DM_Vec[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
				end
			end
		end
	end

    for spin = 1:Nspin
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            hst += 1
            DM_Vec[spin][atom][Rn][ist][jst] = DM[spin][hst]
        end
	end

    return DM_Vec
end


function Set_HNL2HNL_Vec(HNL, system_grid::System_Grid)

    Nspin = length(HNL)
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HNL_Vec = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspin)
    for spin = 1:Nspin
        HNL_Vec[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
        for atom = 1:Natom
            HNL_Vec[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                HNL_Vec[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    HNL_Vec[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end

    
    for spin = 1:Nspin
        counts = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
            counts += 1
            HNL_Vec[spin][atom][Rn][ist][jst] = HNL[spin][counts]
        end
    end


    return HNL_Vec
end


function Set_HVNA2HVNA_Vec(HVNA, system_grid::System_Grid)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    HVNA_Vec = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
    for atom = 1:Natom
        HVNA_Vec[atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            HVNA_Vec[atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                HVNA_Vec[atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
            end
        end
    end

    counts = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        counts += 1
        HVNA_Vec[atom][Rn][ist][jst] = HVNA[counts]
    end


    return HVNA_Vec
end


function Set_H_FNAN(Natom, FNAN, natn, ncn, atv_ijk)

    Total_FNAN = sum(FNAN.+1)
    H_FNAN = zeros(Int32, Total_FNAN)
    OneDatom = zeros(Int32, Total_FNAN)
    OneDjatom = zeros(Int32, Total_FNAN)
    OneDFNAN = zeros(Int32, Total_FNAN)
    num = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        num += 1
        OneDatom[num] = atom
        OneDjatom[num] = natn[atom][Rn]
        OneDFNAN[num] = Rn
    end

    OneDjatom_sort = sortperm(OneDjatom)
    OneDatom = OneDatom[OneDjatom_sort]
    OneDjatom = OneDjatom[OneDjatom_sort]
    OneDFNAN = OneDFNAN[OneDjatom_sort]

    for num = 1:Total_FNAN

        atom = OneDatom[num]
        jatom = OneDjatom[num]
        Rn = OneDFNAN[num]
        cell = ncn[atom][Rn]+1
        m1, m2, m3 = -atv_ijk[cell]

        po = 0
        j = 0
        while po == 0

            katom = natn[jatom][j+1]
            cell2 = ncn[jatom][j+1]+1

            n1, n2, n3 = atv_ijk[cell2]

            if m1==n1 && m2==n2 && m3==n3 && atom == katom
                H_FNAN[num] = j
                po = 1
            end

            j += 1

            if FNAN[jatom] < j && po == 0
                po = 2
            end
        end
    end


    Sum = 0
    MP_FNAN = zeros(Int32, Natom+1)
    for atom = 1:Natom
        MP_FNAN[atom+1] = Sum + FNAN[atom] + 1
        Sum += FNAN[atom]+1
    end

    return OneDatom, OneDjatom, OneDFNAN, H_FNAN, MP_FNAN
end


function Set_DS_VNA_Natom!(atom, OneDatom, OneDFNAN, HVNA_FNAN, MP_FNAN, FNAN, Total_NumOrbs, VNATotal_Num, DS_VNA)
    
    num = MP_FNAN[atom]
    for Rn = 1:FNAN[atom]+1
        atom1 = OneDatom[Rn+num]
        Rm = OneDFNAN[Rn+num]
        Rl = HVNA_FNAN[Rn+num]+1
        NO0 = Total_NumOrbs[atom1]
        for ist = 1:NO0, jst = 1:VNATotal_Num
            DS_VNA[1][end][Rl][ist][jst] = DS_VNA[1][atom1][Rm][ist][jst]
            DS_VNA[2][end][Rl][ist][jst] = DS_VNA[2][atom1][Rm][ist][jst]
            DS_VNA[3][end][Rl][ist][jst] = DS_VNA[3][atom1][Rm][ist][jst]
            DS_VNA[4][end][Rl][ist][jst] = DS_VNA[4][atom1][Rm][ist][jst]
        end
    end
end


function Set_NLP_Natom!(atom, OneDatom, OneDjatom, OneDFNAN, HNL_FNAN, MP_FNAN, FNAN, Total_NumOrbs, NLTotal_Num, VPS_j_dependency, DS_NL)
    
    num = MP_FNAN[atom]
    for Rn = 1:FNAN[atom]+1
        atom1 = OneDatom[Rn+num]
        jatom = OneDjatom[Rn+num]
        Rm = OneDFNAN[Rn+num]
        Rl = HNL_FNAN[Rn+num]+1
        NO0 = Total_NumOrbs[atom1]
        NO1 = NLTotal_Num[jatom]
        for ist = 1:NO0, so = 1:VPS_j_dependency+1, jst = 1:NO1
            DS_NL[1][end][Rl][ist][so][jst] = DS_NL[1][atom1][Rm][ist][so][jst]
            DS_NL[2][end][Rl][ist][so][jst] = DS_NL[2][atom1][Rm][ist][so][jst]
            DS_NL[3][end][Rl][ist][so][jst] = DS_NL[3][atom1][Rm][ist][so][jst]
            DS_NL[4][end][Rl][ist][so][jst] = DS_NL[4][atom1][Rm][ist][so][jst]
        end
    end
end


function dHVNA!(
    where_flag,
    atom, Rn, Rm,
    Hx, Hy, Hz, 
    HVNA2x, HVNA2y, HVNA2z, 
    HVNA3x, HVNA3y, HVNA3z, 
    DS_VNA, HVNA,
    system_grid::System_Grid)

    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv = system_grid.atv
    Total_NumOrbs = system_grid.Total_NumOrbs
    Atom_Cut1 = system_grid.Atom_Cut1
    RMI = system_grid.RMI
    Dis = system_grid.Dis

    ig = natn[atom][Rn]
    Rni = ncn[atom][Rn]+1
    NO0 = Total_NumOrbs[ig]
    Rcut1 = Atom_Cut1[ig]

    jg = natn[atom][Rm]
    Rnj = ncn[atom][Rm]+1
    NO1 = Total_NumOrbs[jg]
    Rcut2 = Atom_Cut1[jg]

    Rcut = Rcut1 + Rcut2
    kl = RMI[atom][Rn][Rm]
    dmp = dampingF(Rcut, Dis[ig][kl+1])
    


    if Rn == 1 && Rm == 1 && where_flag == 0
        for k = 2:FNAN[atom]+1, ist = 1:NO0, jst = 1:NO1
            Hx[ist,jst] += HVNA2x[k][ist][jst]
            Hy[ist,jst] += HVNA2y[k][ist][jst]
            Hz[ist,jst] += HVNA2z[k][ist][jst]
        end
    elseif Rn == Rm && Rn ≠ 1
        for ist = 1:NO0, jst = 1:NO1
            Hx[ist,jst] = -HVNA3x[Rn][ist][jst]
            Hy[ist,jst] = -HVNA3y[Rn][ist][jst]
            Hz[ist,jst] = -HVNA3z[Rn][ist][jst]
        end
    else
        if Rn == 1
            for Rl = 1:FNAN[atom]+1

                kl = RMI[atom][Rm][Rl]
                if kl >= 0 && where_flag == 0

                    jatom = ifelse(jg <= Natom, jg, Natom)
                    for ist = 1:NO0, jst = 1:NO1
                        Sumx = dot(DS_VNA[2][atom][Rl][ist], DS_VNA[1][jatom][kl+1][jst])
                        Sumy = dot(DS_VNA[3][atom][Rl][ist], DS_VNA[1][jatom][kl+1][jst])
                        Sumz = dot(DS_VNA[4][atom][Rl][ist], DS_VNA[1][jatom][kl+1][jst])

                        Hx[ist,jst] += Sumx
                        Hy[ist,jst] += Sumy
                        Hz[ist,jst] += Sumz
                    end
                end
            end

            if Rm == 1
                for ist = 1:NO0, jst = 1:NO1
                    tmpx = Hx[ist,jst] + Hx[jst,ist]
                    Hx[ist,jst] = tmpx
                    Hx[jst,ist] = tmpx

                    tmpy = Hy[ist,jst] + Hy[jst,ist]
                    Hy[ist,jst] = tmpy
                    Hy[jst,ist] = tmpy

                    tmpz = Hz[ist,jst] + Hz[jst,ist]
                    Hz[ist,jst] = tmpz
                    Hz[jst,ist] = tmpz
                end
            elseif where_flag == 1          
                
                if jg <= Natom
                    jatom = jg
                    kl = RMI[atom][Rm][1]
                else
                    jatom = Natom
                    kl = RMI[atom][1][Rm]
                end

                for ist = 1:NO0, jst = 1:NO1
                    Sumx = -dot(DS_VNA[2][jatom][kl+1][jst], DS_VNA[1][atom][1][ist])
                    Sumy = -dot(DS_VNA[3][jatom][kl+1][jst], DS_VNA[1][atom][1][ist])
                    Sumz = -dot(DS_VNA[4][jatom][kl+1][jst], DS_VNA[1][atom][1][ist])

                    Hx[ist,jst] += Sumx
                    Hy[ist,jst] += Sumy
                    Hz[ist,jst] += Sumz
                end
            end
        else
            kl1 = RMI[atom][1][Rn]+1
            kl2 = RMI[atom][1][Rm]+1

            for ist = 1:NO0, jst = 1:NO1
                Sumx = -dot(DS_VNA[2][end][kl1][ist], DS_VNA[1][end][kl2][jst])
                Sumy = -dot(DS_VNA[3][end][kl1][ist], DS_VNA[1][end][kl2][jst])
                Sumz = -dot(DS_VNA[4][end][kl1][ist], DS_VNA[1][end][kl2][jst])
                
                Hx[ist,jst] = Sumx
                Hy[ist,jst] = Sumy
                Hz[ist,jst] = Sumz
            end

            if Rm ≠ 1
                for ist = 1:NO0, jst = 1:NO1
                    Sumx = -dot(DS_VNA[2][end][kl2][jst], DS_VNA[1][end][kl1][ist])
                    Sumy = -dot(DS_VNA[3][end][kl2][jst], DS_VNA[1][end][kl1][ist])
                    Sumz = -dot(DS_VNA[4][end][kl2][jst], DS_VNA[1][end][kl1][ist])

                    Hx[ist,jst] += Sumx
                    Hy[ist,jst] += Sumy
                    Hz[ist,jst] += Sumz
                end
            end
        end
    end


    for ist = 1:NO0, jst = 1:NO1
        Hx[ist,jst] = dmp*Hx[ist,jst]
        Hy[ist,jst] = dmp*Hy[ist,jst]
        Hz[ist,jst] = dmp*Hz[ist,jst]
    end

    if (Rn == 1 && Rm ≠ 1) || (Rn ≠ 1 && Rm == 1)

        kl = ifelse(Rn == 1, Rm, Rn)
        r = Dis[atom][kl]

        if r >= Rcut
            deri_dmp = 0.0
            tmp = 0.0
        else
            deri_dmp = deri_dampingF(Rcut, r)
            tmp = deri_dmp/dmp
        end

        x0 = Gxyz[ig][1] + atv[Rni][1]
        x1 = Gxyz[jg][1] + atv[Rnj][1]
        y0 = Gxyz[ig][2] + atv[Rni][2]
        y1 = Gxyz[jg][2] + atv[Rnj][2]
        z0 = Gxyz[ig][3] + atv[Rni][3]
        z1 = Gxyz[jg][3] + atv[Rnj][3]

        r = ifelse(r<1.0e-10, 1.0e-10, r)

        if Rn==1
            dx = tmp*(x0-x1)/r
            dy = tmp*(y0-y1)/r
            dz = tmp*(z0-z1)/r
        elseif Rm == 1
            dx = tmp*(x1-x0)/r
            dy = tmp*(y1-y0)/r
            dz = tmp*(z1-z0)/r
        end
    
        if Rn == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst] += HVNA[kl][ist][jst]*dx
                Hy[ist,jst] += HVNA[kl][ist][jst]*dy
                Hz[ist,jst] += HVNA[kl][ist][jst]*dz
            end
        elseif Rm == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst] += HVNA[kl][jst][ist]*dx
                Hy[ist,jst] += HVNA[kl][jst][ist]*dy
                Hz[ist,jst] += HVNA[kl][jst][ist]*dz
            end
        end
    end
end


function dHNL!(
    where_flag,
    atom, Rn, Rm, 
    tmpL, NLTotal_Num, VNLE, 
    Hx, Hy, Hz, 
    NLP, HNL,
    system_grid::System_Grid)


    Gxyz = system_grid.Gxyz
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv = system_grid.atv
    RMI = system_grid.RMI
    Dis = system_grid.Dis
    Total_NumOrbs = system_grid.Total_NumOrbs
    Atom_Cut1 = system_grid.Atom_Cut1


    ig = natn[atom][Rn]
    Rni = ncn[atom][Rn]+1
    NO0 = Total_NumOrbs[ig]
    Rcut1 = Atom_Cut1[ig]


    jg = natn[atom][Rm]
    Rnj = ncn[atom][Rm]+1
    NO1 = Total_NumOrbs[jg]
    Rcut2 = Atom_Cut1[jg]

    Rcut = Rcut1 + Rcut2
    kl = RMI[atom][Rn][Rm]
    dmp = dampingF(Rcut, Dis[ig][kl+1])



    if Rn == 1
        jatom = jg
        for Rl = 1:FNAN[atom]+1
            kg = natn[atom][Rl]
            kl = RMI[atom][Rm][Rl]

            if kl >= 0 && where_flag == 0
                for ist = 1:NO0, jst = 1:NO1

                    # @. tmpL[kg] = VNLE[kg]*NLP[1][jatom][kl+1][jst][1]
                    # Sumx = dot(tmpL[kg], NLP[2][atom][Rl][ist][1])
                    # Sumy = dot(tmpL[kg], NLP[3][atom][Rl][ist][1])
                    # Sumz = dot(tmpL[kg], NLP[4][atom][Rl][ist][1])

                    
                    Sumx = 0.0
                    Sumy = 0.0
                    Sumz = 0.0
                    for l = 1:NLTotal_Num[kg]
                        Sumx += VNLE[kg][l]*NLP[1][jatom][kl+1][jst][1][l]*NLP[2][atom][Rl][ist][1][l]
                        Sumy += VNLE[kg][l]*NLP[1][jatom][kl+1][jst][1][l]*NLP[3][atom][Rl][ist][1][l]
                        Sumz += VNLE[kg][l]*NLP[1][jatom][kl+1][jst][1][l]*NLP[4][atom][Rl][ist][1][l]
                    end

                    Hx[ist,jst] += Sumx
                    Hy[ist,jst] += Sumy
                    Hz[ist,jst] += Sumz
                end
            end
        end


        if Rm == 1 
            for ist = 1:NO0, jst = ist:NO1
                tmpx = Hx[ist,jst] + Hx[jst,ist]
                Hx[ist,jst] = tmpx
                Hx[jst,ist] = tmpx

                tmpy = Hy[ist,jst] + Hy[jst,ist]
                Hy[ist,jst] = tmpy
                Hy[jst,ist] = tmpy

                tmpz = Hz[ist,jst] + Hz[jst,ist]
                Hz[ist,jst] = tmpz
                Hz[jst,ist] = tmpz
            end
        elseif where_flag == 1

            jatom = jg
            kg = natn[atom][1]
            kl = RMI[atom][Rm][1]+1

            for ist = 1:NO0, jst = 1:NO1

                #=
                @. tmpL[kg] = VNLE[kg]*NLP[1][atom][1][ist][1]
                Sumx = -dot(tmpL[kg], NLP[2][jatom][kl][jst][1])
                Sumy = -dot(tmpL[kg], NLP[3][jatom][kl][jst][1])
                Sumz = -dot(tmpL[kg], NLP[4][jatom][kl][jst][1])
                =#
                
                Sumx = 0.0
                Sumy = 0.0
                Sumz = 0.0
                for l = 1:NLTotal_Num[kg]
                    Sumx -= VNLE[kg][l]*NLP[1][atom][1][ist][1][l]*NLP[2][jatom][kl][jst][1][l]
                    Sumy -= VNLE[kg][l]*NLP[1][atom][1][ist][1][l]*NLP[3][jatom][kl][jst][1][l]
                    Sumz -= VNLE[kg][l]*NLP[1][atom][1][ist][1][l]*NLP[4][jatom][kl][jst][1][l]
                end

                Hx[ist,jst] += Sumx
                Hy[ist,jst] += Sumy
                Hz[ist,jst] += Sumz
            end
        end
    elseif where_flag == 0

    else
        kg = natn[atom][1]
        kl1 = RMI[atom][1][Rn]+1
        kl2 = RMI[atom][1][Rm]+1

        for ist = 1:NO0, jst = 1:NO1            
            Sumx = 0.0
            Sumy = 0.0
            Sumz = 0.0
            for l = 1:NLTotal_Num[kg]
                Sumx -= VNLE[kg][l]*NLP[1][end][kl2][jst][1][l]*NLP[2][end][kl1][ist][1][l]
                Sumy -= VNLE[kg][l]*NLP[1][end][kl2][jst][1][l]*NLP[3][end][kl1][ist][1][l]
                Sumz -= VNLE[kg][l]*NLP[1][end][kl2][jst][1][l]*NLP[4][end][kl1][ist][1][l]
            end

            Hx[ist,jst] = Sumx
            Hy[ist,jst] = Sumy
            Hz[ist,jst] = Sumz
        end

        if Rm ≠ 1
            for ist = 1:NO0, jst = 1:NO1
                Sumx = 0.0
                Sumy = 0.0
                Sumz = 0.0
                for l = 1:NLTotal_Num[kg]
                    Sumx -= VNLE[kg][l]*NLP[1][end][kl1][ist][1][l]*NLP[2][end][kl2][jst][1][l]
                    Sumy -= VNLE[kg][l]*NLP[1][end][kl1][ist][1][l]*NLP[3][end][kl2][jst][1][l]
                    Sumz -= VNLE[kg][l]*NLP[1][end][kl1][ist][1][l]*NLP[4][end][kl2][jst][1][l]
                end

                Hx[ist,jst] += Sumx
                Hy[ist,jst] += Sumy
                Hz[ist,jst] += Sumz
            end
        end
    end




    for ist = 1:NO0, jst = 1:NO1
        Hx[ist,jst] = dmp*Hx[ist,jst]
        Hy[ist,jst] = dmp*Hy[ist,jst]
        Hz[ist,jst] = dmp*Hz[ist,jst]
    end


    if (Rn == 1 && Rm ≠ 1) || (Rn ≠ 1 && Rm == 1)

        kl = ifelse(Rn == 1, Rm, Rn)

        r = Dis[atom][kl]

        if r >= Rcut
            deri_dmp = 0.0
            tmp = 0.0
        else
            deri_dmp = deri_dampingF(Rcut, r)
            tmp = deri_dmp/dmp
        end

        x0 = Gxyz[ig][1] + atv[Rni][1]
        x1 = Gxyz[jg][1] + atv[Rnj][1]

        y0 = Gxyz[ig][2] + atv[Rni][2]
        y1 = Gxyz[jg][2] + atv[Rnj][2]

        z0 = Gxyz[ig][3] + atv[Rni][3]
        z1 = Gxyz[jg][3] + atv[Rnj][3]

        r = ifelse(r<1.0e-10, 1.0e-10, r)

        if Rn==1 && Rm≠1
            dx = tmp*(x0-x1)/r
            dy = tmp*(y0-y1)/r
            dz = tmp*(z0-z1)/r
        elseif Rn≠1 &&  Rm==1
            dx = tmp*(x1-x0)/r
            dy = tmp*(y1-y0)/r
            dz = tmp*(z1-z0)/r
        end
    

        if Rn == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst] += HNL[kl][ist][jst]*dx
                Hy[ist,jst] += HNL[kl][ist][jst]*dy
                Hz[ist,jst] += HNL[kl][ist][jst]*dz
            end
        elseif Rm == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst] += HNL[kl][jst][ist]*dx
                Hy[ist,jst] += HNL[kl][jst][ist]*dy
                Hz[ist,jst] += HNL[kl][jst][ist]*dz
            end
        end
    end    
end


function dHNL_NC!(
    where_flag,
    atom, Rn, Rm, 
    Atoms_Num_RVPS, Atoms_VNLE, Atoms_VPS_List,
    Hx, Hy, Hz, 
    NLP, HNL, iHNL,
    system_grid::System_Grid)

    Natom = system_grid.Natom
    Gxyz = system_grid.Gxyz
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv = system_grid.atv
    RMI = system_grid.RMI
    Dis = system_grid.Dis
    Total_NumOrbs = system_grid.Total_NumOrbs
    Atom_Cut1 = system_grid.Atom_Cut1


    ig = natn[atom][Rn]
    Rni = ncn[atom][Rn]+1
    NO0 = Total_NumOrbs[ig]
    Rcut1 = Atom_Cut1[ig]

    jg = natn[atom][Rm]
    Rnj = ncn[atom][Rm]+1
    NO1 = Total_NumOrbs[jg]
    Rcut2 = Atom_Cut1[jg]

    Rcut = Rcut1 + Rcut2
    kl = RMI[atom][Rn][Rm]
    dmp = dampingF(Rcut, Dis[ig][kl+1])


    if Rn == 1
        jatom = jg
        for Rl = 1:FNAN[atom]+1
            kg = natn[atom][Rl]
            kl = RMI[atom][Rm][Rl]

            if kl >= 0 && where_flag == 0 
                dHNL_SO!(0, Hx, Hy, Hz, NO0, NO1, 1.0, Atoms_Num_RVPS[kg], Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][jatom][kl+1], NLP[2][atom][Rl], NLP[3][atom][Rl], NLP[4][atom][Rl])

                if Rm == 1
                    dHNL_SO!(1, Hx, Hy, Hz, NO1, NO0, -1.0, Atoms_Num_RVPS[kg], Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][atom][Rl], NLP[2][jatom][kl+1], NLP[3][jatom][kl+1], NLP[4][jatom][kl+1])
                end
            end
        end


        kg = natn[atom][FNAN[atom]+1]
        if where_flag == 1
            jatom = jg
            kg = natn[atom][1]
            kl = RMI[atom][Rm][1]+1
            dHNL_SO!(1, Hx, Hy, Hz, NO1, NO0, -1.0, Atoms_Num_RVPS[kg], -Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][atom][1], NLP[2][jatom][kl], NLP[3][jatom][kl], NLP[4][jatom][kl])
        end
    elseif where_flag == 0
        for k = 1:FNAN[atom]+1
            kg = natn[atom][k]
            kl = RMI[atom][Rn][k]
            if kl >= 0
                dHNL_SO!(1, Hx, Hy, Hz, NO1, NO0, -1.0, Atoms_Num_RVPS[kg], Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][ig][kl+1], NLP[2][jg][k], NLP[3][jg][k], NLP[4][jg][k])
            end
        end
    else

        kg = natn[atom][1]
        kl1 = RMI[atom][1][Rn]+1
        kl2 = RMI[atom][1][Rm]+1
        fill!(Hx, 0.0)
        fill!(Hy, 0.0)
        fill!(Hz, 0.0)

        dHNL_SO!(0, Hx, Hy, Hz, NO0, NO1, 1.0, Atoms_Num_RVPS[kg], -Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][end][kl2], NLP[2][end][kl1], NLP[3][end][kl1], NLP[4][end][kl1])
  
        if Rm ≠ 1
            dHNL_SO!(1, Hx, Hy, Hz, NO1, NO0, -1.0, Atoms_Num_RVPS[kg], -Atoms_VNLE[kg], Atoms_VPS_List[kg], NLP[1][end][kl1], NLP[2][end][kl2], NLP[3][end][kl2], NLP[4][end][kl2])
        end
    end




    for ist = 1:NO0, jst = 1:NO1
        Hx[ist,jst,1] = dmp*Hx[ist,jst,1]
        Hx[ist,jst,2] = dmp*Hx[ist,jst,2]
        Hx[ist,jst,3] = dmp*Hx[ist,jst,3]
        Hy[ist,jst,1] = dmp*Hy[ist,jst,1]
        Hy[ist,jst,2] = dmp*Hy[ist,jst,2]
        Hy[ist,jst,3] = dmp*Hy[ist,jst,3]
        Hz[ist,jst,1] = dmp*Hz[ist,jst,1]
        Hz[ist,jst,2] = dmp*Hz[ist,jst,2]
        Hz[ist,jst,3] = dmp*Hz[ist,jst,3]
    end


    if (Rn == 1 && Rm ≠ 1) || (Rn ≠ 1 && Rm == 1)

        kl = ifelse(Rn == 1, Rm, Rn)
        r = Dis[atom][kl]

        if r >= Rcut
            deri_dmp = 0.0
            tmp = 0.0
        else
            deri_dmp = deri_dampingF(Rcut, r)
            tmp = deri_dmp/dmp
        end

        x0 = Gxyz[ig][1] + atv[Rni][1]
        x1 = Gxyz[jg][1] + atv[Rnj][1]

        y0 = Gxyz[ig][2] + atv[Rni][2]
        y1 = Gxyz[jg][2] + atv[Rnj][2]

        z0 = Gxyz[ig][3] + atv[Rni][3]
        z1 = Gxyz[jg][3] + atv[Rnj][3]

        r = ifelse(r<1.0e-10, 1.0e-10, r)

        if Rn==1 && Rm≠1
            dx = tmp*(x0-x1)/r
            dy = tmp*(y0-y1)/r
            dz = tmp*(z0-z1)/r
        elseif Rn≠1 && Rm==1
            dx = tmp*(x1-x0)/r
            dy = tmp*(y1-y0)/r
            dz = tmp*(z1-z0)/r
        end
    

        if Rn == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst,1] += HNL[1][atom][kl][ist][jst]*dx
                Hx[ist,jst,2] += HNL[2][atom][kl][ist][jst]*dx
                Hx[ist,jst,3] += HNL[3][atom][kl][ist][jst]*dx
                Hy[ist,jst,1] += HNL[1][atom][kl][ist][jst]*dy
                Hy[ist,jst,2] += HNL[2][atom][kl][ist][jst]*dy
                Hy[ist,jst,3] += HNL[3][atom][kl][ist][jst]*dy
                Hz[ist,jst,1] += HNL[1][atom][kl][ist][jst]*dz
                Hz[ist,jst,2] += HNL[2][atom][kl][ist][jst]*dz
                Hz[ist,jst,3] += HNL[3][atom][kl][ist][jst]*dz
            end
        elseif Rm == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst,1] += HNL[1][atom][kl][jst][ist]*dx
                Hx[ist,jst,2] += HNL[2][atom][kl][jst][ist]*dx
                Hx[ist,jst,3] += HNL[3][atom][kl][jst][ist]*dx
                Hy[ist,jst,1] += HNL[1][atom][kl][jst][ist]*dy
                Hy[ist,jst,2] += HNL[2][atom][kl][jst][ist]*dy
                Hy[ist,jst,3] += HNL[3][atom][kl][jst][ist]*dy
                Hz[ist,jst,1] += HNL[1][atom][kl][jst][ist]*dz
                Hz[ist,jst,2] += HNL[2][atom][kl][jst][ist]*dz
                Hz[ist,jst,3] += HNL[3][atom][kl][jst][ist]*dz
            end
        end


        if Rn == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst,1] += im*iHNL[1][atom][kl][ist][jst]*dx
                Hx[ist,jst,2] += im*iHNL[2][atom][kl][ist][jst]*dx
                Hx[ist,jst,3] += im*iHNL[3][atom][kl][ist][jst]*dx
                Hy[ist,jst,1] += im*iHNL[1][atom][kl][ist][jst]*dy
                Hy[ist,jst,2] += im*iHNL[2][atom][kl][ist][jst]*dy
                Hy[ist,jst,3] += im*iHNL[3][atom][kl][ist][jst]*dy
                Hz[ist,jst,1] += im*iHNL[1][atom][kl][ist][jst]*dz
                Hz[ist,jst,2] += im*iHNL[2][atom][kl][ist][jst]*dz
                Hz[ist,jst,3] += im*iHNL[3][atom][kl][ist][jst]*dz
            end
        elseif Rm == 1
            for ist = 1:NO0, jst = 1:NO1
                Hx[ist,jst,1] += im*iHNL[1][atom][kl][jst][ist]*dx
                Hx[ist,jst,2] += im*iHNL[2][atom][kl][jst][ist]*dx
                Hx[ist,jst,3] += im*iHNL[3][atom][kl][jst][ist]*dx
                Hy[ist,jst,1] += im*iHNL[1][atom][kl][jst][ist]*dy
                Hy[ist,jst,2] += im*iHNL[2][atom][kl][jst][ist]*dy
                Hy[ist,jst,3] += im*iHNL[3][atom][kl][jst][ist]*dy
                Hz[ist,jst,1] += im*iHNL[1][atom][kl][jst][ist]*dz
                Hz[ist,jst,2] += im*iHNL[2][atom][kl][jst][ist]*dz
                Hz[ist,jst,3] += im*iHNL[3][atom][kl][jst][ist]*dz
            end
        end
    end    
end


function dHNL_SO!(
    orb_switch,
    Hx, Hy, Hz, NO0, NO1, fugou, 
    Spe_Num_RVPS, Spe_VNLE, Spe_VPS_List,
    NLP0, NLPx, NLPy, NLPz)

    l2 = 0
    PFp = 0.0
    PFm = 0.0
    for ist = 1:NO0, jst = 1:NO1

        Sum0x_re = 0.0
        Sum0x_im = 0.0
        Sum0y_re = 0.0
        Sum0y_im = 0.0
        Sum0z_re = 0.0
        Sum0z_im = 0.0
        Sum1x_re = 0.0
        Sum1x_im = 0.0
        Sum1y_re = 0.0
        Sum1y_im = 0.0
        Sum1z_re = 0.0
        Sum1z_im = 0.0
        Sum2x_re = 0.0
        Sum2x_im = 0.0
        Sum2y_re = 0.0
        Sum2y_im = 0.0
        Sum2z_re = 0.0
        Sum2z_im = 0.0

        l = 1
        for lnum = 1:Spe_Num_RVPS

            ene_p = Spe_VNLE[1,lnum]
            ene_m = Spe_VNLE[2,lnum]

            if Spe_VPS_List[lnum] == 0
                l2 = 0
                PFp = 1.0
                PFm = 0.0
            elseif Spe_VPS_List[lnum] == 1
                l2 = 2
                PFp = 2/3
                PFm = 1/3
            elseif Spe_VPS_List[lnum] == 2
                l2 = 4
                PFp = 3/5
                PFm = 2/5
            elseif Spe_VPS_List[lnum] == 3
                l2 = 6
                PFp = 4/7
                PFm = 3/7
            end

            
            if l2 == 2
                
                tmp = ene_p/3

                # real contribution of l+1/2 to off diagonal up-down matrix
                tmpx = (tmp*NLPx[ist][1][l  ]*NLP0[jst][1][l+2] 
                       -tmp*NLPx[ist][1][l+2]*NLP0[jst][1][l  ])
                tmpy = (tmp*NLPy[ist][1][l  ]*NLP0[jst][1][l+2] 
                       -tmp*NLPy[ist][1][l+2]*NLP0[jst][1][l  ])
                tmpz = (tmp*NLPz[ist][1][l  ]*NLP0[jst][1][l+2] 
                       -tmp*NLPz[ist][1][l+2]*NLP0[jst][1][l  ])
                
                Sum2x_re += fugou*tmpx
                Sum2y_re += fugou*tmpy
                Sum2z_re += fugou*tmpz

                # imaginary contribution of l+1/2 to off diagonal up-down matrix
                tmpx = (-tmp*NLPx[ist][1][l+1]*NLP0[jst][1][l+2] 
                        +tmp*NLPx[ist][1][l+2]*NLP0[jst][1][l+1])
                tmpy = (-tmp*NLPy[ist][1][l+1]*NLP0[jst][1][l+2] 
                        +tmp*NLPy[ist][1][l+2]*NLP0[jst][1][l+1])
                tmpz = (-tmp*NLPz[ist][1][l+1]*NLP0[jst][1][l+2] 
                        +tmp*NLPz[ist][1][l+2]*NLP0[jst][1][l+1])
                
                Sum2x_im += fugou*tmpx
                Sum2y_im += fugou*tmpy
                Sum2z_im += fugou*tmpz

                tmp = ene_m/3

                # real contribution of l-1/2 for to diagonal up-down matrix
                tmpx = (tmp*NLPx[ist][2][l  ]*NLP0[jst][2][l+2] 
                       -tmp*NLPx[ist][2][l+2]*NLP0[jst][2][l  ])
                tmpy = (tmp*NLPy[ist][2][l  ]*NLP0[jst][2][l+2] 
                       -tmp*NLPy[ist][2][l+2]*NLP0[jst][2][l  ])
                tmpz = (tmp*NLPz[ist][2][l  ]*NLP0[jst][2][l+2] 
                       -tmp*NLPz[ist][2][l+2]*NLP0[jst][2][l  ])
                
                Sum2x_re -= fugou*tmpx
                Sum2y_re -= fugou*tmpy
                Sum2z_re -= fugou*tmpz

                # imaginary contribution of l-1/2 to off diagonal up-down matrix
                tmpx = (-tmp*NLPx[ist][2][l+1]*NLP0[jst][2][l+2] 
                        +tmp*NLPx[ist][2][l+2]*NLP0[jst][2][l+1])
                tmpy = (-tmp*NLPy[ist][2][l+1]*NLP0[jst][2][l+2] 
                        +tmp*NLPy[ist][2][l+2]*NLP0[jst][2][l+1])
                tmpz = (-tmp*NLPz[ist][2][l+1]*NLP0[jst][2][l+2] 
                        +tmp*NLPz[ist][2][l+2]*NLP0[jst][2][l+1])
                
                Sum2x_im -= fugou*tmpx
                Sum2y_im -= fugou*tmpy
                Sum2z_im -= fugou*tmpz
            end

            if l2 == 4

                tmp0 = sqrt(3)
                tmp1 = ene_p/5
                tmp2 = tmp0*tmp1

                # real contribution of l+1/2 to off diagonal up-down matrix
                tmpx = (-tmp2*NLPx[ist][1][l  ]*NLP0[jst][1][l+3] 
                        +tmp2*NLPx[ist][1][l+3]*NLP0[jst][1][l  ]
                        +tmp1*NLPx[ist][1][l+1]*NLP0[jst][1][l+3] 
                        -tmp1*NLPx[ist][1][l+3]*NLP0[jst][1][l+1]
                        +tmp1*NLPx[ist][1][l+2]*NLP0[jst][1][l+4] 
                        -tmp1*NLPx[ist][1][l+4]*NLP0[jst][1][l+2])
                tmpy = (-tmp2*NLPy[ist][1][l  ]*NLP0[jst][1][l+3] 
                        +tmp2*NLPy[ist][1][l+3]*NLP0[jst][1][l  ]
                        +tmp1*NLPy[ist][1][l+1]*NLP0[jst][1][l+3] 
                        -tmp1*NLPy[ist][1][l+3]*NLP0[jst][1][l+1]
                        +tmp1*NLPy[ist][1][l+2]*NLP0[jst][1][l+4] 
                        -tmp1*NLPy[ist][1][l+4]*NLP0[jst][1][l+2])
                tmpz = (-tmp2*NLPz[ist][1][l  ]*NLP0[jst][1][l+3] 
                        +tmp2*NLPz[ist][1][l+3]*NLP0[jst][1][l  ]
                        +tmp1*NLPz[ist][1][l+1]*NLP0[jst][1][l+3] 
                        -tmp1*NLPz[ist][1][l+3]*NLP0[jst][1][l+1]
                        +tmp1*NLPz[ist][1][l+2]*NLP0[jst][1][l+4] 
                        -tmp1*NLPz[ist][1][l+4]*NLP0[jst][1][l+2])
                
                Sum2x_re += fugou*tmpx
                Sum2y_re += fugou*tmpy
                Sum2z_re += fugou*tmpz

                # imaginary contribution of l+1/2 to off diagonal up-down matrix
                tmpx = ( tmp2*NLPx[ist][1][l  ]*NLP0[jst][1][l+4] 
                        -tmp2*NLPx[ist][1][l+4]*NLP0[jst][1][l  ]
                        +tmp1*NLPx[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp1*NLPx[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp1*NLPx[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp1*NLPx[ist][1][l+3]*NLP0[jst][1][l+2])
                tmpy = ( tmp2*NLPy[ist][1][l  ]*NLP0[jst][1][l+4] 
                        -tmp2*NLPy[ist][1][l+4]*NLP0[jst][1][l  ]
                        +tmp1*NLPy[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp1*NLPy[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp1*NLPy[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp1*NLPy[ist][1][l+3]*NLP0[jst][1][l+2])
                tmpz = ( tmp2*NLPz[ist][1][l  ]*NLP0[jst][1][l+4] 
                        -tmp2*NLPz[ist][1][l+4]*NLP0[jst][1][l  ]
                        +tmp1*NLPz[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp1*NLPz[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp1*NLPz[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp1*NLPz[ist][1][l+3]*NLP0[jst][1][l+2])
                
                Sum2x_im += fugou*tmpx
                Sum2y_im += fugou*tmpy
                Sum2z_im += fugou*tmpz


                tmp1 = ene_m/5
                tmp2 = tmp0*tmp1

                # real contribution of l-1/2 for to diagonal up-down matrix
                tmpx = (-tmp2*NLPx[ist][2][l  ]*NLP0[jst][2][l+3] 
                        +tmp2*NLPx[ist][2][l+3]*NLP0[jst][2][l  ]
                        +tmp1*NLPx[ist][2][l+1]*NLP0[jst][2][l+3] 
                        -tmp1*NLPx[ist][2][l+3]*NLP0[jst][2][l+1]
                        +tmp1*NLPx[ist][2][l+2]*NLP0[jst][2][l+4] 
                        -tmp1*NLPx[ist][2][l+4]*NLP0[jst][2][l+2])
                tmpy = (-tmp2*NLPy[ist][2][l  ]*NLP0[jst][2][l+3] 
                        +tmp2*NLPy[ist][2][l+3]*NLP0[jst][2][l  ]
                        +tmp1*NLPy[ist][2][l+1]*NLP0[jst][2][l+3] 
                        -tmp1*NLPy[ist][2][l+3]*NLP0[jst][2][l+1]
                        +tmp1*NLPy[ist][2][l+2]*NLP0[jst][2][l+4] 
                        -tmp1*NLPy[ist][2][l+4]*NLP0[jst][2][l+2])
                tmpz = (-tmp2*NLPz[ist][2][l  ]*NLP0[jst][2][l+3] 
                        +tmp2*NLPz[ist][2][l+3]*NLP0[jst][2][l  ]
                        +tmp1*NLPz[ist][2][l+1]*NLP0[jst][2][l+3] 
                        -tmp1*NLPz[ist][2][l+3]*NLP0[jst][2][l+1]
                        +tmp1*NLPz[ist][2][l+2]*NLP0[jst][2][l+4] 
                        -tmp1*NLPz[ist][2][l+4]*NLP0[jst][2][l+2])
                
                Sum2x_re -= fugou*tmpx
                Sum2y_re -= fugou*tmpy
                Sum2z_re -= fugou*tmpz

                # imaginary contribution of l-1/2 to off diagonal up-down matrix
                tmpx = ( tmp2*NLPx[ist][2][l  ]*NLP0[jst][2][l+4] 
                        -tmp2*NLPx[ist][2][l+4]*NLP0[jst][2][l  ]
                        +tmp1*NLPx[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp1*NLPx[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp1*NLPx[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp1*NLPx[ist][2][l+3]*NLP0[jst][2][l+2])
                tmpy = ( tmp2*NLPy[ist][2][l  ]*NLP0[jst][2][l+4] 
                        -tmp2*NLPy[ist][2][l+4]*NLP0[jst][2][l  ]
                        +tmp1*NLPy[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp1*NLPy[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp1*NLPy[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp1*NLPy[ist][2][l+3]*NLP0[jst][2][l+2])
                tmpz = ( tmp2*NLPz[ist][2][l  ]*NLP0[jst][2][l+4] 
                        -tmp2*NLPz[ist][2][l+4]*NLP0[jst][2][l  ]
                        +tmp1*NLPz[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp1*NLPz[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp1*NLPz[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp1*NLPz[ist][2][l+3]*NLP0[jst][2][l+2])
                
                Sum2x_im -= fugou*tmpx
                Sum2y_im -= fugou*tmpy
                Sum2z_im -= fugou*tmpz
            end

            if l2 == 6
                
                tmp0 = sqrt(6)
                tmp1 = sqrt(3/2)
                tmp2 = sqrt(5/2)

                tmp3 = ene_p/7
                tmp4 = tmp1*tmp3
                tmp5 = tmp2*tmp3
                tmp6 = tmp0*tmp3

                # real contribution of l+1/2 to off diagonal up-down matrix
                
                tmpx = (-tmp6*NLPx[ist][1][l  ]*NLP0[jst][1][l+1] 
                        +tmp6*NLPx[ist][1][l+1]*NLP0[jst][1][l  ]
                        -tmp5*NLPx[ist][1][l+1]*NLP0[jst][1][l+3] 
                        +tmp5*NLPx[ist][1][l+3]*NLP0[jst][1][l+1]
                        -tmp5*NLPx[ist][1][l+2]*NLP0[jst][1][l+4] 
                        +tmp5*NLPx[ist][1][l+4]*NLP0[jst][1][l+2]
                        -tmp4*NLPx[ist][1][l+3]*NLP0[jst][1][l+5] 
                        +tmp4*NLPx[ist][1][l+5]*NLP0[jst][1][l+3]
                        -tmp4*NLPx[ist][1][l+4]*NLP0[jst][1][l+6] 
                        +tmp4*NLPx[ist][1][l+6]*NLP0[jst][1][l+4])
                tmpy = (-tmp6*NLPy[ist][1][l  ]*NLP0[jst][1][l+1] 
                        +tmp6*NLPy[ist][1][l+1]*NLP0[jst][1][l  ]
                        -tmp5*NLPy[ist][1][l+1]*NLP0[jst][1][l+3] 
                        +tmp5*NLPy[ist][1][l+3]*NLP0[jst][1][l+1]
                        -tmp5*NLPy[ist][1][l+2]*NLP0[jst][1][l+4] 
                        +tmp5*NLPy[ist][1][l+4]*NLP0[jst][1][l+2]
                        -tmp4*NLPy[ist][1][l+3]*NLP0[jst][1][l+5] 
                        +tmp4*NLPy[ist][1][l+5]*NLP0[jst][1][l+3]
                        -tmp4*NLPy[ist][1][l+4]*NLP0[jst][1][l+6] 
                        +tmp4*NLPy[ist][1][l+6]*NLP0[jst][1][l+4])
                tmpz = (-tmp6*NLPz[ist][1][l  ]*NLP0[jst][1][l+1] 
                        +tmp6*NLPz[ist][1][l+1]*NLP0[jst][1][l  ]
                        -tmp5*NLPz[ist][1][l+1]*NLP0[jst][1][l+3] 
                        +tmp5*NLPz[ist][1][l+3]*NLP0[jst][1][l+1]
                        -tmp5*NLPz[ist][1][l+2]*NLP0[jst][1][l+4] 
                        +tmp5*NLPz[ist][1][l+4]*NLP0[jst][1][l+2]
                        -tmp4*NLPz[ist][1][l+3]*NLP0[jst][1][l+5] 
                        +tmp4*NLPz[ist][1][l+5]*NLP0[jst][1][l+3]
                        -tmp4*NLPz[ist][1][l+4]*NLP0[jst][1][l+6] 
                        +tmp4*NLPz[ist][1][l+6]*NLP0[jst][1][l+4])
                
                Sum2x_re += fugou*tmpx
                Sum2y_re += fugou*tmpy
                Sum2z_re += fugou*tmpz

                # imaginary contribution of l+1/2 to off diagonal up-down matrix
                tmpx = ( tmp6*NLPx[ist][1][l  ]*NLP0[jst][1][l+2] 
                        -tmp6*NLPx[ist][1][l+2]*NLP0[jst][1][l  ]
                        +tmp5*NLPx[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp5*NLPx[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp5*NLPx[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp5*NLPx[ist][1][l+3]*NLP0[jst][1][l+2]
                        +tmp4*NLPx[ist][1][l+3]*NLP0[jst][1][l+6] 
                        -tmp4*NLPx[ist][1][l+6]*NLP0[jst][1][l+3]
                        -tmp4*NLPx[ist][1][l+4]*NLP0[jst][1][l+5] 
                        +tmp4*NLPx[ist][1][l+5]*NLP0[jst][1][l+4])
                tmpy = ( tmp6*NLPy[ist][1][l  ]*NLP0[jst][1][l+2] 
                        -tmp6*NLPy[ist][1][l+2]*NLP0[jst][1][l  ]
                        +tmp5*NLPy[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp5*NLPy[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp5*NLPy[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp5*NLPy[ist][1][l+3]*NLP0[jst][1][l+2]
                        +tmp4*NLPy[ist][1][l+3]*NLP0[jst][1][l+6] 
                        -tmp4*NLPy[ist][1][l+6]*NLP0[jst][1][l+3]
                        -tmp4*NLPy[ist][1][l+4]*NLP0[jst][1][l+5] 
                        +tmp4*NLPy[ist][1][l+5]*NLP0[jst][1][l+4])
                tmpz = ( tmp6*NLPz[ist][1][l  ]*NLP0[jst][1][l+2] 
                        -tmp6*NLPz[ist][1][l+2]*NLP0[jst][1][l  ]
                        +tmp5*NLPz[ist][1][l+1]*NLP0[jst][1][l+4] 
                        -tmp5*NLPz[ist][1][l+4]*NLP0[jst][1][l+1]
                        -tmp5*NLPz[ist][1][l+2]*NLP0[jst][1][l+3] 
                        +tmp5*NLPz[ist][1][l+3]*NLP0[jst][1][l+2]
                        +tmp4*NLPz[ist][1][l+3]*NLP0[jst][1][l+6] 
                        -tmp4*NLPz[ist][1][l+6]*NLP0[jst][1][l+3]
                        -tmp4*NLPz[ist][1][l+4]*NLP0[jst][1][l+5] 
                        +tmp4*NLPz[ist][1][l+5]*NLP0[jst][1][l+4])
                
                Sum2x_im += fugou*tmpx
                Sum2y_im += fugou*tmpy
                Sum2z_im += fugou*tmpz

                # real contribution of l-1/2 for to diagonal up-down matrix
                tmp3 = ene_m/7
                tmp4 = tmp1*tmp3
                tmp5 = tmp2*tmp3
                tmp6 = tmp0*tmp3

                tmpx = (-tmp6*NLPx[ist][2][l  ]*NLP0[jst][2][l+1] 
                        +tmp6*NLPx[ist][2][l+1]*NLP0[jst][2][l  ]
                        -tmp5*NLPx[ist][2][l+1]*NLP0[jst][2][l+3] 
                        +tmp5*NLPx[ist][2][l+3]*NLP0[jst][2][l+1]
                        -tmp5*NLPx[ist][2][l+2]*NLP0[jst][2][l+4] 
                        +tmp5*NLPx[ist][2][l+4]*NLP0[jst][2][l+2]
                        -tmp4*NLPx[ist][2][l+3]*NLP0[jst][2][l+5] 
                        +tmp4*NLPx[ist][2][l+5]*NLP0[jst][2][l+3]
                        -tmp4*NLPx[ist][2][l+4]*NLP0[jst][2][l+6] 
                        +tmp4*NLPx[ist][2][l+6]*NLP0[jst][2][l+4])
                tmpy = (-tmp6*NLPy[ist][2][l  ]*NLP0[jst][2][l+1] 
                        +tmp6*NLPy[ist][2][l+1]*NLP0[jst][2][l  ]
                        -tmp5*NLPy[ist][2][l+1]*NLP0[jst][2][l+3] 
                        +tmp5*NLPy[ist][2][l+3]*NLP0[jst][2][l+1]
                        -tmp5*NLPy[ist][2][l+2]*NLP0[jst][2][l+4] 
                        +tmp5*NLPy[ist][2][l+4]*NLP0[jst][2][l+2]
                        -tmp4*NLPy[ist][2][l+3]*NLP0[jst][2][l+5] 
                        +tmp4*NLPy[ist][2][l+5]*NLP0[jst][2][l+3]
                        -tmp4*NLPy[ist][2][l+4]*NLP0[jst][2][l+6] 
                        +tmp4*NLPy[ist][2][l+6]*NLP0[jst][2][l+4])
                tmpz = (-tmp6*NLPz[ist][2][l  ]*NLP0[jst][2][l+1] 
                        +tmp6*NLPz[ist][2][l+1]*NLP0[jst][2][l  ]
                        -tmp5*NLPz[ist][2][l+1]*NLP0[jst][2][l+3] 
                        +tmp5*NLPz[ist][2][l+3]*NLP0[jst][2][l+1]
                        -tmp5*NLPz[ist][2][l+2]*NLP0[jst][2][l+4] 
                        +tmp5*NLPz[ist][2][l+4]*NLP0[jst][2][l+2]
                        -tmp4*NLPz[ist][2][l+3]*NLP0[jst][2][l+5] 
                        +tmp4*NLPz[ist][2][l+5]*NLP0[jst][2][l+3]
                        -tmp4*NLPz[ist][2][l+4]*NLP0[jst][2][l+6] 
                        +tmp4*NLPz[ist][2][l+6]*NLP0[jst][2][l+4])
                
                Sum2x_re -= fugou*tmpx
                Sum2y_re -= fugou*tmpy
                Sum2z_re -= fugou*tmpz

                # imaginary contribution of l-1/2 to off diagonal up-down matrix
                tmpx = ( tmp6*NLPx[ist][2][l  ]*NLP0[jst][2][l+2] 
                        -tmp6*NLPx[ist][2][l+2]*NLP0[jst][2][l  ]
                        +tmp5*NLPx[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp5*NLPx[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp5*NLPx[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp5*NLPx[ist][2][l+3]*NLP0[jst][2][l+2]
                        +tmp4*NLPx[ist][2][l+3]*NLP0[jst][2][l+6] 
                        -tmp4*NLPx[ist][2][l+6]*NLP0[jst][2][l+3]
                        -tmp4*NLPx[ist][2][l+4]*NLP0[jst][2][l+5] 
                        +tmp4*NLPx[ist][2][l+5]*NLP0[jst][2][l+4])
                tmpy = ( tmp6*NLPy[ist][2][l  ]*NLP0[jst][2][l+2] 
                        -tmp6*NLPy[ist][2][l+2]*NLP0[jst][2][l  ]
                        +tmp5*NLPy[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp5*NLPy[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp5*NLPy[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp5*NLPy[ist][2][l+3]*NLP0[jst][2][l+2]
                        +tmp4*NLPy[ist][2][l+3]*NLP0[jst][2][l+6] 
                        -tmp4*NLPy[ist][2][l+6]*NLP0[jst][2][l+3]
                        -tmp4*NLPy[ist][2][l+4]*NLP0[jst][2][l+5] 
                        +tmp4*NLPy[ist][2][l+5]*NLP0[jst][2][l+4])
                tmpz = ( tmp6*NLPz[ist][2][l  ]*NLP0[jst][2][l+2] 
                        -tmp6*NLPz[ist][2][l+2]*NLP0[jst][2][l  ]
                        +tmp5*NLPz[ist][2][l+1]*NLP0[jst][2][l+4] 
                        -tmp5*NLPz[ist][2][l+4]*NLP0[jst][2][l+1]
                        -tmp5*NLPz[ist][2][l+2]*NLP0[jst][2][l+3] 
                        +tmp5*NLPz[ist][2][l+3]*NLP0[jst][2][l+2]
                        +tmp4*NLPz[ist][2][l+3]*NLP0[jst][2][l+6] 
                        -tmp4*NLPz[ist][2][l+6]*NLP0[jst][2][l+3]
                        -tmp4*NLPz[ist][2][l+4]*NLP0[jst][2][l+5] 
                        +tmp4*NLPz[ist][2][l+5]*NLP0[jst][2][l+4])
                
                Sum2x_im -= fugou*tmpx
                Sum2y_im -= fugou*tmpy
                Sum2z_im -= fugou*tmpz
            end
            


            if l2 == 2
                
                tmp = ene_p/3

                tmpx = (tmp*NLPx[ist][1][l  ]*NLP0[jst][1][l+1] 
                       -tmp*NLPx[ist][1][l+1]*NLP0[jst][1][l  ])
                tmpy = (tmp*NLPy[ist][1][l  ]*NLP0[jst][1][l+1] 
                       -tmp*NLPy[ist][1][l+1]*NLP0[jst][1][l  ])
                tmpz = (tmp*NLPz[ist][1][l  ]*NLP0[jst][1][l+1] 
                       -tmp*NLPz[ist][1][l+1]*NLP0[jst][1][l  ])
                
                # contribution of l+1/2 for up spin
                Sum0x_im += -fugou*tmpx
                Sum0y_im += -fugou*tmpy
                Sum0z_im += -fugou*tmpz

                # contribution of l+1/2 for down spin
                Sum1x_im += fugou*tmpx
                Sum1y_im += fugou*tmpy
                Sum1z_im += fugou*tmpz


                tmp = ene_m/3

                tmpx = (tmp*NLPx[ist][2][l  ]*NLP0[jst][2][l+1] 
                       -tmp*NLPx[ist][2][l+1]*NLP0[jst][2][l  ])
                tmpy = (tmp*NLPy[ist][2][l  ]*NLP0[jst][2][l+1] 
                       -tmp*NLPy[ist][2][l+1]*NLP0[jst][2][l  ])
                tmpz = (tmp*NLPz[ist][2][l  ]*NLP0[jst][2][l+1] 
                       -tmp*NLPz[ist][2][l+1]*NLP0[jst][2][l  ])
                
                # contribution of l-1/2 for up spin
                Sum0x_im += fugou*tmpx
                Sum0y_im += fugou*tmpy
                Sum0z_im += fugou*tmpz

                # contribution of l-1/2 for down spin
                Sum1x_im += -fugou*tmpx
                Sum1y_im += -fugou*tmpy
                Sum1z_im += -fugou*tmpz
            end

            if l2 == 4

                tmp1 = ene_p*1/5
                tmp2 = ene_p*2/5

                tmpx = ( tmp2*NLPx[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp2*NLPx[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp1*NLPx[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp1*NLPx[ist][1][l+4]*NLP0[jst][1][l+3])
                tmpy = ( tmp2*NLPy[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp2*NLPy[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp1*NLPy[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp1*NLPy[ist][1][l+4]*NLP0[jst][1][l+3])
                tmpz = ( tmp2*NLPz[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp2*NLPz[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp1*NLPz[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp1*NLPz[ist][1][l+4]*NLP0[jst][1][l+3])
                
                # contribution of l+1/2 for up spin
                Sum0x_im += -fugou*tmpx
                Sum0y_im += -fugou*tmpy
                Sum0z_im += -fugou*tmpz

                # contribution of l+1/2 for down spin
                Sum1x_im += fugou*tmpx
                Sum1y_im += fugou*tmpy
                Sum1z_im += fugou*tmpz


                tmp1 = ene_m*1/5
                tmp2 = ene_m*2/5

                tmpx = ( tmp2*NLPx[ist][2][l+1]*NLP0[jst][2][l+2] 
                        -tmp2*NLPx[ist][2][l+2]*NLP0[jst][2][l+1]
                        +tmp1*NLPx[ist][2][l+3]*NLP0[jst][2][l+4] 
                        -tmp1*NLPx[ist][2][l+4]*NLP0[jst][2][l+3])
                tmpy = ( tmp2*NLPy[ist][2][l+1]*NLP0[jst][2][l+2] 
                        -tmp2*NLPy[ist][2][l+2]*NLP0[jst][2][l+1]
                        +tmp1*NLPy[ist][2][l+3]*NLP0[jst][2][l+4] 
                        -tmp1*NLPy[ist][2][l+4]*NLP0[jst][2][l+3])
                tmpz = ( tmp2*NLPz[ist][2][l+1]*NLP0[jst][2][l+2] 
                        -tmp2*NLPz[ist][2][l+2]*NLP0[jst][2][l+1]
                        +tmp1*NLPz[ist][2][l+3]*NLP0[jst][2][l+4] 
                        -tmp1*NLPz[ist][2][l+4]*NLP0[jst][2][l+3])
                
                # contribution of l-1/2 for up spin
                Sum0x_im += fugou*tmpx
                Sum0y_im += fugou*tmpy
                Sum0z_im += fugou*tmpz

                # contribution of l-1/2 for down spin
                Sum1x_im += -fugou*tmpx
                Sum1y_im += -fugou*tmpy
                Sum1z_im += -fugou*tmpz
            end

            if l2 == 6

                tmp1 = ene_p/7
                tmp2 = ene_p*2/7
                tmp3 = ene_p*3/7

                tmpx = ( tmp1*NLPx[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPx[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPx[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPx[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPx[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPx[ist][1][l+6]*NLP0[jst][1][l+5])
                tmpy = ( tmp1*NLPy[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPy[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPy[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPy[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPy[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPy[ist][1][l+6]*NLP0[jst][1][l+5])
                tmpz = ( tmp1*NLPz[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPz[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPz[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPz[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPz[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPz[ist][1][l+6]*NLP0[jst][1][l+5])
                
                # contribution of l+1/2 for up spin
                Sum0x_im += -fugou*tmpx
                Sum0y_im += -fugou*tmpy
                Sum0z_im += -fugou*tmpz

                # contribution of l+1/2 for down spin
                Sum1x_im += fugou*tmpx
                Sum1y_im += fugou*tmpy
                Sum1z_im += fugou*tmpz

                tmp1 = ene_m/7
                tmp2 = ene_m*2/7
                tmp3 = ene_m*3/7

                tmpx = ( tmp1*NLPx[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPx[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPx[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPx[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPx[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPx[ist][1][l+6]*NLP0[jst][1][l+5])
                tmpy = ( tmp1*NLPy[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPy[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPy[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPy[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPy[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPy[ist][1][l+6]*NLP0[jst][1][l+5])
                tmpz = ( tmp1*NLPz[ist][1][l+1]*NLP0[jst][1][l+2] 
                        -tmp1*NLPz[ist][1][l+2]*NLP0[jst][1][l+1]
                        +tmp2*NLPz[ist][1][l+3]*NLP0[jst][1][l+4] 
                        -tmp2*NLPz[ist][1][l+4]*NLP0[jst][1][l+3]
                        +tmp3*NLPz[ist][1][l+5]*NLP0[jst][1][l+6] 
                        -tmp3*NLPz[ist][1][l+6]*NLP0[jst][1][l+5])
                
                # contribution of l-1/2 for up spin
                Sum0x_im += fugou*tmpx
                Sum0y_im += fugou*tmpy
                Sum0z_im += fugou*tmpz

                # contribution of l-1/2 for down spin
                Sum1x_im += -fugou*tmpx
                Sum1y_im += -fugou*tmpy
                Sum1z_im += -fugou*tmpz
            end
            

            for lll = 0:l2

                # VNL for j=l+1/2
                tmpx = PFp*ene_p*NLPx[ist][1][l]*NLP0[jst][1][l]
                tmpy = PFp*ene_p*NLPy[ist][1][l]*NLP0[jst][1][l]
                tmpz = PFp*ene_p*NLPz[ist][1][l]*NLP0[jst][1][l]
                
                Sum0x_re += tmpx
                Sum0y_re += tmpy
                Sum0z_re += tmpz

                Sum1x_re += tmpx
                Sum1y_re += tmpy
                Sum1z_re += tmpz

                # VNL for j=l-1/2
                tmpx = PFm*ene_m*NLPx[ist][2][l]*NLP0[jst][2][l]
                tmpy = PFm*ene_m*NLPy[ist][2][l]*NLP0[jst][2][l]
                tmpz = PFm*ene_m*NLPz[ist][2][l]*NLP0[jst][2][l]
                
                Sum0x_re += tmpx
                Sum0y_re += tmpy
                Sum0z_re += tmpz

                Sum1x_re += tmpx
                Sum1y_re += tmpy
                Sum1z_re += tmpz

                l += 1
            end
        end
        

        if orb_switch == 0
            Hx[ist,jst,1] += Sum0x_re + im*Sum0x_im
            Hx[ist,jst,2] += Sum1x_re + im*Sum1x_im
            Hx[ist,jst,3] += Sum2x_re + im*Sum2x_im

            Hy[ist,jst,1] += Sum0y_re + im*Sum0y_im
            Hy[ist,jst,2] += Sum1y_re + im*Sum1y_im
            Hy[ist,jst,3] += Sum2y_re + im*Sum2y_im

            Hz[ist,jst,1] += Sum0z_re + im*Sum0z_im
            Hz[ist,jst,2] += Sum1z_re + im*Sum1z_im
            Hz[ist,jst,3] += Sum2z_re + im*Sum2z_im
        else
            Hx[jst,ist,1] += Sum0x_re + im*Sum0x_im
            Hx[jst,ist,2] += Sum1x_re + im*Sum1x_im
            Hx[jst,ist,3] += Sum2x_re + im*Sum2x_im

            Hy[jst,ist,1] += Sum0y_re + im*Sum0y_im
            Hy[jst,ist,2] += Sum1y_re + im*Sum1y_im
            Hy[jst,ist,3] += Sum2y_re + im*Sum2y_im

            Hz[jst,ist,1] += Sum0z_re + im*Sum0z_im
            Hz[jst,ist,2] += Sum1z_re + im*Sum1z_im
            Hz[jst,ist,3] += Sum2z_re + im*Sum2z_im
        end
    end
end


function _Calc_Force3_8(NO0, NO1, NumOLG, GridListAtom, GListTAtoms1, GListTAtoms2, dOrbs_Gridx, dOrbs_Gridy, dOrbs_Gridz, Orbs_Grid, Vpot_Grid, DM)

    Sumx = 0.0
    Sumy = 0.0
    Sumz = 0.0
	for Nog = 1:8:NumOLG-7

		Nc0 = GListTAtoms1[Nog]+1
		Nc1 = GListTAtoms1[Nog+1]+1
		Nc2 = GListTAtoms1[Nog+2]+1
		Nc3 = GListTAtoms1[Nog+3]+1
		Nc4 = GListTAtoms1[Nog+4]+1
		Nc5 = GListTAtoms1[Nog+5]+1
		Nc6 = GListTAtoms1[Nog+6]+1
		Nc7 = GListTAtoms1[Nog+7]+1

		MN0 = GridListAtom[Nc0]+1
		MN1 = GridListAtom[Nc1]+1
		MN2 = GridListAtom[Nc2]+1
		MN3 = GridListAtom[Nc3]+1
		MN4 = GridListAtom[Nc4]+1
		MN5 = GridListAtom[Nc5]+1
		MN6 = GridListAtom[Nc6]+1
		MN7 = GridListAtom[Nc7]+1

		Nh0 = GListTAtoms2[Nog]+1
		Nh1 = GListTAtoms2[Nog+1]+1
		Nh2 = GListTAtoms2[Nog+2]+1
		Nh3 = GListTAtoms2[Nog+3]+1
		Nh4 = GListTAtoms2[Nog+4]+1
		Nh5 = GListTAtoms2[Nog+5]+1
		Nh6 = GListTAtoms2[Nog+6]+1
		Nh7 = GListTAtoms2[Nog+7]+1
		
		temp0 = Vpot_Grid[MN0]
		temp1 = Vpot_Grid[MN1]
		temp2 = Vpot_Grid[MN2]
		temp3 = Vpot_Grid[MN3]
		temp4 = Vpot_Grid[MN4]
		temp5 = Vpot_Grid[MN5]
		temp6 = Vpot_Grid[MN6]
		temp7 = Vpot_Grid[MN7]

		for ist = 1:NO0
			Sum0x = temp0 * dOrbs_Gridx[Nc0][ist]
			Sum1x = temp1 * dOrbs_Gridx[Nc1][ist]
			Sum2x = temp2 * dOrbs_Gridx[Nc2][ist]
			Sum3x = temp3 * dOrbs_Gridx[Nc3][ist]
			Sum4x = temp4 * dOrbs_Gridx[Nc4][ist]
			Sum5x = temp5 * dOrbs_Gridx[Nc5][ist]
			Sum6x = temp6 * dOrbs_Gridx[Nc6][ist]
			Sum7x = temp7 * dOrbs_Gridx[Nc7][ist]
            
            Sum0y = temp0 * dOrbs_Gridy[Nc0][ist]
			Sum1y = temp1 * dOrbs_Gridy[Nc1][ist]
			Sum2y = temp2 * dOrbs_Gridy[Nc2][ist]
			Sum3y = temp3 * dOrbs_Gridy[Nc3][ist]
			Sum4y = temp4 * dOrbs_Gridy[Nc4][ist]
			Sum5y = temp5 * dOrbs_Gridy[Nc5][ist]
			Sum6y = temp6 * dOrbs_Gridy[Nc6][ist]
			Sum7y = temp7 * dOrbs_Gridy[Nc7][ist]

            Sum0z = temp0 * dOrbs_Gridz[Nc0][ist]
			Sum1z = temp1 * dOrbs_Gridz[Nc1][ist]
			Sum2z = temp2 * dOrbs_Gridz[Nc2][ist]
			Sum3z = temp3 * dOrbs_Gridz[Nc3][ist]
			Sum4z = temp4 * dOrbs_Gridz[Nc4][ist]
			Sum5z = temp5 * dOrbs_Gridz[Nc5][ist]
			Sum6z = temp6 * dOrbs_Gridz[Nc6][ist]
			Sum7z = temp7 * dOrbs_Gridz[Nc7][ist]
			for jst = 1:NO1
                DM_tmp = DM[ist][jst]
                orb0 = Orbs_Grid[Nh0][jst]*DM_tmp
                orb1 = Orbs_Grid[Nh1][jst]*DM_tmp
                orb2 = Orbs_Grid[Nh2][jst]*DM_tmp
                orb3 = Orbs_Grid[Nh3][jst]*DM_tmp
                orb4 = Orbs_Grid[Nh4][jst]*DM_tmp
                orb5 = Orbs_Grid[Nh5][jst]*DM_tmp
                orb6 = Orbs_Grid[Nh6][jst]*DM_tmp
                orb7 = Orbs_Grid[Nh7][jst]*DM_tmp
 
				Sumx += Sum0x * orb0
				Sumx += Sum1x * orb1
				Sumx += Sum2x * orb2
				Sumx += Sum3x * orb3
				Sumx += Sum4x * orb4
				Sumx += Sum5x * orb5
				Sumx += Sum6x * orb6
				Sumx += Sum7x * orb7

                Sumy += Sum0y * orb0
				Sumy += Sum1y * orb1
				Sumy += Sum2y * orb2
				Sumy += Sum3y * orb3
				Sumy += Sum4y * orb4
				Sumy += Sum5y * orb5
				Sumy += Sum6y * orb6
				Sumy += Sum7y * orb7

                Sumz += Sum0z * orb0
				Sumz += Sum1z * orb1
				Sumz += Sum2z * orb2
				Sumz += Sum3z * orb3
				Sumz += Sum4z * orb4
				Sumz += Sum5z * orb5
				Sumz += Sum6z * orb6
				Sumz += Sum7z * orb7
			end
		end
	end


	Nog1 = 8*div(NumOLG, 8)
	rem_NumOLG = rem(NumOLG, 8)
    for Nog = 1:rem_NumOLG
        Nc = GListTAtoms1[Nog+Nog1]+1
		MN = GridListAtom[Nc]+1
		Nh = GListTAtoms2[Nog+Nog1]+1

        temp = Vpot_Grid[MN]
        for ist = 1:NO0
            Sum1x = temp * dOrbs_Gridx[Nc][ist]
            Sum1y = temp * dOrbs_Gridy[Nc][ist]
            Sum1z = temp * dOrbs_Gridz[Nc][ist]
            for jst = 1:NO1
                DM_tmp = DM[ist][jst]
                Orbs_temp = Orbs_Grid[Nh][jst]
                Sumx += Sum1x * Orbs_temp * DM_tmp
                Sumy += Sum1y * Orbs_temp * DM_tmp
                Sumz += Sum1z * Orbs_temp * DM_tmp
            end
        end
    end

    return Sumx, Sumy, Sumz
end