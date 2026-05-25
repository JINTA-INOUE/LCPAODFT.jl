@timeit timer "Set_AdenPCC_Grid" function Set_AdenPCC_Grid(SpinPol::AbstractString, Init_Atoms_Nspin, Init_Atoms_Angle, pao::Vector{PAO}, pspot::Vector{Pspot}, ucell::UCell)

    system_grid = ucell.system_grid
    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    atv = system_grid.atv
    Gxyz = system_grid.Gxyz
    Grid_Origin = system_grid.Grid_Origin
    Ngrid = system_grid.Ngrid
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)
    
    GridN_Atom = ucell.GridN_Atom
    GridListAtom = ucell.GridListAtom
    CellListAtom = ucell.CellListAtom


    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end


    gLatvecs = zeros(Float64, 3, 3)
    gLatvecs[1,:] = Latvecs[1,:]/Ngrid1
	gLatvecs[2,:] = Latvecs[2,:]/Ngrid2
	gLatvecs[3,:] = Latvecs[3,:]/Ngrid3
         

    ADensity_Grid = zeros(Float64, NN)
    PCCDensity_Grid = zeros(Float64, NN)

    Density_Grid = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        Density_Grid[spin] = zeros(Float64, NN)
    end


    for atom = 1:Natom

        spe = atom2spe[atom]
        Spe_Num_Mesh_PAO = pao[spe].Spe_Num_Mesh_PAO
        Spe_PAO_RV = pao[spe].Spe_PAO_RV
        Spe_PAO_XV = pao[spe].Spe_PAO_XV
        Spe_Atomic_Den = pao[spe].Spe_Atomic_Den

        Spe_Num_Mesh_VPS = pspot[spe].Spe_Num_Mesh_VPS
        Spe_VPS_RV = pspot[spe].Spe_VPS_RV
        Spe_VPS_XV = pspot[spe].Spe_VPS_XV
        Spe_Atomic_PCC = pspot[spe].Spe_Atomic_PCC
        
        Nu = Init_Atoms_Nspin[atom][1]
        Nd = Init_Atoms_Nspin[atom][2]
        Nele = Nu + Nd
        ocupcy_u = ifelse(Nele>1e-15, Nu/Nele, 0.0)
        ocupcy_p = ifelse(Nele>1e-15, Nd/Nele, 0.0)

        theta = Init_Atoms_Angle[atom][1]
        phi = Init_Atoms_Angle[atom][2]
        

        @inbounds for xyz = 1:GridN_Atom[atom]
            GNc = GridListAtom[atom][xyz]
            GRc = CellListAtom[atom][xyz]+1

            n1 = div(GNc, Ngrid2*Ngrid3)
            n2 = div(GNc - n1*Ngrid2*Ngrid3, Ngrid3)
            n3 = GNc - n1*Ngrid2*Ngrid3 - n2*Ngrid3
            
            Cx = n1*gLatvecs[1,1] + n2*gLatvecs[2,1] + n3*gLatvecs[3,1] + Grid_Origin[1]
            Cy = n1*gLatvecs[1,2] + n2*gLatvecs[2,2] + n3*gLatvecs[3,2] + Grid_Origin[2]
            Cz = n1*gLatvecs[1,3] + n2*gLatvecs[2,3] + n3*gLatvecs[3,3] + Grid_Origin[3]

            dx = Cx + atv[GRc][1] - Gxyz[atom][1]
            dy = Cy + atv[GRc][2] - Gxyz[atom][2]
            dz = Cz + atv[GRc][3] - Gxyz[atom][3]
            x = 0.5*log(dx^2 + dy^2 + dz^2)

            ADen = KumoF(Spe_Num_Mesh_PAO, x, Spe_PAO_XV, Spe_PAO_RV, Spe_Atomic_Den)
            PCCDen = KumoF(Spe_Num_Mesh_VPS, x, Spe_VPS_XV, Spe_VPS_RV, Spe_Atomic_PCC)

            ADensity_Grid[GNc+1] += 0.5*ADen
            PCCDensity_Grid[GNc+1] += 0.5*PCCDen

            if SpinPol == "off"
                Density_Grid[1][GNc+1] += 0.5 * ADen
            elseif SpinPol == "on"
                Density_Grid[1][GNc+1] += ocupcy_u * ADen
                Density_Grid[2][GNc+1] += ocupcy_p * ADen
            elseif SpinPol == "nc"
                mag = (ocupcy_u-ocupcy_p)*ADen
                Density_Grid[1][GNc+1] += ADen
                Density_Grid[2][GNc+1] += mag*sin(theta)*cos(phi)
                Density_Grid[3][GNc+1] += mag*sin(theta)*sin(phi)
                Density_Grid[4][GNc+1] += mag*cos(theta)
            end
        end
    end


    if SpinPol == "nc"
        @inbounds for i = 1:NN
            rho  = Density_Grid[1][i]
            magx = Density_Grid[2][i]
            magy = Density_Grid[3][i]
            magz = Density_Grid[4][i]

            Density_Grid[1][i] =  0.5*(rho + magz)
            Density_Grid[2][i] =  0.5*(rho - magz)
            Density_Grid[3][i] =  0.5*magx
            Density_Grid[4][i] = -0.5*magy
        end
    end


    return ADensity_Grid, PCCDensity_Grid, Density_Grid
end
