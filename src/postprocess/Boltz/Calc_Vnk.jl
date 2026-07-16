@timeit timer "Calc_Vnk!" function Calc_Vnk!(material::CWF_model, kpoints::KPoints, Enk, Cnk)

    Latvecs = material.Latvecs
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    ChemP = material.ChemP
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR
    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh


    Nkpt = kpoints.Nkpt
    kpts = kpoints.MPI_kpts
    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    cartesian_cell = Vector{Vector{Float64}}(undef, NCell)
    for cell = 1:NCell
        cartesian_cell[cell] = zeros(Float64, 3)
        cartesian_cell[cell][1] = Latvecs[1,1]*cell_list_ijk[cell][1] + Latvecs[2,1]*cell_list_ijk[cell][2] + Latvecs[3,1]*cell_list_ijk[cell][3]
        cartesian_cell[cell][2] = Latvecs[1,2]*cell_list_ijk[cell][1] + Latvecs[2,2]*cell_list_ijk[cell][2] + Latvecs[3,2]*cell_list_ijk[cell][3]
        cartesian_cell[cell][3] = Latvecs[1,3]*cell_list_ijk[cell][1] + Latvecs[2,3]*cell_list_ijk[cell][2] + Latvecs[3,3]*cell_list_ijk[cell][3]
        cartesian_cell[cell][1] = cartesian_cell[cell][1]/Ang_to_bohr
        cartesian_cell[cell][2] = cartesian_cell[cell][2]/Ang_to_bohr
        cartesian_cell[cell][3] = cartesian_cell[cell][3]/Ang_to_bohr
    end



    for spin = 1:spinsize, cell = 1:NCell, ist = 1:Nwann, jst = 1:Nwann
        HmnR[jst,ist,cell,spin] = HmnR[jst,ist,cell,spin]*eV2Hartree
    end


    for spin = 1:spinsize, cell = 1:NCell, ist = 1:Nwann
        if cell_list_ijk[cell] == [0,0,0]
            HmnR[ist,ist,cell,spin] = HmnR[ist,ist,cell,spin] - ChemP*eV2Hartree
        end
    end
    

    Vnk = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, spinsize)
    for spin = 1:spinsize
        Vnk[spin] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Vnk[spin][ik] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Vnk[spin][ik][jk] = Vector{Vector{Vector{Float64}}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Vnk[spin][ik][jk][kk] = Vector{Vector{Float64}}(undef, 3)
                    for xyz = 1:3
                        Vnk[spin][ik][jk][kk][xyz] = zeros(Float64, Nwann)
                    end
                end
            end
        end
    end

    


    Hmnk_x = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_y = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_z = zeros(ComplexF64, Nwann, Nwann)
    Htemp1 = zeros(ComplexF64, Nwann, Nwann)
    Htemp2 = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, k = 1:Nkpt

        fill!(Hmnk_x, 0.0)
        fill!(Hmnk_y, 0.0)
        fill!(Hmnk_z, 0.0)
            
        for cell = 1:NCell
            kRn = dot(kpts[k], cell_list_ijk[cell])
            phase = cispi(2*kRn)

            tmpx = -im*cartesian_cell[cell][1]
            tmpy = -im*cartesian_cell[cell][2]
            tmpz = -im*cartesian_cell[cell][3]

            @inbounds for ist = 1:Nwann, jst = 1:Nwann
                Hmnk_x[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpx
                Hmnk_y[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpy
                Hmnk_z[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpz
            end
        end

        ik, jk, kk = kindex[k,:]

        mul!(Htemp1, Cnk[spin][k]', Hmnk_x)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][1] = real(diag(Htemp2))
        mul!(Htemp1, Cnk[spin][k]', Hmnk_y)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][2] = real(diag(Htemp2))
        mul!(Htemp1, Cnk[spin][k]', Hmnk_z)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][3] = real(diag(Htemp2))
    end


    # change unit (Hartree -> eV)
    @inbounds for spin = 1:spinsize, ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3, μ = 1:Nwann
        Enk[spin][ik][jk][kk][μ] = Enk[spin][ik][jk][kk][μ]*eV2Hartree - ChemP*eV2Hartree
    end


    return Vnk
end


@inline @timeit timer "HS_matrix_iR!" function HS_matrix_iR_Vnk!(Sx, Sy, Sz, Hx, Hy, Hz, OLP, Hks, Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, kpts::Vector{Float64}, Rvec)
    ka, kb, kc = kpts
    fill!(Sx, 0.0)
    fill!(Hx, 0.0)
    fill!(Sy, 0.0)
    fill!(Hy, 0.0)
    fill!(Sz, 0.0)
    fill!(Hz, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        phase_x = cispi(2*kRn)*Rvec[cell][1]*im
        phase_y = cispi(2*kRn)*Rvec[cell][2]*im
        phase_z = cispi(2*kRn)*Rvec[cell][3]*im
        @inbounds for ist = 1:NO0, jst = 1:NO1
            Sx[Anum+ist, Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_x
            Hx[Anum+ist, Bnum+jst] += Hks[atom][Rn][ist][jst]*phase_x
            Sy[Anum+ist, Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_y
            Hy[Anum+ist, Bnum+jst] += Hks[atom][Rn][ist][jst]*phase_y
            Sz[Anum+ist, Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_z
            Hz[Anum+ist, Bnum+jst] += Hks[atom][Rn][ist][jst]*phase_z
        end
    end
end


@inline @timeit timer "HS_matrix_iR!" function HS_matrix_iR_Vnk!(Sx, Sy, Sz, OLP, Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, kpts::Vector{Float64}, Rvec)
    ka, kb, kc = kpts
    fill!(Sx, 0.0)
    fill!(Sy, 0.0)
    fill!(Sz, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        phase_x = cispi(2*kRn)*Rvec[cell][1]*im
        phase_y = cispi(2*kRn)*Rvec[cell][2]*im
        phase_z = cispi(2*kRn)*Rvec[cell][3]*im
        @inbounds for ist = 1:NO0, jst = 1:NO1
            Sx[Anum+ist,Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_x
            Sy[Anum+ist,Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_y
            Sz[Anum+ist,Bnum+jst] += OLP[atom][Rn][ist][jst]*phase_z
        end
    end
end


@inline @timeit timer "HS_matrix_iR!" function HS_matrix_NC_iR_Vnk!(
    Hx_uu, Hx_dd, Hx_ud, Hx_du,
    Hy_uu, Hy_dd, Hy_ud, Hy_du,
    Hz_uu, Hz_dd, Hz_ud, Hz_du, 
    Hks, iHks, 
    Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, 
    kpts::Vector{Float64}, Rvec)
    
    ka, kb, kc = kpts
    
    fill!(Hx_uu, 0.0)
    fill!(Hx_dd, 0.0)
    fill!(Hx_ud, 0.0)
    fill!(Hx_du, 0.0)
    fill!(Hy_uu, 0.0)
    fill!(Hy_dd, 0.0)
    fill!(Hy_ud, 0.0)
    fill!(Hy_du, 0.0)
    fill!(Hz_uu, 0.0)
    fill!(Hz_dd, 0.0)
    fill!(Hz_ud, 0.0)
    fill!(Hz_du, 0.0)
    Hks1 = Hks[1]
    Hks2 = Hks[2]
    Hks3 = Hks[3]
    Hks4 = Hks[4]
    iHks1 = iHks[1]
    iHks2 = iHks[2]
    iHks3 = iHks[3]
    
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
        kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
        phase_x = cispi(2*kRn)*Rvec[cell][1]*im
        phase_y = cispi(2*kRn)*Rvec[cell][2]*im
        phase_z = cispi(2*kRn)*Rvec[cell][3]*im
        _Hks_uu = Hks1[atom][Rn]
        _Hks_dd = Hks2[atom][Rn]
        _Hks_ud = Hks3[atom][Rn]
        _iHks_ud1 = Hks4[atom][Rn]
        _iHks_uu = iHks1[atom][Rn]
        _iHks_dd = iHks2[atom][Rn]
        _iHks_ud2 = iHks3[atom][Rn]
        @inbounds for ist = 1:NO0, jst = 1:NO1
            Hks_uu = _Hks_uu[ist][jst]
            Hks_dd = _Hks_dd[ist][jst]
            Hks_ud = _Hks_ud[ist][jst]
            iHks_ud1 = _iHks_ud1[ist][jst]
            iHks_uu = _iHks_uu[ist][jst]
            iHks_dd = _iHks_dd[ist][jst]
            iHks_ud2 = _iHks_ud2[ist][jst]
            Hx_uu[Anum+ist,Bnum+jst] += (Hks_uu + im*iHks_uu)*phase_x
            Hx_dd[Anum+ist,Bnum+jst] += (Hks_dd + im*iHks_dd)*phase_x
            Hx_ud[Anum+ist,Bnum+jst] += (Hks_ud + im*(iHks_ud1 + iHks_ud2))*phase_x
            Hx_du[Anum+ist,Bnum+jst] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_x
            # Hx_du[Bnum+jst,Anum+ist] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_x
            Hy_uu[Anum+ist,Bnum+jst] += (Hks_uu + im*iHks_uu)*phase_y
            Hy_dd[Anum+ist,Bnum+jst] += (Hks_dd + im*iHks_dd)*phase_y
            Hy_ud[Anum+ist,Bnum+jst] += (Hks_ud + im*(iHks_ud1 + iHks_ud2))*phase_y
            # Hy_du[Bnum+jst,Anum+ist] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_y
            Hy_du[Anum+ist,Bnum+jst] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_y
            Hz_uu[Anum+ist,Bnum+jst] += (Hks_uu + im*iHks_uu)*phase_z
            Hz_dd[Anum+ist,Bnum+jst] += (Hks_dd + im*iHks_dd)*phase_z
            Hz_ud[Anum+ist,Bnum+jst] += (Hks_ud + im*(iHks_ud1 + iHks_ud2))*phase_z
            # Hz_du[Bnum+jst,Anum+ist] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_z
            Hz_du[Anum+ist,Bnum+jst] += (Hks_ud - im*(iHks_ud1 + iHks_ud2))*phase_z
        end
    end
end


@timeit timer "Calc_Vnk!" function Calc_Vnk!(material::LCPAO_model, kpoints::KPoints, Enk, Cnk)

    TCpyCell = material.TCpyCell
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Latvecs = material.Latvecs
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol=="nc", 2*fsize, fsize)
    ChemP = material.ChemP

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh


    cartesian_cell = Vector{Vector{Float64}}(undef, TCpyCell)
    for cell = 1:TCpyCell
        cartesian_cell[cell] = zeros(Float64, 3)
        cartesian_cell[cell][1] = Latvecs[1,1]*atv_ijk[cell][1] + Latvecs[2,1]*atv_ijk[cell][2] + Latvecs[3,1]*atv_ijk[cell][3]
        cartesian_cell[cell][2] = Latvecs[1,2]*atv_ijk[cell][1] + Latvecs[2,2]*atv_ijk[cell][2] + Latvecs[3,2]*atv_ijk[cell][3]
        cartesian_cell[cell][3] = Latvecs[1,3]*atv_ijk[cell][1] + Latvecs[2,3]*atv_ijk[cell][2] + Latvecs[3,3]*atv_ijk[cell][3]
        cartesian_cell[cell][1] = cartesian_cell[cell][1]/Ang_to_bohr
        cartesian_cell[cell][2] = cartesian_cell[cell][2]/Ang_to_bohr
        cartesian_cell[cell][3] = cartesian_cell[cell][3]/Ang_to_bohr
    end


    
    # change unit
    @inbounds for spin = 1:spinsize, ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3, μ = 1:Nfsize
        Enk[spin][ik][jk][kk][μ] = Enk[spin][ik][jk][kk][μ]*eV2Hartree
    end
    

    Vnk = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, spinsize)
    for spin = 1:spinsize
        Vnk[spin] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Vnk[spin][ik] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Vnk[spin][ik][jk] = Vector{Vector{Vector{Float64}}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Vnk[spin][ik][jk][kk] = Vector{Vector{Float64}}(undef, 3)
                    for xyz = 1:3
                        Vnk[spin][ik][jk][kk][xyz] = zeros(Float64, Nfsize)
                    end
                end
            end
        end
    end

    if SpinPol ∈ ("off", "on")
        Calc_Vnk_Col!(material, kpoints, cartesian_cell, Enk, Cnk, Vnk)
    elseif SpinPol == "nc"
        Calc_Vnk_NonCol!(material, kpoints, cartesian_cell, Enk, Cnk, Vnk)
    end



    # Band shift
    @inbounds for spin = 1:spinsize, ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3, μ = 1:Nfsize
        Enk[spin][ik][jk][kk][μ] = Enk[spin][ik][jk][kk][μ] - ChemP*eV2Hartree
    end


    return Vnk
end


function Calc_Vnk_Col!(material::LCPAO_model, kpoints::KPoints, cartesian_cell, Enk, Cnk, Vnk)

    Nspin = material.Nspin
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    fsize = sum(Total_NumOrbs)
    MP = material.MP
    OLP = material.OLP
    Hks = material.Hks

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    kpts = kpoints.MPI_kpts

    for spin = 1:Nspin, atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        Hks[spin][atom][Rn][ist][jst] = Hks[spin][atom][Rn][ist][jst]*eV2Hartree
    end


    Hx = zeros(ComplexF64, fsize, fsize)
    Hy = zeros(ComplexF64, fsize, fsize)
    Hz = zeros(ComplexF64, fsize, fsize)
    Sx = zeros(ComplexF64, fsize, fsize)
    Sy = zeros(ComplexF64, fsize, fsize)
    Sz = zeros(ComplexF64, fsize, fsize)
    Cnk_temp1 = zeros(ComplexF64, fsize, fsize)
    Cnk_temp2 = zeros(ComplexF64, fsize, fsize)
    Htemp1 = zeros(ComplexF64, fsize, fsize)
    Htemp2 = zeros(ComplexF64, fsize, fsize)


    for spin = 1:spinsize
        k = 0
        for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            k += 1
            HS_matrix_iR_Vnk!(Sx, Sy, Sz, Hx, Hy, Hz, OLP, Hks[spin], Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, kpts[k], cartesian_cell) 
            
            @. Cnk_temp1 = Cnk[spin][k]'
            @. Cnk_temp2 = Cnk[spin][k]
                
            mul!(Htemp1, Cnk_temp1, Hx)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            Vnk[spin][ik][jk][kk][1] = real(diag(Htemp2))
            mul!(Htemp1, Cnk_temp1, Hy)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            Vnk[spin][ik][jk][kk][2] = real(diag(Htemp2))
            mul!(Htemp1, Cnk_temp1, Hz)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            Vnk[spin][ik][jk][kk][3] = real(diag(Htemp2))

            mul!(Htemp1, Cnk_temp1, Sx)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            @inbounds for μ = 1:fsize
                Vnk[spin][ik][jk][kk][1][μ] -= real(Htemp2[μ,μ])*Enk[spin][ik][jk][kk][μ]
            end
            mul!(Htemp1, Cnk_temp1, Sy)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            @inbounds for μ = 1:fsize
                Vnk[spin][ik][jk][kk][2][μ] -= real(Htemp2[μ,μ])*Enk[spin][ik][jk][kk][μ]
            end
            mul!(Htemp1, Cnk_temp1, Sz)
            mul!(Htemp2, Htemp1, Cnk_temp2)
            @inbounds for μ = 1:fsize
                Vnk[spin][ik][jk][kk][3][μ] -= real(Htemp2[μ,μ])*Enk[spin][ik][jk][kk][μ]
            end
        end
    end
end


function Calc_Vnk_NonCol!(material::LCPAO_model, kpoints::KPoints, cartesian_cell, Enk, Cnk, Vnk)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize
    ChemP = material.ChemP
    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks

    kmesh = kpoints.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    kpts = kpoints.MPI_kpts
    

    for spin = 1:4, atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        Hks[spin][atom][Rn][ist][jst] = Hks[spin][atom][Rn][ist][jst]*eV2Hartree
    end

    for spin = 1:3, atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[natn[atom][Rn]]
        iHks[spin][atom][Rn][ist][jst] = iHks[spin][atom][Rn][ist][jst]*eV2Hartree
    end


    Hx_uu = zeros(ComplexF64, fsize, fsize)
    Hx_dd = zeros(ComplexF64, fsize, fsize)
    Hx_ud = zeros(ComplexF64, fsize, fsize)
    Hx_du = zeros(ComplexF64, fsize, fsize)
    Hy_uu = zeros(ComplexF64, fsize, fsize)
    Hy_dd = zeros(ComplexF64, fsize, fsize)
    Hy_ud = zeros(ComplexF64, fsize, fsize)
    Hy_du = zeros(ComplexF64, fsize, fsize)
    Hz_uu = zeros(ComplexF64, fsize, fsize)
    Hz_dd = zeros(ComplexF64, fsize, fsize)
    Hz_ud = zeros(ComplexF64, fsize, fsize)
    Hz_du = zeros(ComplexF64, fsize, fsize)
    Sx = zeros(ComplexF64, fsize, fsize)
    Sy = zeros(ComplexF64, fsize, fsize)
    Sz = zeros(ComplexF64, fsize, fsize)
    Cnk_up_temp1 = zeros(ComplexF64, Nfsize, fsize)
    Cnk_dn_temp1 = zeros(ComplexF64, Nfsize, fsize)
    Cnk_up_temp2 = zeros(ComplexF64, fsize, Nfsize)
    Cnk_dn_temp2 = zeros(ComplexF64, fsize, Nfsize)
    Htemp1 = zeros(ComplexF64, Nfsize, fsize)
    Htemp2 = zeros(ComplexF64, Nfsize, Nfsize)

    
    k = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        k += 1
        HS_matrix_NC_iR_Vnk!(Hx_uu, Hx_dd, Hx_ud, Hx_du, Hy_uu, Hy_dd, Hy_ud, Hy_du, Hz_uu, Hz_dd, Hz_ud, Hz_du, 
                             Hks, iHks, 
                             Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, 
                             kpts[k], cartesian_cell) 
        HS_matrix_iR_Vnk!(Sx, Sy, Sz, OLP, Natom, FNAN, natn, ncn, atv_ijk, Total_NumOrbs, MP, kpts[k], cartesian_cell) 

        @. Cnk_up_temp1 = @views(Cnk[1][k][1:fsize,:]')
        @. Cnk_dn_temp1 = @views(Cnk[1][k][fsize+1:Nfsize,:]')
        @. Cnk_up_temp2 = @views(Cnk[1][k][1:fsize,:])
        @. Cnk_dn_temp2 = @views(Cnk[1][k][fsize+1:Nfsize,:])
        
        # up-up part
        mul!(Htemp1, Cnk_up_temp1, Hx_uu)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][1] = real(diag(Htemp2))
        mul!(Htemp1, Cnk_up_temp1, Hy_uu)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][2] = real(diag(Htemp2))
        mul!(Htemp1, Cnk_up_temp1, Hz_uu)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][3] = real(diag(Htemp2))

        mul!(Htemp1, Cnk_up_temp1, Sx)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][1][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end
        mul!(Htemp1, Cnk_up_temp1, Sy)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][2][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end
        mul!(Htemp1, Cnk_up_temp1, Sz)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][3][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end

        # dn-dn part
        mul!(Htemp1, Cnk_dn_temp1, Hx_dd)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][1] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_dn_temp1, Hy_dd)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][2] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_dn_temp1, Hz_dd)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][3] += real(diag(Htemp2))

        mul!(Htemp1, Cnk_dn_temp1, Sx)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][1][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end
        mul!(Htemp1, Cnk_dn_temp1, Sy)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][2][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end
        mul!(Htemp1, Cnk_dn_temp1, Sz)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        @inbounds for μ = 1:Nfsize
            Vnk[1][ik][jk][kk][3][μ] -= real(Htemp2[μ,μ])*Enk[1][ik][jk][kk][μ]
        end

        # up-dn part
        mul!(Htemp1, Cnk_up_temp1, Hx_ud)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][1] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_up_temp1, Hy_ud)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][2] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_up_temp1, Hz_ud)
        mul!(Htemp2, Htemp1, Cnk_dn_temp2)
        Vnk[1][ik][jk][kk][3] += real(diag(Htemp2))

        # dn-up part
        mul!(Htemp1, Cnk_dn_temp1, Hx_du)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][1] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_dn_temp1, Hy_du)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][2] += real(diag(Htemp2))
        mul!(Htemp1, Cnk_dn_temp1, Hz_du)
        mul!(Htemp2, Htemp1, Cnk_up_temp2)
        Vnk[1][ik][jk][kk][3] += real(diag(Htemp2))
    end
end
