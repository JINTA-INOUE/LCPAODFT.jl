"""
    Crystal_DFT!

    1. Cluster_DFT_Collinear!
        solve eigen value problem Hc = ϵSc
    2. Calc_Band_Energy!
        calculation ε_n (electron.En) and fermi energy (ChemP)
    3. Calc_DM_Cluster_Collinear!
        calculation Density matrix (DM)

Mandatory arguments:

- `Ham`: an instance of `Hamiltonian`
- `electron`: an instance of `Electron`
- `Hks`: Hamiltonian matrix
- `DM`: Density matrix
"""
@timeit timer "Crystal_DFT" function Crystal_DFT!(Ham::Hamiltonian, system_grid::System_Grid, electron::CrystalBloch, Hks, DM)

    SpinPol = Ham.SpinPol

    if SpinPol == "off"

        Crystal_DFT_Collinear_nonpol!(Ham.OLP, Hks[1], electron, system_grid)
        Calc_Band_Energy!(electron)
        Calc_DM_Crystal_Collinear_nopol!(electron, system_grid, DM)

    elseif SpinPol == "on"

        Crystal_DFT_Collinear_pol!(Ham.OLP, Hks, electron, system_grid)
        Calc_Band_Energy!(electron )
        Calc_DM_Crystal_Collinear_pol!(electron, system_grid, DM)

    elseif SpinPol == "nc"

        Crystal_DFT_NonCollinear!(Ham.OLP, Hks, Ham.iHNL, electron, system_grid)   
        Calc_Band_Energy!(electron )
        Calc_DM_Crystal_NonCollinear!(electron, system_grid, DM)
        
    else
        println("SpinPol = ", SpinPol)
        error("please check SpinPol")
    end
end


function HS_matrix!(S, H, OLP::Vector{Vector{Vector{Vector{Float64}}}}, Hks::Vector{Vector{Vector{Vector{Float64}}}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    fill!(S, 0.0)
    fill!(H, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            S[Anum+ist,Bnum+jst] += OLP[atom][Rn][ist][jst]*ex
            H[Anum+ist,Bnum+jst] += Hks[atom][Rn][ist][jst]*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Vector{Vector{Vector{Float64}}}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            HS[Anum+ist,Bnum+jst] += A[atom][Rn][ist][jst]*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Vector{Vector{Vector{Float64}}}}, B::Vector{Vector{Vector{Vector{Float64}}}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            HS[Anum+ist,Bnum+jst] += (A[atom][Rn][ist][jst]+im*B[atom][Rn][ist][jst])*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Vector{Vector{Vector{Float64}}}}, B::Vector{Vector{Vector{Vector{Float64}}}}, C::Vector{Vector{Vector{Vector{Float64}}}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            HS[Anum+ist,Bnum+jst] += (A[atom][Rn][ist][jst]+im*(B[atom][Rn][ist][jst]+C[atom][Rn][ist][jst]))*ex
        end
    end
end


function HS_matrix_NC!(tmpH, H, Hks::Vector{Vector{Vector{Vector{Vector{Float64}}}}}, iHks::Vector{Vector{Vector{Vector{Vector{Float64}}}}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    fsize = sum(Total_NumOrbs)
    
    HS_matrix!(tmpH, Hks[1], iHks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[1:fsize,1:fsize] = tmpH

    HS_matrix!(tmpH, Hks[2], iHks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[fsize+1:end,fsize+1:end] = tmpH

    HS_matrix!(tmpH, Hks[3], Hks[4], iHks[3], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[1:fsize,fsize+1:end] = tmpH
end


function HS_matrix!(S, H, OLP::Vector{Float64}, Hks::Vector{Float64}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    hst = 0
    fill!(S, 0.0)
    fill!(H, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            S[Anum+ist,Bnum+jst] += OLP[hst]*ex
            H[Anum+ist,Bnum+jst] += Hks[hst]*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Float64}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    hst = 0
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            HS[Anum+ist,Bnum+jst] += A[hst]*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Float64}, B::Vector{Float64}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    hst = 0
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            HS[Anum+ist,Bnum+jst] += (A[hst]+im*B[hst])*ex
        end
    end
end


function HS_matrix!(HS, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    hst = 0
    fill!(HS, 0.0)
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        NO0 = Total_NumOrbs[atom]
        Anum = MP[atom]
        jatom = natn[atom][Rn]
        cell = ncn[atom][Rn]+1
        NO1 = Total_NumOrbs[jatom]
        kRn = dot(kpts, atv_ijk[cell])
        Bnum = MP[jatom]
        ex = cispi(2*kRn)
        @inbounds for ist = 1:NO0, jst = 1:NO1
            hst += 1
            HS[Anum+ist,Bnum+jst] += (A[hst]+im*(B[hst]+C[hst]))*ex
        end
    end
end


function HS_matrix_NC!(tmpH, H, Hks::Vector{Vector{Float64}}, iHNL::Vector{Vector{Float64}}, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts::Vector{Float64})
    fsize = sum(Total_NumOrbs)
    
    HS_matrix!(tmpH, Hks[1], iHNL[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. @views H[1:fsize,1:fsize] = tmpH

    HS_matrix!(tmpH, Hks[2], iHNL[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. @views H[fsize+1:end,fsize+1:end] = tmpH

    HS_matrix!(tmpH, Hks[3], Hks[4], iHNL[3], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. @views H[1:fsize,fsize+1:end] = tmpH
end


@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Collinear_nonpol!( 
    OLP, Hks, 
    electron::CrystalBloch, 
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    fsize = sum(Total_NumOrbs)
    
    S = zeros(ComplexF64, fsize, fsize)
    H = zeros(ComplexF64, fsize, fsize)

    @inbounds for ik = 1:MPI_Nkpt
        HS_matrix!(S, H, OLP, Hks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[1,:,knum+ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Collinear_pol!( 
    OLP, Hks, 
    electron::CrystalBloch, 
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    fsize = sum(Total_NumOrbs)
    
    S = zeros(ComplexF64, fsize, fsize)
    H = zeros(ComplexF64, fsize, fsize)

    @inbounds for ik = 1:MPI_Nkpt
        HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[1,:,knum+ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
        HS_matrix!(H, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[2,:,knum+ik], Cnk[2][ik] = eigen(Hermitian(H), Hermitian(S))
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Crystal_DFT_NonCollinear" function Crystal_DFT_NonCollinear!( 
    OLP, Hks, iHNL,
    electron::CrystalBloch, 
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    fsize = sum(Total_NumOrbs)

    tmpH = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, 2*fsize, 2*fsize)
    H = zeros(ComplexF64, 2*fsize, 2*fsize)
    
    @inbounds for ik = 1:MPI_Nkpt
        HS_matrix_NC!(tmpH, H, Hks, iHNL, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        @. @views S[1:fsize, 1:fsize] = tmpH
        @. @views S[fsize+1:end, fsize+1:end] = tmpH
        Enk[1,:,knum+ik], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


function Calc_ChemP(Spindeg, TZ, kweight, Beta, Enk; loopmax=2000, max_x=60.0)
    
    spinsize, Nfsize, Nkpt = size(Enk)
    All_Nkpt = sum(kweight)
     
    ChemP = 0.0
    ChemP_max = 20.0
    ChemP_min = -20.0
    loopN = 0

    while loopN < loopmax
        
        loopN += 1
        ChemP = 0.5*(ChemP_min + ChemP_max)
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:Nfsize
                
            x = (Enk[spin,μ,ik] - ChemP)*Beta*0.2
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF*kweight[ik]
        end

        Num_state = Spindeg*Num_state/All_Nkpt
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-14
            break
        end
    end


    ChemP_max = 20.0
    ChemP_min = -20.0
    loopN = 0
    
    while loopN < loopmax
        
        loopN += 1
        if loopN ≠ 1
            ChemP = 0.5*(ChemP_min + ChemP_max)
        end
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:Nfsize
                
            x = (Enk[spin,μ,ik] - ChemP)*Beta
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF*kweight[ik]
        end

        Num_state = Spindeg*Num_state/All_Nkpt
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-14
            break
        end
    end

    
    return ChemP
end


function Calc_Band_Energy!(electron::CrystalBloch)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Spindeg = electron.Spindeg
    TotalZ = electron.TotalZ
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    Enk = electron.Enk
    FF = electron.FF
    AllNkpt = electron.kpoints.AllNkpt
    All_kweight = electron.kpoints.All_kweight

    spinsize, Nfsize, Nkpt = size(Enk)

    ChemP = Calc_ChemP(Spindeg, TotalZ, All_kweight, Beta, Enk)
    electron.ChemP = ChemP


    Eele0 = 0.0
    @inbounds for ik = 1:Nkpt, μ = 1:Nfsize, spin = 1:spinsize 
            
        x = (Enk[spin,μ,ik]-ChemP)*Beta
        if x <= -max_x
            x = -max_x
         end

        if x >= max_x
            x = max_x
        end
            
        FermiF = 1/(1 + exp(x))
        FF[spin,μ,ik] = FermiF
        Eele0 += FermiF*Enk[spin,μ,ik]*All_kweight[ik]
    end

    electron.FF = FF
    electron.Eele = Spindeg*Eele0/AllNkpt
end


# pre calculation of fermi function
function Calc_fnkCnk!(electron::CrystalBloch)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    SpinPol = electron.SpinPol
    Cnk = electron.Cnk
    FF = electron.FF
    Nfsize = size(FF,2)

    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kweight = electron.kpoints.MPI_kweight
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]
    
    if SpinPol == "off"
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[1,μ,ik+knum]*MPI_kweight[ik])
            @views mul!(tmpCnk, F, Cnk[1][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk
        end
    elseif SpinPol == "on"
        tmpCnk_up = zeros(ComplexF64, Nfsize)
        tmpCnk_dn = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F_up = sqrt(FF[1,μ,ik+knum]*MPI_kweight[ik])
            F_dn = sqrt(FF[2,μ,ik+knum]*MPI_kweight[ik])
            @views mul!(tmpCnk_up, F_up, Cnk[1][ik][:,μ])
            @views mul!(tmpCnk_dn, F_dn, Cnk[2][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk_up
            @. @views Cnk[2][ik][:,μ] = tmpCnk_dn
        end
    elseif SpinPol == "nc"
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[1,μ,ik+knum]*MPI_kweight[ik])
            @views mul!(tmpCnk, F, Cnk[1][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk
        end
    else
        error("please check SpinPol")
    end

    electron.Cnk = Cnk
end


@timeit timer "Calc_DM_Crystal_Collinear" function Calc_DM_Crystal_Collinear_nopol!(electron::CrystalBloch, system_grid::System_Grid, DM)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    Cnk = electron.Cnk
    fsize = sum(Total_NumOrbs)

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts

    Calc_fnkCnk!(electron)

    ctemp = zeros(ComplexF64, fsize, fsize)
    fill!(DM[1], 0.0)

    for ik = 1:MPI_Nkpt
        @. ctemp = Cnk[1][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            coskRn = cos(2*pi*kRn)/AllNkpt
            sinkRn = sin(2*pi*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                tmp_re = 0.0
                tmp_im = 0.0
                @inbounds for μ = 1:fsize
                    a = ctemp[Anum+ist,μ]
                    b = ctemp[Bnum+jst,μ]
                    tmp_re += real(a)*real(b) + imag(a)*imag(b)
                    tmp_im += imag(a)*real(b) - real(a)*imag(b)
                end
                DMst += 1
                DM[1][DMst] += tmp_re*coskRn + tmp_im*sinkRn
            end
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Crystal_Collinear" function Calc_DM_Crystal_Collinear_pol!(electron::CrystalBloch, system_grid::System_Grid, DM)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    Cnk = electron.Cnk
    fsize = sum(Total_NumOrbs)

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts

    Calc_fnkCnk!(electron)


    ctemp_up = zeros(ComplexF64, fsize, fsize)
    ctemp_dn = zeros(ComplexF64, fsize, fsize)

    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)

    for ik = 1:MPI_Nkpt
        @. ctemp_up = Cnk[1][ik]
        @. ctemp_dn = Cnk[2][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                tmp_up = ComplexF64(0.0, 0.0)
                tmp_dn = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:fsize
                    tmp_up += conj(ctemp_up[Anum+ist,μ]) * ctemp_up[Bnum+jst,μ]
                    tmp_dn += conj(ctemp_dn[Anum+ist,μ]) * ctemp_dn[Bnum+jst,μ]
                end
                DMst += 1
                DM[1][DMst] += real(tmp_up*ex)
                DM[2][DMst] += real(tmp_dn*ex)
            end
        end
    end


    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Crystal_NonCollinear" function Calc_DM_Crystal_NonCollinear!(electron::CrystalBloch, system_grid::System_Grid, DM)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    Cnk = electron.Cnk
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts


    Calc_fnkCnk!(electron)

    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)
    fill!(DM[3], 0.0)
    fill!(DM[4], 0.0)
    

    ctemp = zeros(ComplexF64, Nfsize, Nfsize)

    for ik = 1:MPI_Nkpt
        @. ctemp = Cnk[1][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                cc2 = ComplexF64(0.0, 0.0)
                cc3 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:Nfsize
                    cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ]
                    cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ]
                    cc3 += conj(ctemp[Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ]
                end
                DMst += 1
                DM[1][DMst] += real(cc1*ex)
                DM[2][DMst] += real(cc2*ex)
                DM[3][DMst] += real(cc3*ex)
                DM[4][DMst] += imag(cc3*ex)
            end
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
    MPI.Allreduce!(DM[3], MPI.SUM, comm)
    MPI.Allreduce!(DM[4], MPI.SUM, comm)
end


@timeit timer "Calc_iDM_Crystal_NonCollinear" function Calc_iDM_Crystal_NonCollinear!(electron::CrystalBloch, system_grid::System_Grid, iDM)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    Cnk = electron.Cnk
    fsize = sum(Total_NumOrbs)

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts

    
    fill!(iDM[1], 0.0)
    fill!(iDM[2], 0.0)

    
    ctemp = zeros(ComplexF64, 2*fsize, 2*fsize)
    for ik = 1:MPI_Nkpt
        @. ctemp = Cnk[1][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                cc2 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:2*fsize
                    cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ]
                    cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ]
                end
                DMst += 1
                iDM[1][DMst] += imag(cc1*ex)
                iDM[2][DMst] += imag(cc2*ex)
            end
        end
    end


    MPI.Allreduce!(iDM[1], MPI.SUM, comm)
    MPI.Allreduce!(iDM[2], MPI.SUM, comm)
end


# calculation Energy Density matrix when using Force calculation
function Calc_EDM(electron::CrystalBloch, system_grid::System_Grid, occ_flag::Bool)

    SpinPol = electron.SpinPol
    Total_Hsize = system_grid.Total_Hsize
    spinsize = ifelse(SpinPol=="off", 1, 2)

    EDM = Vector{Vector{Float64}}(undef, spinsize)
    for spin = 1:spinsize
        EDM[spin] = zeros(Float64, Total_Hsize)
    end

    if occ_flag
        Calc_fnkCnk!(electron)
    end

    if SpinPol ∈ ("off", "on")
        Calc_EDM_Collinear!(EDM, electron, system_grid)
    elseif SpinPol == "nc"
        Calc_EDM_NonCollinear!(EDM, electron, system_grid)
    else
        error("please check")
    end


    return EDM
end


function Calc_EDM_Collinear!(EDM, electron::CrystalBloch, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    SpinPol = electron.SpinPol
    spinsize = ifelse(SpinPol=="off", 1, 2)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    fsize = sum(Total_NumOrbs)

    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]

    
    ctemp = zeros(ComplexF64, fsize, fsize)
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        @. ctemp = Cnk[spin][ik]
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt
            
            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:fsize
                    cc1 += conj(ctemp[Anum+ist,μ])*ctemp[Bnum+jst,μ]*Enk[spin,μ,ik+knum]
                end
                hst += 1
                EDM[spin][hst] += real(cc1*ex)
            end
        end
    end

    for spin = 1:spinsize
        MPI.Allreduce!(EDM[spin], MPI.SUM, comm)
    end
end


function Calc_EDM_NonCollinear!(EDM, electron::CrystalBloch, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize

    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = electron.kpoints.AllNkpt
    MPI_Nkpt = electron.kpoints.MPI_Nkpt
    MPI_kpts = electron.kpoints.MPI_kpts
    MPkpts = electron.kpoints.MPkpts
    knum = MPkpts[myrank+1]
    
    
    ctemp = zeros(ComplexF64, Nfsize, Nfsize)
    for ik = 1:MPI_Nkpt
        @. ctemp = Cnk[1][ik]
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = MPI_kpts[ik][1]*atv_ijk[cell][1] + MPI_kpts[ik][2]*atv_ijk[cell][2] + MPI_kpts[ik][3]*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                cc2 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:Nfsize
                    cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ] * Enk[1,μ,ik+knum]
                    cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ] * Enk[1,μ,ik+knum]
                end
                hst += 1
                EDM[1][hst] += real(cc1*ex)
                EDM[2][hst] += real(cc2*ex)
            end
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
    MPI.Allreduce!(EDM[2], MPI.SUM, comm)
end
