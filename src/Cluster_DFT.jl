@timeit timer "Cluster_DFT" function Cluster_DFT!(
    Ham::Hamiltonian, 
    system_grid::System_Grid, 
    electron::ClusterBloch,
    Hks, DM)

    SpinPol = Ham.SpinPol

    if SpinPol == "off"
        Cluster_DFT_Collinear_nonpol!(Ham.OLP, Hks[1], electron, system_grid)
        Calc_Band_Energy!(electron)
        Calc_DM_Cluster_Collinear_nopol!(electron, system_grid, DM)
    elseif SpinPol == "on"
        Cluster_DFT_Collinear_pol!(Ham.OLP, Hks, electron, system_grid)
        Calc_Band_Energy!(electron)
        Calc_DM_Cluster_Collinear_pol!(electron, system_grid, DM)
    elseif SpinPol == "nc"
        Cluster_DFT_NonCollinear!(Ham.OLP, Hks, Ham.iHNL, electron, system_grid)   
        Calc_Band_Energy!(electron)
        Calc_DM_Cluster_NonCollinear!(electron, system_grid, DM)
    else
        println("SpinPol = ", SpinPol)
        error("please check SpinPol")
    end
end


@timeit timer "Cluster_DFT_Collinear" function Cluster_DFT_Collinear_nonpol!( 
    OLP, Hks, 
    electron::ClusterBloch, 
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)

    S = zeros(Float64, fsize, fsize)
    H = zeros(Float64, fsize, fsize)

    HS_matrix!(S, H, OLP, Hks, Natom, Total_NumOrbs, MP, FNAN, natn)
    Enk[:,1], Cnk[1] = eigen(Symmetric(H), Symmetric(S))
    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Cluster_DFT_Collinear" function Cluster_DFT_Collinear_pol!( 
    OLP, Hks, 
    electron::ClusterBloch, 
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    
    S = zeros(Float64, fsize, fsize)
    H = zeros(Float64, fsize, fsize)

    HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn)
    Enk[:,1], Cnk[1] = eigen(Symmetric(H), Symmetric(S))
    HS_matrix!(H, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn)
    Enk[:,2], Cnk[2] = eigen(Symmetric(H), Symmetric(S))

    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Cluster_DFT_NonCollinear" function Cluster_DFT_NonCollinear!( 
    OLP, Hks, iHNL,
    electron::ClusterBloch,
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)
    
    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)

    tmpH = zeros(Float64, fsize, fsize)
    S = zeros(ComplexF64, 2*fsize, 2*fsize)
    H = zeros(ComplexF64, 2*fsize, 2*fsize)
    
    HS_matrix_NC!(H, Hks, iHNL, Natom, Total_NumOrbs, MP, FNAN, natn)
    HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn)
    @. @views S[1:fsize, 1:fsize] = tmpH
    @. @views S[fsize+1:end, fsize+1:end] = tmpH
    Enk[:,1], Cnk[1] = eigen(Hermitian(H), Hermitian(S))
    electron.Enk = Enk
    electron.Cnk = Cnk
end


# pre calculation of fermi function
function Calc_fnkCnk!(electron::ClusterBloch)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    SpinPol = electron.SpinPol
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk
    FF = electron.FF
    
    if SpinPol == "off"
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for μ = 1:Nfsize
            F = sqrt(FF[μ,1])
            @views mul!(tmpCnk, F, Cnk[1][:,μ])
            @. @views Cnk[1][:,μ] = tmpCnk
        end
    elseif SpinPol == "on"
        tmpCnk_up = zeros(ComplexF64, Nfsize)
        tmpCnk_dn = zeros(ComplexF64, Nfsize)
        @inbounds for μ = 1:Nfsize
            F_up = sqrt(FF[μ,1])
            F_dn = sqrt(FF[μ,2])
            @views mul!(tmpCnk_up, F_up, Cnk[1][:,μ])
            @views mul!(tmpCnk_dn, F_dn, Cnk[2][:,μ])
            @. @views Cnk[1][:,μ] = tmpCnk_up
            @. @views Cnk[2][:,μ] = tmpCnk_dn
        end
    elseif SpinPol == "nc"
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for μ = 1:Nfsize
            F = sqrt(FF[μ,1])
            @views mul!(tmpCnk, F, Cnk[1][:,μ])
            @. @views Cnk[1][:,μ] = tmpCnk
        end
    else
        error("please check SpinPol")
    end

    electron.Cnk = Cnk
end


@timeit timer "Calc_DM_Cluster_Collinear" function Calc_DM_Cluster_Collinear_nopol!(
    electron::ClusterBloch, 
    system_grid::System_Grid, 
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Cnk = electron.Cnk

    Calc_fnkCnk!(electron)

    fill!(DM[1], 0.0)

    ctemp = Cnk[1]
    DMst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1
        
        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
        
        for ist = 1:NO0, jst = 1:NO1
            tmp_re = 0.0
            @inbounds for μ = 1:fsize
                a = ctemp[Anum+ist,μ]
                b = ctemp[Bnum+jst,μ]
                tmp_re += real(a)*real(b) + imag(a)*imag(b)
            end
            DMst += 1
            DM[1][DMst] += tmp_re
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Cluster_Collinear" function Calc_DM_Cluster_Collinear_pol!(
    electron::ClusterBloch, 
    system_grid::System_Grid, 
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Cnk = electron.Cnk

    Calc_fnkCnk!(electron)


    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)

    ctemp_up = Cnk[1]
    ctemp_dn = Cnk[2]
    DMst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        for ist = 1:NO0, jst = 1:NO1
            tmp_up = ComplexF64(0.0, 0.0)
            tmp_dn = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:fsize
                tmp_up += conj(ctemp_up[Anum+ist,μ]) * ctemp_up[Bnum+jst,μ]
                tmp_dn += conj(ctemp_dn[Anum+ist,μ]) * ctemp_dn[Bnum+jst,μ]
            end
            DMst += 1
            DM[1][DMst] += real(tmp_up)
            DM[2][DMst] += real(tmp_dn)
        end
    end


    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Cluster_NonCollinear" function Calc_DM_Cluster_NonCollinear!(
    electron::ClusterBloch, 
    system_grid::System_Grid,
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk


    Calc_fnkCnk!(electron)

    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)
    fill!(DM[3], 0.0)
    fill!(DM[4], 0.0)
    

    ctemp = Cnk[1]
    DMst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

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
            DM[1][DMst] += real(cc1)
            DM[2][DMst] += real(cc2)
            DM[3][DMst] += real(cc3)
            DM[4][DMst] += imag(cc3)
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
    MPI.Allreduce!(DM[3], MPI.SUM, comm)
    MPI.Allreduce!(DM[4], MPI.SUM, comm)
end


@timeit timer "Calc_iDM_Cluster_NonCollinear" function Calc_iDM_Cluster_NonCollinear!(
    electron::ClusterBloch, 
    system_grid::System_Grid, 
    iDM)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    if isnothing(iDM)
        return
    end

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk
    
    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        fill!(iDM[1][atom][Rn][ist], 0.0)
        fill!(iDM[2][atom][Rn][ist], 0.0)
    end

    
    ctemp = Cnk[1]
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        iDM1 = iDM[1][atom][Rn]
        iDM2 = iDM[2][atom][Rn]
        for ist = 1:NO0, jst = 1:NO1
            cc1 = ComplexF64(0.0, 0.0)
            cc2 = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:Nfsize
                cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ]
                cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ]
            end
            iDM1[ist][jst] += imag(cc1)
            iDM2[ist][jst] += imag(cc2)
        end
    end


    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(iDM[1][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iDM[2][atom][Rn][ist], MPI.SUM, comm)
    end
end


# calculation Energy Density matrix when using Force calculation
function Calc_EDM(electron::ClusterBloch, kpoints::KPoints, system_grid::System_Grid, occ_flag::Bool)

    SpinPol = electron.SpinPol
    Total_Hsize = system_grid.Total_Hsize
    Nspin = ifelse(SpinPol=="off", 1, 2)

    EDM = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        EDM[spin] = zeros(Float64, Total_Hsize)
    end

    if occ_flag
        Calc_fnkCnk!(electron)
    end

    if SpinPol == "off"
        Calc_EDM_Collinear_nonpol!(EDM, electron, system_grid)
    elseif SpinPol == "on"
        Calc_EDM_Collinear_pol!(EDM, electron, system_grid)
    elseif SpinPol == "nc"
        Calc_EDM_NonCollinear!(EDM, electron, system_grid)
    else
        error("please check")
    end


    return EDM
end


function Calc_EDM_Collinear_nonpol!(EDM, electron::ClusterBloch, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    ctemp = Cnk[1]
    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
            
        for ist = 1:NO0, jst = 1:NO1
            cc1 = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:fsize
                cc1 += conj(ctemp[Anum+ist,μ])*ctemp[Bnum+jst,μ]*Enk[μ,1]
            end
            hst += 1
            EDM[1][hst] += real(cc1)
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
end


function Calc_EDM_Collinear_pol!(EDM, electron::ClusterBloch, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    ctemp_up = Cnk[1]
    ctemp_dn = Cnk[2]
    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]
            
        for ist = 1:NO0, jst = 1:NO1
            cc_up = ComplexF64(0.0, 0.0)
            cc_dn = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:fsize
                cc_up += conj(ctemp_up[Anum+ist,μ])*ctemp_up[Bnum+jst,μ]*Enk[μ,1]
                cc_dn += conj(ctemp_dn[Anum+ist,μ])*ctemp_dn[Bnum+jst,μ]*Enk[μ,2]
            end
            hst += 1
            EDM[1][hst] += real(cc_up)
            EDM[2][hst] += real(cc_dn)
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
    MPI.Allreduce!(EDM[2], MPI.SUM, comm)
end


function Calc_EDM_NonCollinear!(EDM, electron::ClusterBloch, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Enk = electron.Enk
    Cnk = electron.Cnk    
    
    ctemp = Cnk[1]
    hst = 0
    for atom = 1:Natom, Rn = 1:FNAN[atom]+1

        jatom = natn[atom][Rn]
        NO0 = Total_NumOrbs[atom]
        NO1 = Total_NumOrbs[jatom]
        Anum = MP[atom]
        Bnum = MP[jatom]

        for ist = 1:NO0, jst = 1:NO1
            cc1 = ComplexF64(0.0, 0.0)
            cc2 = ComplexF64(0.0, 0.0)
            @inbounds for μ = 1:Nfsize
                cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ] * Enk[μ,1]
                cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ] * Enk[μ,1]
            end
            hst += 1
            EDM[1][hst] += real(cc1)
            EDM[2][hst] += real(cc2)
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
    MPI.Allreduce!(EDM[2], MPI.SUM, comm)
end
