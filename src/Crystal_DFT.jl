@timeit timer "Crystal_DFT" function Crystal_DFT!(
    Ham::Hamiltonian, 
    system_grid::System_Grid, 
    electron::CrystalBloch, 
    kpoints::KPoints,
    Hks, DM)

    SpinPol = Ham.SpinPol

    if SpinPol == "off"
        Crystal_DFT_Collinear_nonpol!(Ham.OLP, Hks[1], electron, kpoints, system_grid)
        Calc_Band_Energy!(electron, kpoints)
        Calc_DM_Crystal_Collinear_nopol!(electron, kpoints, system_grid, DM)
    elseif SpinPol == "on"
        Crystal_DFT_Collinear_pol!(Ham.OLP, Hks, electron, kpoints, system_grid)
        Calc_Band_Energy!(electron, kpoints)
        Calc_DM_Crystal_Collinear_pol!(electron, kpoints, system_grid, DM)
    elseif SpinPol == "nc"
        Crystal_DFT_NonCollinear!(Ham.OLP, Hks, Ham.iHNL, electron, kpoints, system_grid)   
        Calc_Band_Energy!(electron, kpoints)
        Calc_DM_Crystal_NonCollinear!(electron, kpoints, system_grid, DM)
    else
        println("SpinPol = ", SpinPol)
        error("please check SpinPol")
    end
end


@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Collinear_nonpol!( 
    OLP, Hks, 
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)

    S = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        H = Cnk[1][ik]
        HS_matrix!(S, H, OLP, Hks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        decomposition = eigen!(Hermitian(H), Hermitian(S))
        @views Enk[:,knum+ik,1] .= decomposition.values
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Collinear_pol!( 
    OLP, Hks, 
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    

    S_base = zeros(ComplexF64, fsize, fsize)
    S_work = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        H_up = Cnk[1][ik]
        HS_matrix!(S_base, H_up, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        copyto!(S_work, S_base)
        decomposition_up = eigen!(Hermitian(H_up), Hermitian(S_work))
        @views Enk[:,knum+ik,1] .= decomposition_up.values

        H_dn = Cnk[2][ik]
        HS_matrix!(H_dn, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        copyto!(S_work, S_base)
        decomposition_dn = eigen!(Hermitian(H_dn), Hermitian(S_work))
        @views Enk[:,knum+ik,2] .= decomposition_dn.values
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


@timeit timer "Crystal_DFT_NonCollinear" function Crystal_DFT_NonCollinear!( 
    OLP, Hks, iHNL,
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    fsize = sum(Total_NumOrbs)
    
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)

    tmpH = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, 2*fsize, 2*fsize)
    
    for ik = 1:MPI_Nkpt
        H = Cnk[1][ik]
        HS_matrix_NC!(H, Hks, iHNL, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        fill!(S, 0.0)
        @. @views S[1:fsize, 1:fsize] = tmpH
        @. @views S[fsize+1:end, fsize+1:end] = tmpH
        decomposition = eigen!(Hermitian(H), Hermitian(S))
        @views Enk[:,knum+ik,1] .= decomposition.values
    end
    MPI.Allreduce!(Enk, MPI.SUM, comm)
    electron.Enk = Enk
    electron.Cnk = Cnk
end


# pre calculation of fermi function
function Calc_fnkCnk!(electron::CrystalBloch, kpoints::KPoints)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    SpinPol = electron.SpinPol
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk
    FF = electron.FF

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kweight = kpoints.MPI_kweight
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]
    
    if SpinPol == "off"
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            @views rmul!(Cnk[1][ik][:,μ], F)
        end
    elseif SpinPol == "on"
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F_up = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            F_dn = sqrt(FF[μ,ik+knum,2]*MPI_kweight[ik])
            @views rmul!(Cnk[1][ik][:,μ], F_up)
            @views rmul!(Cnk[2][ik][:,μ], F_dn)
        end
    elseif SpinPol == "nc"
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            @views rmul!(Cnk[1][ik][:,μ], F)
        end
    else
        error("please check SpinPol")
    end

    electron.Cnk = Cnk
end


@timeit timer "Calc_DM_Crystal_Collinear" function Calc_DM_Crystal_Collinear_nopol!(
    electron::CrystalBloch, 
    kpoints::KPoints, 
    system_grid::System_Grid, 
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    Calc_fnkCnk!(electron, kpoints)

    fill!(DM[1], 0.0)
    DM_dense = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        mul!(DM_dense, ctemp, adjoint(ctemp))
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            @inbounds for ist = 1:NO0, jst = 1:NO1
                DMst += 1
                matrix_element = DM_dense[Bnum+jst,Anum+ist]
                DM[1][DMst] += real(matrix_element*ex)
            end
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Crystal_Collinear" function Calc_DM_Crystal_Collinear_pol!(
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid, 
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    Calc_fnkCnk!(electron, kpoints)


    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)
    DM_dense_up = zeros(ComplexF64, fsize, fsize)
    DM_dense_dn = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp_up = Cnk[1][ik]
        ctemp_dn = Cnk[2][ik]
        mul!(DM_dense_up, ctemp_up, adjoint(ctemp_up))
        mul!(DM_dense_dn, ctemp_dn, adjoint(ctemp_dn))
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            @inbounds for ist = 1:NO0, jst = 1:NO1
                DMst += 1
                matrix_up = DM_dense_up[Bnum+jst,Anum+ist]
                matrix_dn = DM_dense_dn[Bnum+jst,Anum+ist]
                DM[1][DMst] += real(matrix_up*ex)
                DM[2][DMst] += real(matrix_dn*ex)
            end
        end
    end


    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
end


@timeit timer "Calc_DM_Crystal_NonCollinear" function Calc_DM_Crystal_NonCollinear!(
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid,
    DM::Vector{Vector{Float64}})
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    Calc_fnkCnk!(electron, kpoints)

    fill!(DM[1], 0.0)
    fill!(DM[2], 0.0)
    fill!(DM[3], 0.0)
    fill!(DM[4], 0.0)
    DM_dense = zeros(ComplexF64, Nfsize, Nfsize)
    

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        mul!(DM_dense, ctemp, adjoint(ctemp))
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            @inbounds for ist = 1:NO0, jst = 1:NO1
                DMst += 1
                matrix_uu = DM_dense[Bnum+jst,Anum+ist]
                matrix_dd = DM_dense[fsize+Bnum+jst,fsize+Anum+ist]
                matrix_ud = DM_dense[fsize+Bnum+jst,Anum+ist]
                DM[1][DMst] += real(matrix_uu*ex)
                DM[2][DMst] += real(matrix_dd*ex)
                DM[3][DMst] += real(matrix_ud*ex)
                DM[4][DMst] += imag(matrix_ud*ex)
            end
        end
    end

    MPI.Allreduce!(DM[1], MPI.SUM, comm)
    MPI.Allreduce!(DM[2], MPI.SUM, comm)
    MPI.Allreduce!(DM[3], MPI.SUM, comm)
    MPI.Allreduce!(DM[4], MPI.SUM, comm)
end


@timeit timer "Calc_iDM_Crystal_NonCollinear" function Calc_iDM_Crystal_NonCollinear!(
    electron::CrystalBloch, 
    kpoints::KPoints,
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
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts


    iDM_dense = zeros(ComplexF64, Nfsize, Nfsize)

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        fill!(iDM[1][atom][Rn][ist], 0.0)
        fill!(iDM[2][atom][Rn][ist], 0.0)
    end

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        mul!(iDM_dense, ctemp, adjoint(ctemp))
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            iDM1 = iDM[1][atom][Rn]
            iDM2 = iDM[2][atom][Rn]
            @inbounds for ist = 1:NO0, jst = 1:NO1
                matrix_uu = iDM_dense[Bnum+jst,Anum+ist]
                matrix_dd = iDM_dense[fsize+Bnum+jst,fsize+Anum+ist]
                iDM1[ist][jst] += imag(matrix_uu*ex)
                iDM2[ist][jst] += imag(matrix_dd*ex)
            end
        end
    end

    @inbounds for atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom]
        MPI.Allreduce!(iDM[1][atom][Rn][ist], MPI.SUM, comm)
        MPI.Allreduce!(iDM[2][atom][Rn][ist], MPI.SUM, comm)
    end
end


# calculation Energy Density matrix when using Force calculation
function Calc_EDM(electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid, occ_flag::Bool)

    SpinPol = electron.SpinPol
    Total_Hsize = system_grid.Total_Hsize
    Nspin = ifelse(SpinPol=="off", 1, 2)

    EDM = Vector{Vector{Float64}}(undef, Nspin)
    for spin = 1:Nspin
        EDM[spin] = zeros(Float64, Total_Hsize)
    end

    if occ_flag
        Calc_fnkCnk!(electron, kpoints)
    end

    if SpinPol == "off"
        Calc_EDM_Collinear_nonpol!(EDM, electron, kpoints, system_grid)
    elseif SpinPol == "on"
        Calc_EDM_Collinear_pol!(EDM, electron, kpoints, system_grid)
    elseif SpinPol == "nc"
        Calc_EDM_NonCollinear!(EDM, electron, kpoints, system_grid)
    else
        error("please check")
    end


    return EDM
end


function Calc_EDM_Collinear_nonpol!(EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]
    weighted_C = zeros(ComplexF64, fsize, fsize)
    EDM_dense = zeros(ComplexF64, fsize, fsize)

    
    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        @views weighted_C .= ctemp .* transpose(Enk[:,ik+knum,1])
        mul!(EDM_dense, weighted_C, adjoint(ctemp))
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt
            
            @inbounds for ist = 1:NO0, jst = 1:NO1
                hst += 1
                matrix_element = EDM_dense[Bnum+jst,Anum+ist]
                EDM[1][hst] += real(matrix_element*ex)
            end
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
end


function Calc_EDM_Collinear_pol!(EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]
    weighted_C_up = zeros(ComplexF64, fsize, fsize)
    weighted_C_dn = zeros(ComplexF64, fsize, fsize)
    EDM_dense_up = zeros(ComplexF64, fsize, fsize)
    EDM_dense_dn = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp_up = Cnk[1][ik]
        ctemp_dn = Cnk[2][ik]
        @views weighted_C_up .= ctemp_up .* transpose(Enk[:,ik+knum,1])
        @views weighted_C_dn .= ctemp_dn .* transpose(Enk[:,ik+knum,2])
        mul!(EDM_dense_up, weighted_C_up, adjoint(ctemp_up))
        mul!(EDM_dense_dn, weighted_C_dn, adjoint(ctemp_dn))
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt
            
            @inbounds for ist = 1:NO0, jst = 1:NO1
                hst += 1
                matrix_up = EDM_dense_up[Bnum+jst,Anum+ist]
                matrix_dn = EDM_dense_dn[Bnum+jst,Anum+ist]
                EDM[1][hst] += real(matrix_up*ex)
                EDM[2][hst] += real(matrix_dn*ex)
            end
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
    MPI.Allreduce!(EDM[2], MPI.SUM, comm)
end


function Calc_EDM_NonCollinear!(EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    fsize = electron.fsize
    Nfsize = electron.Nfsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]
    weighted_C = zeros(ComplexF64, Nfsize, Nfsize)
    EDM_dense = zeros(ComplexF64, Nfsize, Nfsize)
    
    
    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        @views weighted_C .= ctemp .* transpose(Enk[:,ik+knum,1])
        mul!(EDM_dense, weighted_C, adjoint(ctemp))
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
            ex = cispi(2*kRn)/AllNkpt

            @inbounds for ist = 1:NO0, jst = 1:NO1
                hst += 1
                matrix_uu = EDM_dense[Bnum+jst,Anum+ist]
                matrix_dd = EDM_dense[fsize+Bnum+jst,fsize+Anum+ist]
                EDM[1][hst] += real(matrix_uu*ex)
                EDM[2][hst] += real(matrix_dd*ex)
            end
        end
    end

    MPI.Allreduce!(EDM[1], MPI.SUM, comm)
    MPI.Allreduce!(EDM[2], MPI.SUM, comm)
end
