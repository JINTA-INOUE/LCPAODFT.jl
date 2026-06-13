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


function HS_matrix!(
    S, H, 
    OLP::Vector{Vector{Vector{Vector{Float64}}}}, 
    Hks::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
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


function HS_matrix!(
    HS, 
    A::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
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


function HS_matrix!(
    HS, 
    A::Vector{Vector{Vector{Vector{Float64}}}}, 
    B::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
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


function HS_matrix!(
    HS, 
    A::Vector{Vector{Vector{Vector{Float64}}}}, 
    B::Vector{Vector{Vector{Vector{Float64}}}}, 
    C::Vector{Vector{Vector{Vector{Float64}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
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


function HS_matrix!(
    S, H, 
    OLP::Vector{Float64}, 
    Hks::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    fill!(S, 0.0)
    fill!(H, 0.0)
    
    hst = 0
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


function HS_matrix!(
    HS, 
    A::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    fill!(HS, 0.0)
    
    hst = 0
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


function HS_matrix!(
    HS, 
    A::Vector{Float64}, 
    B::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    fill!(HS, 0.0)

    hst = 0
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


function HS_matrix!(
    HS, 
    A::Vector{Float64}, 
    B::Vector{Float64}, 
    C::Vector{Float64}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})

    fill!(HS, 0.0)

    hst = 0
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


function HS_matrix_NC!(
    tmpH, H, 
    Hks::Union{Vector{Vector{Float64}},Vector{Vector{Vector{Vector{Vector{Float64}}}}}}, 
    iHks::Union{Vector{Vector{Float64}},Vector{Vector{Vector{Vector{Vector{Float64}}}}}}, 
    Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, 
    kpts::Vector{Float64})
    
    fsize = sum(Total_NumOrbs)
    
    HS_matrix!(tmpH, Hks[1], iHks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[1:fsize,1:fsize] = tmpH

    HS_matrix!(tmpH, Hks[2], iHks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[fsize+1:end,fsize+1:end] = tmpH

    HS_matrix!(tmpH, Hks[3], Hks[4], iHks[3], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts)
    @. H[1:fsize,fsize+1:end] = tmpH
end


@timeit timer "Crystal_DFT_Collinear" function Crystal_DFT_Collinear_nonpol!( 
    OLP, Hks, 
    electron::CrystalBloch, 
    kpoints::KPoints,
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

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    fsize = sum(Total_NumOrbs)

    S = zeros(ComplexF64, fsize, fsize)
    H = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        HS_matrix!(S, H, OLP, Hks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[:,knum+ik,1], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
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
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)
    fsize = sum(Total_NumOrbs)
    
    S = zeros(ComplexF64, fsize, fsize)
    H = zeros(ComplexF64, fsize, fsize)

    for ik = 1:MPI_Nkpt
        HS_matrix!(S, H, OLP, Hks[1], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[:,knum+ik,1], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
        HS_matrix!(H, Hks[2], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        Enk[:,knum+ik,2], Cnk[2][ik] = eigen(Hermitian(H), Hermitian(S))
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
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk
    
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    fill!(Enk, 0.0)

    tmpH = zeros(ComplexF64, fsize, fsize)
    S = zeros(ComplexF64, 2*fsize, 2*fsize)
    H = zeros(ComplexF64, 2*fsize, 2*fsize)
    
    for ik = 1:MPI_Nkpt
        HS_matrix_NC!(tmpH, H, Hks, iHNL, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, MPI_kpts[ik])
        @. @views S[1:fsize, 1:fsize] = tmpH
        @. @views S[fsize+1:end, fsize+1:end] = tmpH
        Enk[:,knum+ik,1], Cnk[1][ik] = eigen(Hermitian(H), Hermitian(S))
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
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            @views mul!(tmpCnk, F, Cnk[1][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk
        end
    elseif SpinPol == "on"
        tmpCnk_up = zeros(ComplexF64, Nfsize)
        tmpCnk_dn = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F_up = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            F_dn = sqrt(FF[μ,ik+knum,2]*MPI_kweight[ik])
            @views mul!(tmpCnk_up, F_up, Cnk[1][ik][:,μ])
            @views mul!(tmpCnk_dn, F_dn, Cnk[2][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk_up
            @. @views Cnk[2][ik][:,μ] = tmpCnk_dn
        end
    elseif SpinPol == "nc"
        tmpCnk = zeros(ComplexF64, Nfsize)
        @inbounds for ik = 1:MPI_Nkpt, μ = 1:Nfsize
            F = sqrt(FF[μ,ik+knum,1]*MPI_kweight[ik])
            @views mul!(tmpCnk, F, Cnk[1][ik][:,μ])
            @. @views Cnk[1][ik][:,μ] = tmpCnk
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
    DM)
    
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

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
        DMst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1

            Anum = MP[atom]
            NO0 = Total_NumOrbs[atom]
            cell = ncn[atom][Rn]+1
            jatom = natn[atom][Rn]
            Bnum = MP[jatom]
            NO1 = Total_NumOrbs[jatom]
            kRn = ka*atv_ijk[cell][1] + kb*atv_ijk[cell][2] + kc*atv_ijk[cell][3]
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


@timeit timer "Calc_DM_Crystal_Collinear" function Calc_DM_Crystal_Collinear_pol!(
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid, 
    DM)
    
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

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp_up = Cnk[1][ik]
        ctemp_dn = Cnk[2][ik]
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


@timeit timer "Calc_DM_Crystal_NonCollinear" function Calc_DM_Crystal_NonCollinear!(
    electron::CrystalBloch, 
    kpoints::KPoints,
    system_grid::System_Grid,
    DM)
    
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
    

    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
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
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts

    
    fill!(iDM[1], 0.0)
    fill!(iDM[2], 0.0)

    
    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
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

    if SpinPol ∈ ("off", "on")
        Calc_EDM_Collinear!(EDM, electron, kpoints, system_grid)
    elseif SpinPol == "nc"
        Calc_EDM_NonCollinear!(EDM, electron, kpoints, system_grid)
    else
        error("please check")
    end


    return EDM
end


function Calc_EDM_Collinear!(EDM, electron::CrystalBloch, kpoints::KPoints, system_grid::System_Grid)
    
    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)

    Natom = system_grid.Natom
    Total_NumOrbs = system_grid.Total_NumOrbs
    MP = system_grid.MP
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    ncn = system_grid.ncn
    atv_ijk = system_grid.atv_ijk

    spinsize = electron.spinsize
    fsize = electron.fsize
    Enk = electron.Enk
    Cnk = electron.Cnk

    AllNkpt = kpoints.AllNkpt
    MPI_Nkpt = kpoints.MPI_Nkpt
    MPI_kpts = kpoints.MPI_kpts
    MPkpts = kpoints.MPkpts
    knum = MPkpts[myrank+1]

    
    for spin = 1:spinsize, ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[spin][ik]
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
            
            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:fsize
                    cc1 += conj(ctemp[Anum+ist,μ])*ctemp[Bnum+jst,μ]*Enk[μ,ik+knum,spin]
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
    
    
    for ik = 1:MPI_Nkpt
        ka, kb, kc = MPI_kpts[ik]
        ctemp = Cnk[1][ik]
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

            for ist = 1:NO0, jst = 1:NO1
                cc1 = ComplexF64(0.0, 0.0)
                cc2 = ComplexF64(0.0, 0.0)
                @inbounds for μ = 1:Nfsize
                    cc1 += conj(ctemp[Anum+ist,μ]) * ctemp[Bnum+jst,μ] * Enk[μ,ik+knum,1]
                    cc2 += conj(ctemp[fsize+Anum+ist,μ]) * ctemp[fsize+Bnum+jst,μ] * Enk[μ,ik+knum,1]
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
