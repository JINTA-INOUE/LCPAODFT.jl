function Calc_EVec_Col!(spin, ik, Natom, Total_NumOrbs, MP, S, SD, Cnk, EVec)

    fsize = sum(Total_NumOrbs)
    for μ = 1:fsize
        fill!(SD, 0.0f0)
        for atom = 1:Natom, jatom = 1:Natom
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]
            Anum = MP[atom]
            Bnum = MP[jatom]

            for ist = 1:NO0, jst = 1:NO1
                tmp = conj(Cnk[Anum+ist,μ])*Cnk[Bnum+jst,μ]
                SD[Anum+ist] += real(tmp*S[Anum+ist,Bnum+jst])
            end
        end

        for ist = 1:fsize
            EVec[μ,ist,spin,ik] = SD[ist]
        end
    end
end


function Calc_EVec_NonCol!(ik, Natom, Total_NumOrbs, MP, Angle_spin, S, SD, Cnk, EVec)


    fsize = sum(Total_NumOrbs)
    Nfsize = 2*fsize
    for μ = 1:Nfsize
        fill!(SD, 0.0f0)
        for atom = 1:Natom, jatom = 1:Natom
            sit = sin(Angle_spin[atom][1])
            cot = cos(Angle_spin[atom][1])
            sip = sin(Angle_spin[atom][2])
            cop = cos(Angle_spin[atom][2])
            NO0 = Total_NumOrbs[atom]
            NO1 = Total_NumOrbs[jatom]
            Anum = MP[atom]
            Bnum = MP[jatom]

            for ist = 1:NO0, jst = 1:NO1
                svalue = S[Anum+ist,Bnum+jst]
                tmp_uu = real(conj(Cnk[Anum+ist,μ])*Cnk[Bnum+jst,μ]*svalue)
                tmp_dd = real(conj(Cnk[Anum+ist+fsize,μ])*Cnk[Bnum+jst+fsize,μ]*svalue)
                tmp_ud = conj(Cnk[Anum+ist,μ]) * Cnk[Bnum+jst+fsize,μ]*svalue
                tmp_ud_real = real(tmp_ud)
                tmp_ud_imag = imag(tmp_ud)

                SD[Anum+ist] += 0.5*(tmp_uu + tmp_dd) + 0.5*cot*(tmp_uu - tmp_dd) + (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
                SD[Anum+ist+fsize] += 0.5*(tmp_uu + tmp_dd) - 0.5*cot*(tmp_uu - tmp_dd) - (tmp_ud_real*cop - tmp_ud_imag*sip)*sit
            end
        end

        for ist = 1:fsize
            EVec[μ,ist,1,ik] = SD[ist]
            EVec[μ,ist,2,ik] = SD[ist + fsize]
        end
    end
end


function Calc_EVec_CWF!(SpinPol, spin, ik, Nwann, Cnk, EVec)

    if SpinPol in ("off", "on")
        @inbounds for ist = 1:Nwann, μ = 1:Nwann
            EVec[μ,ist,spin,ik] = abs2(Cnk[ist,μ])
        end
    elseif SpinPol == "nc"
        Nwann_half = div(Nwann, 2)
        @inbounds for ist = 1:Nwann_half, μ = 1:Nwann
            EVec[μ,ist,1,ik] = abs2(Cnk[ist,μ])
            EVec[μ,ist,2,ik] = abs2(Cnk[ist+Nwann_half,μ])
        end
    end
end


function _setup_band_kpath(kpath::Vector{Vector{Float64}}, Nk::Integer)

    length(kpath) >= 2 || error("kpath must contain at least two points")
    all(length(kpt) == 3 for kpt in kpath) ||
        error("each k point must contain three coordinates")
    Nk >= 2 || error("Nk must be at least 2")

    Nkpath = length(kpath) - 1
    kpath_Nk = fill(Int(Nk), Nkpath)
    kpath_start = kpath[1:end-1]
    kpath_end = kpath[2:end]

    all_kpts = Vector{Vector{Float64}}(undef, sum(kpath_Nk))
    global_k = 0
    for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]
        global_k += 1
        fraction = (ipath - 1) / (kpath_Nk[ik] - 1)
        all_kpts[global_k] = [
            kpath_start[ik][axis] +
            (kpath_end[ik][axis] - kpath_start[ik][axis]) * fraction
            for axis = 1:3
        ]
    end

    return Nkpath, kpath_Nk, kpath_start, kpath_end, all_kpts
end


function _split_band_kpoints(all_kpts, comm)

    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    MPI_krange = split_evenly(1:length(all_kpts), nprocs)
    local_range = MPI_krange[myrank + 1]

    return MPI_krange, all_kpts[local_range]
end


function _gather_band_energies(local_Enk, MPI_krange, Nk_total, comm)

    myrank = MPI.Comm_rank(comm)
    Nstate, spinsize, _ = size(local_Enk)
    recvcounts = [Nstate*spinsize*length(krange) for krange in MPI_krange]

    if myrank == 0
        Enk = Array{Float64}(undef, Nstate, spinsize, Nk_total)
        MPI.Gatherv!(local_Enk, MPI.VBuffer(Enk, recvcounts), comm; root=0)
        return Enk
    else
        MPI.Gatherv!(local_Enk, nothing, comm; root=0)
        return nothing
    end
end



function _gather_evec(local_EVec, MPI_krange, Nk_total, comm)

    myrank = MPI.Comm_rank(comm)
    Nstate, Norb, EVec_spinsize, _ = size(local_EVec)
    recvcounts = [Nstate*Norb*EVec_spinsize*length(krange) for krange in MPI_krange]

    if myrank == 0
        EVec = Array{Float32}(undef, Nstate, Norb, EVec_spinsize, Nk_total)
        MPI.Gatherv!(local_EVec, MPI.VBuffer(EVec, recvcounts), comm; root=0)
        return EVec
    else
        MPI.Gatherv!(local_EVec, nothing, comm; root=0)
        return nothing
    end
end
