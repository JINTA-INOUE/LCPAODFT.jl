function Calc_EVec_Dos(material::LCPAO_model, kpoints::DosKPoints, Cnk,
                       iemin, iemax)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    fsize = sum(material.Total_NumOrbs)
    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    neg = iemax - iemin + 1

    # (band, orbital, projected spin, rank-local k point)
    EVec = zeros(Float32, neg, fsize, spinsize, MPI_Nkpt)

    if SpinPol in ("off", "on")
        Get_EVec_Collinear!(material, kpoints, Cnk, EVec, iemin, iemax)
    elseif SpinPol == "nc"
        Get_EVec_NonCollinear!(material, kpoints, Cnk, EVec, iemin, iemax)
    else
        error("unsupported spin polarization: $SpinPol")
    end

    return EVec
end


function Get_EVec_Collinear!(material::LCPAO_model, kpoints::DosKPoints,
                              Cnk, EVec, iemin, iemax)

    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    fsize = sum(Total_NumOrbs)
    spinsize = ifelse(SpinPol == "off", 1, 2)
    OLP = material.OLP

    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    MPI_kpts = kpoints.MPI_kpts
    S = zeros(ComplexF64, fsize, fsize)
    SD = zeros(Float32, fsize)

    for spin = 1:spinsize, local_k = 1:MPI_Nkpt
        C = Cnk[spin][local_k]

        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn,
                   atv_ijk, MPI_kpts[local_k])
        @inbounds for band = iemin:iemax
            fill!(SD, 0.0f0)
            for atom = 1:Natom, jatom = 1:Natom
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                Anum = MP[atom]
                Bnum = MP[jatom]

                for ist = 1:NO0, jst = 1:NO1
                    tmp = conj(C[Anum + ist, band]) * C[Bnum + jst, band]
                    SD[Anum + ist] += real(tmp * S[Anum + ist, Bnum + jst])
                end
            end

            window_band = band - iemin + 1
            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                orbital = MP[atom] + ist
                EVec[window_band, orbital, spin, local_k] = SD[orbital]
            end
        end
    end
end


function Get_EVec_NonCollinear!(material::LCPAO_model, kpoints::DosKPoints,
                                 Cnk, EVec, iemin, iemax)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    fsize = sum(Total_NumOrbs)
    OLP = material.OLP
    Angle_spin = Get_Angle_spin(material)

    sin_theta = sin.([angle[1] for angle in Angle_spin])
    cos_theta = cos.([angle[1] for angle in Angle_spin])
    sin_phi = sin.([angle[2] for angle in Angle_spin])
    cos_phi = cos.([angle[2] for angle in Angle_spin])

    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    MPI_kpts = kpoints.MPI_kpts
    S = zeros(ComplexF64, fsize, fsize)
    SD = zeros(Float32, 2 * fsize)

    for local_k = 1:MPI_Nkpt
        C = Cnk[1][local_k]

        HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn,
                   atv_ijk, MPI_kpts[local_k])
        @inbounds for band = iemin:iemax
            fill!(SD, 0.0f0)
            for atom = 1:Natom, jatom = 1:Natom
                sit = sin_theta[atom]
                cot = cos_theta[atom]
                sip = sin_phi[atom]
                cop = cos_phi[atom]
                NO0 = Total_NumOrbs[atom]
                NO1 = Total_NumOrbs[jatom]
                Anum = MP[atom]
                Bnum = MP[jatom]

                for ist = 1:NO0, jst = 1:NO1
                    svalue = S[Anum + ist, Bnum + jst]
                    tmp_uu = real(conj(C[Anum + ist, band]) *
                                  C[Bnum + jst, band] * svalue)
                    tmp_dd = real(conj(C[Anum + ist + fsize, band]) *
                                  C[Bnum + jst + fsize, band] * svalue)
                    tmp_ud = conj(C[Anum + ist, band]) *
                             C[Bnum + jst + fsize, band] * svalue
                    tmp_ud_real = real(tmp_ud)
                    tmp_ud_imag = imag(tmp_ud)

                    SD[Anum + ist] +=
                        0.5 * (tmp_uu + tmp_dd) +
                        0.5 * cot * (tmp_uu - tmp_dd) +
                        (tmp_ud_real * cop - tmp_ud_imag * sip) * sit
                    SD[Anum + ist + fsize] +=
                        0.5 * (tmp_uu + tmp_dd) -
                        0.5 * cot * (tmp_uu - tmp_dd) -
                        (tmp_ud_real * cop - tmp_ud_imag * sip) * sit
                end
            end

            window_band = band - iemin + 1
            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                orbital = MP[atom] + ist
                EVec[window_band, orbital, 1, local_k] = SD[orbital]
                EVec[window_band, orbital, 2, local_k] = SD[orbital + fsize]
            end
        end
    end
end


function Calc_EVec_Dos(material::CWF_model, kpoints::DosKPoints, Cnk,
                       iemin, iemax)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)
    cnk_spinsize = ifelse(SpinPol == "on", 2, 1)
    gsize = Int(material.gsize)
    MPI_Nkpt = Int(kpoints.MPI_Nkpt)
    neg = iemax - iemin + 1

    EVec = zeros(Float32, neg, gsize, spinsize, MPI_Nkpt)

    if SpinPol in ("off", "on")
        for spin = 1:cnk_spinsize, local_k = 1:MPI_Nkpt
            C = Cnk[spin][local_k]
            @inbounds for ist = 1:gsize, band = iemin:iemax
                EVec[band - iemin + 1, ist, spin, local_k] =
                    abs2(C[ist, band])
            end
        end
    elseif SpinPol == "nc"
        for local_k = 1:MPI_Nkpt
            C = Cnk[1][local_k]
            @inbounds for ist = 1:gsize, band = iemin:iemax
                window_band = band - iemin + 1
                EVec[window_band, ist, 1, local_k] = abs2(C[ist, band])
                EVec[window_band, ist, 2, local_k] =
                    abs2(C[ist + gsize, band])
            end
        end
    else
        error("unsupported spin polarization: $SpinPol")
    end

    return EVec
end
