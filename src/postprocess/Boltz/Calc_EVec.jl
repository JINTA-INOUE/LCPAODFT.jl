@timeit timer "Calc_EVec" function Calc_EVec(
    boltz_setup::Boltz_Setup, kpoints::BoltzKPoints, Cnk)

    material = boltz_setup.material
    boltz_setup.decomp || error("Calc_EVec requires decomp=true")
    boltz_setup.mat_type == "CWF" ||
        error("transport decomposition is supported only for CWF models")

    Nwann = Int(material.Ngsize)
    spinsize = material.SpinPol == "on" ? 2 : 1
    EVec = zeros(Float32, Nwann, Nwann, spinsize, kpoints.MPI_Nkpt)

    for spin = 1:spinsize, local_k = 1:kpoints.MPI_Nkpt
        Cn = Cnk[spin][local_k]
        @inbounds for orbital = 1:Nwann, band = 1:Nwann
            coefficient = Cn[orbital, band]
            weight = real(conj(coefficient) * coefficient)
            EVec[orbital, band, spin, local_k] = Float32(weight)
        end
    end
    return EVec
end
