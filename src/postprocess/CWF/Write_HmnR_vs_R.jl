function Write_HmnR_vs_R(cwf_setup::CWF_Setup, filepath::String)

    cwf_model = Load_CWF_model(filepath)

    spinsize = cwf_model.spinsize
    Latvecs = cwf_model.Latvecs
    Ngsize = cwf_model.Ngsize
    NCell = cwf_model.NCell
    cell_list_ijk = cwf_model.cell_list_ijk
    HmnR = cwf_model.HmnR
    Wannier_Guide = cwf_setup.Wannier_Guide
    filename = cwf_setup.filename


    for spin = 1:spinsize
        data = open(filename*"HmnR_vs_R$spin.dat", "w")
        for cell = 1:NCell
            l1, l2, l3 = cell_list_ijk[cell]
            Rx = l1*Latvecs[1,1] + l2*Latvecs[2,1] + l3*Latvecs[3,1]
            Ry = l1*Latvecs[1,2] + l2*Latvecs[2,2] + l3*Latvecs[3,2]
            Rz = l1*Latvecs[1,3] + l2*Latvecs[2,3] + l3*Latvecs[3,3]
            for m = 1:Ngsize, n = 1:Ngsize
                xm, ym, zm = Wannier_Guide[m]
                xn, yn, zn = Wannier_Guide[n]

                distance = (xn-xm)^2 + (yn-ym)^2 +(zn-zm)^2 + Rx^2 + Ry^2 + Rz^2
                distance = sqrt(distance)/Ang_to_bohr
                @printf(data, "%5.12f  %5.12f  %5.12f  %5.12f\n", distance, log10(abs(real(HmnR[n,m,cell,spin]))), log10(abs(imag(HmnR[n,m,cell,spin]))), log10(abs(HmnR[n,m,cell,spin])))
            end
        end
        close(data)
    end
end