function Write_win(filename::String, kmesh, BANDNUM, WANNUM, material::LCPAO_model)

    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Natom = material.Natom
    Latvecs = material.Latvecs
    Recvecs = material.Recvecs
    Gxyz_AU = material.Gxyz
    Atoms_symbol = material.Atoms_symbol
    kmesh1, kmesh2, kmesh3 = kmesh


    Gxyz_frac = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Gxyz_frac[atom] = zeros(Float64, 3)
    end

    for atom = 1:Natom
        Gxyz_frac[atom][1] = dot(Gxyz_AU[atom], Recvecs[1,:])*0.5/pi
        Gxyz_frac[atom][2] = dot(Gxyz_AU[atom], Recvecs[2,:])*0.5/pi
        Gxyz_frac[atom][3] = dot(Gxyz_AU[atom], Recvecs[3,:])*0.5/pi
                
        for i = 1:3
            tmp = floor(Int64, Gxyz_frac[atom][i])
            if Gxyz_frac[atom][i] > 1.0
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]-tmp)
            elseif Gxyz_frac[atom][i] < -1e-13
                Gxyz_frac[atom][i] = abs(Gxyz_frac[atom][i]+abs(tmp)+1)
            end
        end
    end


    for spin = 1:spinsize
        if SpinPol == "on"
            data = open(filename*"_$spin.win", "w")
        else
            data = open(filename*".win", "w")
        end

        @printf(data, "num_bands %d\n", BANDNUM)
        @printf(data, "num_wann %d\n", WANNUM)
        @printf(data, "\n")
        if SpinPol == "nc"
            @printf(data, "spinors = true\n")
            @printf(data, "\n")
        end
        @printf(data, "dis_num_iter 0\n")
        @printf(data, "num_iter 0\n")
        @printf(data, "! dis_conv_tol = 1.0e-12\n")
        @printf(data, "! conv_tol = 1.0e-6\n")
        @printf(data, "\n")
        @printf(data, "write_rmn  = false\n")
        @printf(data, "write_r2mn = false\n")
        @printf(data, "write_hr = false\n")
        @printf(data, "write_u_matrices = false\n")
        @printf(data, "\n")
        @printf(data, "! iprint = 1\n")
        @printf(data, "! kmesh_tol = 1.0e-5\n")
        @printf(data, "\n")
        @printf(data, "bands_plot F\n")
        @printf(data, "bands_num_points 200\n")
        @printf(data, "bands_plot_format gnuplot\n")
        @printf(data, "!begin kpoint_path\n")
        @printf(data, "!end kpoint_path\n")
        @printf(data, "\n")
        @printf(data, "begin unit_cell_cart\n")
        @printf(data, "Ang\n")
        for i = 1:3
            @printf(data, "%18.12f %18.12f %18.12f\n", Latvecs[i,1]/Ang_to_bohr, Latvecs[i,2]/Ang_to_bohr, Latvecs[i,3]/Ang_to_bohr)
        end
        @printf(data, "end unit_cell_cart\n")
        @printf(data, "\n")
        @printf(data, "begin atoms_frac\n")
        for atom = 1:Natom
            @printf(data, "%4s %18.14f %18.14f %18.14f\n", Atoms_symbol[atom], Gxyz_frac[atom][1], Gxyz_frac[atom][2], Gxyz_frac[atom][3])
        end
        @printf(data, "end atoms_frac\n")
        @printf(data, "\n")
        @printf(data, "mp_grid %d %d %d\n", kmesh1, kmesh2, kmesh3)
        @printf(data, "\n")
        @printf(data, "begin_kpoints\n")
        for i = 0:kmesh1-1, j = 0:kmesh2-1, k = 0:kmesh3-1
            k1 = ifelse(kmesh1==1, 0.0, i/kmesh1)
            k2 = ifelse(kmesh2==1, 0.0, j/kmesh2)
            k3 = ifelse(kmesh3==1, 0.0, k/kmesh3)
            @printf(data, "%18.14f %18.14f %18.14f\n", k1, k2, k3)
        end
        @printf(data, "end_kpoints\n")
        close(data)
    end
end