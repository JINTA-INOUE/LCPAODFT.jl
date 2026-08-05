function Write_BANDDAT(filename::String, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, ChemP, Recvecs)

    size(Enk, 1) == Nfsize || error("unexpected number of band energies")
    size(Enk, 2) >= spin || error("requested spin channel is unavailable")
    size(Enk, 3) == sum(kpath_Nk) || error("unexpected number of k points")

    open("$filename.BANDDAT$(spin)", "w") do Band_file
        for μ = 1:Nfsize
            Sum = 0.0
            global_k = 0
            for ik = 1:Nkpath
                k1_tmp = kpath_start[ik][1]
                k2_tmp = kpath_start[ik][2]
                k3_tmp = kpath_start[ik][3]
                klen = 0.0
                for ipath = 1:kpath_Nk[ik]
                    global_k += 1
                    fraction = (ipath - 1) / (kpath_Nk[ik] - 1)
                    k1 = kpath_start[ik][1] +
                         (kpath_end[ik][1] - kpath_start[ik][1]) * fraction
                    k2 = kpath_start[ik][2] +
                         (kpath_end[ik][2] - kpath_start[ik][2]) * fraction
                    k3 = kpath_start[ik][3] +
                         (kpath_end[ik][3] - kpath_start[ik][3]) * fraction

                    klen = norm(
                        (k1 - k1_tmp) * Recvecs[1, :] +
                        (k2 - k2_tmp) * Recvecs[2, :] +
                        (k3 - k3_tmp) * Recvecs[3, :],
                    )
                    energy = (Enk[μ, spin, global_k] - ChemP) * eV2Hartree
                    @printf(Band_file, "%5.12f %5.12f\n", klen + Sum, energy)
                end
                Sum += klen
                print(Band_file, "\n\n")
            end
        end
    end
end


function Write_GNUBAND(filename::String, spinsize, Nkpath, kpath, kname, Recvecs)

    x = zeros(Float64, Nkpath+1)
    Sum = 0.0
    for ik = 1:Nkpath
        k1 = kpath[ik+1][1] - kpath[ik][1]
        k2 = kpath[ik+1][2] - kpath[ik][2]
        k3 = kpath[ik+1][3] - kpath[ik][3]
        Sum += norm(k1*Recvecs[1,:] + k2*Recvecs[2,:] + k3*Recvecs[3,:])
        x[ik+1] = Sum
    end

    
    gnu_file = open("$filename.plt", "w")
    println(gnu_file, "# set terminal postscript eps enhanced color")
    println(gnu_file, "set size 1,1")
    println(gnu_file, "")
    println(gnu_file, "set xlabel font \"Arial,15\"")
    println(gnu_file, "set ylabel font \"Arial,15\"")
    println(gnu_file, "set tics font \"Arial, 12\"")
    println(gnu_file, "set ylabel 'Energy (eV)'")
    println(gnu_file, "")
    println(gnu_file, "set xlabel offset 0,0")
    println(gnu_file, "set ylabel offset -2,0")
    println(gnu_file, "set lmargin 12")
    println(gnu_file, "set bmargin 2")
    println(gnu_file, "")
    println(gnu_file, "# set ytics 3")
    println(gnu_file, "set title font \"Arial, 20\"")
    println(gnu_file, "unset key")
    println(gnu_file, "")
    for ik = 1:Nkpath+1
        println(gnu_file, "x$(ik) = $(x[ik])")
    end
    println(gnu_file, "")
    println(gnu_file, "ymin = -15.0")
    println(gnu_file, "ymax = 15.0")
    println(gnu_file, "")
    println(gnu_file, "set zeroaxis lt 0")
    println(gnu_file, "")
    println(gnu_file, "set xra [0.000000:x$(Nkpath+1)]")
    println(gnu_file, "set yra [ymin:ymax]")
    print(gnu_file, "set xtics (")
    for ik = 1:Nkpath+1
        if kname[ik]=="G"
            print(gnu_file, "\"{/Symbol G}\" x$ik,")
        else
            print(gnu_file, "\"$(kname[ik])\" x$ik,")
        end
    end
    println(gnu_file, ")")
    for ik = 1:Nkpath-1
        println(gnu_file, "set arrow $ik nohead from x$(ik+1), ymin to x$(ik+1), ymax")
    end
    println(gnu_file, "set arrow $Nkpath nohead from 0, ymin to 0, x$(Nkpath+1), 0")
    println(gnu_file, "set size 1")
    println(gnu_file, "set origin 0,0")
    println(gnu_file, "plot \"$filename.BANDDAT1\" using 1:2 with lines linewidth 3")
    if spinsize == 2
        println(gnu_file, "# replot \"$filename.BANDDAT2\" using 1:2 with lines linewidth 3")
    end
    println(gnu_file, "set term pdf size 5in, 4in")
    println(gnu_file, "set output \"$filename.pdf\"")
    println(gnu_file, "replot")
    close(gnu_file)
end
