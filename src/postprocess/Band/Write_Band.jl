function Write_BANDDAT(spinsize, filename, kpath_start, kpath_end, kpath_Nk, Nkpath, Nfsize, Enk, ChemP, rtv)

    k1 = 0.0
    k2 = 0.0
    k3 = 0.0
    klen = 0.0
    for spin = 1:spinsize
        Band_file = open("$filename.BANDDAT$(spin)", "w")
        for μ = 1:Nfsize
            Sum = 0.0
            for ik = 1:Nkpath
                k1_tmp = kpath_start[ik][1]
                k2_tmp = kpath_start[ik][2]
                k3_tmp = kpath_start[ik][3]
                for ipath = 1:kpath_Nk[ik]
                    k1 = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
                    k2 = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
                    k3 = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)
                
                    klen = norm((k1-k1_tmp)*rtv[1,:] + (k2-k2_tmp)*rtv[2,:] + (k3-k3_tmp)*rtv[3,:])
                    @printf(Band_file, "%lf %5.12f\n", klen+Sum, (Enk[ik][ipath][spin][μ]-ChemP)*eV2Hartree)
                end
                Sum += klen
                print(Band_file, "\n\n")
            end
        end
        close(Band_file)
    end
end


function Write_GNUBAND(spinsize, filename, Nkpath, kpath, kname, rtv)

    x = zeros(Float64, Nkpath+1)
    Sum = 0.0
    for ik = 1:Nkpath
        k1 = kpath[ik+1][1] - kpath[ik][1]
        k2 = kpath[ik+1][2] - kpath[ik][2]
        k3 = kpath[ik+1][3] - kpath[ik][3]
        Sum += norm(k1*rtv[1,:] + k2*rtv[2,:] + k3*rtv[3,:])
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
    if kname[Nkpath+1]=="G"
        println(gnu_file, "\"{/Symbol G}\" x$(Nkpath+1))")
    else
        println(gnu_file, "\"$(kname[Nkpath+1])\" x$(Nkpath+1))")
    end
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