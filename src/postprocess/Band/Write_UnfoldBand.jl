function Write_UnfoldBANDDAT(filename::String, SpinPol, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nstate, Enk, EVec, ChemP, Recvecs)

    if SpinPol == "nc"
        Nstate_half = div(Nstate, 2)
    end

    open("$(filename)_unfold.BANDDAT$(spin)", "w") do Band_file
        for μ = 1:Nstate
            Sum = 0.0
            global_k = 0
            for ik = 1:Nkpath
                k1_tmp, k2_tmp, k3_tmp = kpath_start[ik]
                klen = 0.0
                for ipath = 1:kpath_Nk[ik]
                    global_k += 1
                    fraction = (ipath - 1)/(kpath_Nk[ik] - 1)
                    k1 = kpath_start[ik][1] + (kpath_end[ik][1] - kpath_start[ik][1])*fraction
                    k2 = kpath_start[ik][2] + (kpath_end[ik][2] - kpath_start[ik][2])*fraction
                    k3 = kpath_start[ik][3] + (kpath_end[ik][3] - kpath_start[ik][3])*fraction

                    klen = norm((k1 - k1_tmp)*Recvecs[1,:] + (k2 - k2_tmp)*Recvecs[2,:] + (k3 - k3_tmp)*Recvecs[3,:])
                    energy = (Enk[μ,spin,global_k] - ChemP)*eV2Hartree
                    @printf(Band_file, "%5.12f %5.12f", klen + Sum, energy)

                    if SpinPol ∈ ("off", "on")
                        for ist = 1:Nstate
                            @printf(Band_file, " %5.12f", EVec[μ,ist,spin,global_k])
                        end
                    elseif SpinPol == "nc"
                        for ist = 1:Nstate_half
                            @printf(Band_file, " %5.12f", EVec[μ,ist,1,global_k] + EVec[μ,ist,2,global_k])
                        end
                    end
                    @printf(Band_file, "\n")
                end
                Sum += klen
                print(Band_file, "\n\n")
            end
        end
    end
end


function Write_GNUUNFOLDBAND(filename::String, spinsize, Nkpath, kpath, kname, Recvecs)

    x = zeros(Float64, Nkpath+1)
    Sum = 0.0
    for ik = 1:Nkpath
        k1 = kpath[ik+1][1] - kpath[ik][1]
        k2 = kpath[ik+1][2] - kpath[ik][2]
        k3 = kpath[ik+1][3] - kpath[ik][3]
        Sum += norm(k1*Recvecs[1,:] + k2*Recvecs[2,:] + k3*Recvecs[3,:])
        x[ik+1] = Sum
    end

    
    gnu_file = open("$(filename)_unfold.plt", "w")
    println(gnu_file, "# set encoding iso_8859_1")
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
    println(gnu_file, "")
    println(gnu_file, "# set palette defined ( 0 \'#bebebe\', 0.25 \'#ffee00\', 0.5 \'#ff7000\', 0.75 \'#ee0000\', 1 \'#7f0000\')")
    println(gnu_file, "set palette defined ( 0.0 '#000090', 0.125 '#000fff',0.25 '#0090ff',0.375 '#0fffee',0.5 '#90ff70',0.625 '#ffee00',0.75 '#ff7000',0.875 '#ee0000',1.0 '#7f0000')")
    println(gnu_file, "set cbrange [0:1]")
    println(gnu_file, "")
    println(gnu_file, "# plot \"$(filename)_unfold.BANDDAT1\" using 1:2 with lines linewidth 3")
    println(gnu_file, "plot \"$(filename)_unfold.BANDDAT1\" using 1:2:(\$3) w lp lw 2 pt 7 ps 1 lc palette")
    if spinsize == 2
        println(gnu_file, "# replot \"$(filename)_unfold.BANDDAT2\" using 1:2 with lines linewidth 3")
    end
    println(gnu_file, "set term pdf size 5in, 4in")
    println(gnu_file, "set output \"$(filename)_unfold.pdf\"")
    println(gnu_file, "replot")
    close(gnu_file)
end
