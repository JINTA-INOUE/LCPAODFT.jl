function Write_Dos_Tetrahedron(filename::String, SpinPol, Dos_N, DosE, Dos, ssum)

    Dos_data = open("$(filename).DOS.Tetrahedron", "w")
    if SpinPol ∈ ("on", "nc")
        for ie = 1:Dos_N
            @printf(Dos_data, "%lf %lf %lf %lf %lf\n", DosE[ie]*eV2Hartree, Dos[ie,1], -Dos[ie,2], ssum[ie,1], ssum[ie,2])
        end
    elseif SpinPol == "off"
        for ie = 1:Dos_N
            @printf(Dos_data, "%lf %lf %lf\n", DosE[ie]*eV2Hartree, Dos[ie,1]*2, ssum[ie,1]*2)
        end
    end
    close(Dos_data)
end


function Write_Dos_gnuplot(filename::String, SpinPol::String, Dos_Erange)

    gnu_file = open("$(filename).DOS.plt", "w")
    println(gnu_file, "# set terminal postscript eps enhanced color")
    println(gnu_file, "set size 1,1")
    println(gnu_file, "")
    println(gnu_file, "set xlabel font \"Arial,15\"")
    println(gnu_file, "set ylabel font \"Arial,15\"")
    println(gnu_file, "set tics font \"Arial, 10\"")
    println(gnu_file, "set xlabel 'Energy (eV)'")
    println(gnu_file, "set ylabel 'Dos (eV^{-1})'")
    println(gnu_file, "")
    println(gnu_file, "set xlabel offset 0,0")
    println(gnu_file, "set ylabel offset -2,0")
    println(gnu_file, "set lmargin 12")
    println(gnu_file, "set bmargin 2")
    println(gnu_file, "")
    # println(gnu_file, "set ytics 0.2")
    println(gnu_file, "set title font \"Arial, 20\"")
    println(gnu_file, "# unset key")
    println(gnu_file, "set grid")
    println(gnu_file, "")
    println(gnu_file, "# ymax = 1.0")
    println(gnu_file, "")
    println(gnu_file, "set xra [$(Dos_Erange[1]*eV2Hartree):$(Dos_Erange[2]*eV2Hartree)]")
    println(gnu_file, "# set yra [0.0:ymax]")
    println(gnu_file, "set size 1")
    println(gnu_file, "plot \"$(filename).DOS.Tetrahedron\" using 1:2 with lines linewidth 3")
    if SpinPol ∈ ("on", "nc")
        println(gnu_file, "replot \"$(filename).DOS.Tetrahedron\" using 1:3 with lines linewidth 3")
    end
    println(gnu_file, "")
    println(gnu_file, "set term pdf size 5in, 4in")
    println(gnu_file, "set output \"$(filename).DOS.pdf\"")
    println(gnu_file, "replot")
    close(gnu_file)
end



function Write_PDos_gnuplot(filename::String, Natom, Dos_Erange)

    gnu_file = open("$(filename).PDOS.plt", "w")
    println(gnu_file, "# set terminal postscript eps enhanced color")
    println(gnu_file, "set size 1,1")
    println(gnu_file, "")
    println(gnu_file, "set xlabel font \"Arial,15\"")
    println(gnu_file, "set ylabel font \"Arial,15\"")
    println(gnu_file, "set tics font \"Arial, 10\"")
    println(gnu_file, "set xlabel 'Energy (eV)'")
    println(gnu_file, "set ylabel 'Dos (eV^{-1})'")
    println(gnu_file, "")
    println(gnu_file, "set xlabel offset 0,0")
    println(gnu_file, "set ylabel offset -2,0")
    println(gnu_file, "set lmargin 12")
    println(gnu_file, "set bmargin 2")
    println(gnu_file, "")
    println(gnu_file, "# set ytics 1.0")
    println(gnu_file, "set title font \"Arial, 20\"")
    println(gnu_file, "# unset key")
    println(gnu_file, "set grid")
    println(gnu_file, "")
    println(gnu_file, "# ymin = 0.0")
    println(gnu_file, "# ymax = 5.0")
    println(gnu_file, "")
    println(gnu_file, "# set xra [$(Dos_Erange[1]*eV2Hartree):$(Dos_Erange[2]*eV2Hartree)]")
    println(gnu_file, "# set yra [ymin:ymax]")
    println(gnu_file, "set size 1")
    println(gnu_file, "plot \"$(filename).PDOS.Tetrahedron.atom1\" using 1:2 with lines linewidth 3 title \"\"")
    for atom = 2:Natom
        println(gnu_file, "# replot \"$(filename).PDOS.Tetrahedron.atom$(atom)\" using 1:2 with lines linewidth 3 title \"\"")
    end
    println(gnu_file, "")
    println(gnu_file, "set term pdf size 5in, 4in")
    println(gnu_file, "set output \"$(filename).PDOS_Tetrahedron_atom1.pdf\"")
    println(gnu_file, "replot")
    close(gnu_file)
end
