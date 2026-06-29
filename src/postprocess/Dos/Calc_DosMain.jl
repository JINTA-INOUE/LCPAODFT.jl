function Calc_DosMain(filename::String, material::Union{LCPAO_model,CWF_model}, Enk, neg, kmesh, Dos_Erange; de_Dos=0.01)

    println("<DosMain>  Generate Density of State using Tetrahedron method")
    
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)


    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    DosEmin, DosEmax = Dos_Erange
    de = de_Dos/eV2Hartree
    Dos_N = floor(Int64, (DosEmax-DosEmin)/de)  # length of Dos

    DosE = zeros(Float64, Dos_N)
    Dos = zeros(Float64, Dos_N, spinsize)
    for ie = 1:Dos_N
        DosE[ie] = DosEmin + (DosEmax-DosEmin)*(ie-1)/(Dos_N-1)
    end

    cell_e = zeros(Float64, 8)
    tetra_e = zeros(Float64, 4)
    tetra_id = [1 2 3 6; 2 3 4 6; 3 4 6 8; 1 3 5 6; 3 5 6 7; 3 6 7 8]



    # tetrahedron
    for spin = 1:spinsize, ieg = 1:neg
        for ik = 1:Nkpt
            i = kindex[ik,1]-1
            j = kindex[ik,2]-1
            k = kindex[ik,3]-1

            for i_in = 0:1, j_in = 0:1, k_in = 0:1
                ii = mod(i+i_in,kmesh1)+1
                jj = mod(j+j_in,kmesh2)+1
                kk = mod(k+k_in,kmesh3)+1
                cell_e[4*i_in+2*j_in+k_in+1] = Enk[spin][ii][jj][kk][ieg]
            end

            for itetra = 1:6
                for ic = 1:4
                    tetra_e[ic] = cell_e[tetra_id[itetra,ic]]
                end
                

                OrderE0!(tetra_e, 4)
                
                x = (tetra_e[1]-DosEmin)/(DosEmax-DosEmin)*(Dos_N-1)-1
                iemin = trunc(Int, x)
                x = (tetra_e[4]-DosEmin)/(DosEmax-DosEmin)*(Dos_N-1)+1
                iemax = trunc(Int, x)

                if iemin < 0
                    iemin = 0
                end
                if iemax >= Dos_N
                    iemax = Dos_N - 1
                end
                if 0 <= iemin < Dos_N && 0 <= iemax < Dos_N
                    for ie = iemin:iemax
                        result = ATM_Dos(tetra_e, DosE[ie])
                        Dos[ie,spin] += result
                    end
                end
            end
        end
    end


    h = (DosEmax - DosEmin)/(Dos_N-1) * eV2Hartree
    ssum = zeros(Float64, Dos_N, spinsize)

    factor = 1/Nkpt/6/eV2Hartree

    # Normalize
    for spin = 1:spinsize
        for ie = 1:Dos_N
            Dos[ie,spin] = Dos[ie,spin]*factor
        end

        # Calc Total Dos
        for q = 1:Dos_N
            s1 = 0.0
            s2 = 0.0
            for ie = 2:2:q-1
                s1 += Dos[ie,spin]
            end
            for ie = 3:2:q-1
                s2 += Dos[ie,spin]
            end
            ssum[q,spin] = (Dos[1,spin] + 4*s1 + 2*s2 + Dos[q,spin])*h/3
        end
    end


    Write_Dos_Tetrahedron(filename, SpinPol, Dos_N, DosE, Dos, ssum)
    Write_Dos_gnuplot(filename, SpinPol, Dos_Erange)
end