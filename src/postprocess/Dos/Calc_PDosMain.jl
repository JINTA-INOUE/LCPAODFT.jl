function Calc_PDosMain(filename::String, material::LCPAO_model, Enk, EVec, neg, kmesh, Dos_Erange)

    println("<PDosMain>  Generate Projecter Density of State using Tetrahedron method")

    Natom = material.Natom
    Nspecies = material.Nspecies
    atom2spe = material.atom2spe
    Atoms_pao = material.Atoms_pao
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    spinsize = ifelse(SpinPol == "off", 1, 2)


    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    kindex = zeros(Int64, Nkpt, 3)
    kindex2 = zeros(Int64, kmesh1, kmesh2, kmesh3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
        kindex2[ik,jk,kk] = kp
    end


    DosEmin, DosEmax = Dos_Erange
    de = 0.01/eV2Hartree
    Dos_N = floor(Int64, (DosEmax-DosEmin)/de)  # length of Dos

    DosE = zeros(Float64, Dos_N)
    Dos = zeros(Float64, Dos_N, spinsize)
    for ie = 1:Dos_N
        DosE[ie] = DosEmin + (DosEmax-DosEmin)*(ie-1)/(Dos_N-1)
    end


    Spe_orb = Get_Atoms_data(Atoms_pao)[3]
    Spe_MaxL_Basis = zeros(Int64, Nspecies)
    Spe_Num_Basis = Vector{Vector{Int64}}(undef, Nspecies)
    Atom_Num_Basis = Vector{Vector{Int64}}(undef, Natom)
    for spe = 1:Nspecies
        Spe_Num_Basis[spe] = zeros(Int64, 4)
    end

    for atom = 1:Natom
        Atom_Num_Basis[atom] = zeros(Int64, 4)
    end


    for spe = 1:Nspecies
        Spe_MaxL_Basis[spe], Num_Basis = get_ialpha_index(Spe_orb[spe])
        for l = 1:Spe_MaxL_Basis[spe]+1
            Spe_Num_Basis[spe][l] = Num_Basis[l]
        end
    end

    for atom = 1:Natom
        spe = atom2spe[atom]
        for l = 1:Spe_MaxL_Basis[spe]+1
            Atom_Num_Basis[atom][l] = Spe_Num_Basis[spe][l]
        end
    end


    Spe_Num_Relation = Vector{Vector{Int64}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_Num_Relation[spe] = zeros(Int64, maximum(Total_NumOrbs))
        id = 0
        for L = 0:3
            if Spe_Num_Basis[spe][L+1] > 0
                for _ = 0:Spe_Num_Basis[spe][L+1]-1, l = 0:2*L
                    id += 1
                    Spe_Num_Relation[spe][id] = 100*L + l + 1
                end
            end
        end
    end


    Dos = Vector{Vector{Vector{Vector{Float64}}}}(undef, spinsize)
    for spin = 1:spinsize
        Dos[spin] = Vector{Vector{Vector{Float64}}}(undef, Natom)
        for atom = 1:Natom
            Dos[spin][atom] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                Dos[spin][atom][ist] = zeros(Float64, Dos_N)
            end
        end
    end



    cell_e = zeros(Float64, 8)
    cell_a = zeros(Float64, Natom, maximum(Total_NumOrbs), 8)
    tetra_e = zeros(Float64, 4)
    tetra_a = zeros(Float64, 4)
    tetra_id = [1 2 3 6; 2 3 4 6; 3 4 6 8; 1 3 5 6; 3 5 6 7; 3 6 7 8]



    # tetrahedron
    for spin = 1:spinsize, ieg = 1:neg
        for ik = 1:Nkpt
            i = kindex[ik,1]-1
            j = kindex[ik,2]-1
            k = kindex[ik,3]-1

            for i_in = 0:1, j_in = 0:1, k_in = 0:1
                ii = mod(i+i_in, kmesh1)+1
                jj = mod(j+j_in, kmesh2)+1
                kk = mod(k+k_in, kmesh3)+1
                cell_e[4*i_in+2*j_in+k_in+1] = Enk[spin][ii][jj][kk][ieg]
                
                kp = kindex2[ii,jj,kk]
                for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                    rval = EVec[spin][kp][atom][ist][ieg]
                    cell_a[atom,ist,4*i_in+2*j_in+k_in+1] = rval
                end
            end


            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                for itetra = 1:6
                    for ic = 1:4
                        tetra_e[ic] = cell_e[tetra_id[itetra,ic]]
                        tetra_a[ic] = cell_a[atom,ist,tetra_id[itetra,ic]]
                    end
                    
    
                    OrderE!(tetra_e, tetra_a, 4)
                    
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
                            result = ATM_Spectrum(tetra_e, tetra_a, DosE[ie+1])
                            Dos[spin][atom][ist][ie+1] += result
                        end
                    end
                end
            end
        end
    end



    # Normalize
    for spin = 1:spinsize, atom = 1:Natom, ist = 1:Total_NumOrbs[atom], ie = 1:Dos_N
        Dos[spin][atom][ist][ie] = Dos[spin][atom][ist][ie]/Nkpt/6/eV2Hartree
    end

    
    Calc_PDos_Atom_proj(filename, material, Dos_Erange, DosE, Dos)
    Calc_PDos_Orbital_proj(filename, material, Spe_Num_Relation, Spe_Num_Basis, Dos_Erange, DosE, Dos)
    Write_PDos_gnuplot(filename, Natom, Dos_Erange)
end


function Calc_PDosMain(filename::String, material::CWF_model, Enk, EVec, neg, kmesh, Dos_Erange)

    println("<PDosMain>  Generate Projecter Density of State using Tetrahedron method")

    gsize = material.gsize
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol == "off", 1, 2)


    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    kindex = zeros(Int64, Nkpt, 3)
    kindex2 = zeros(Int64, kmesh1, kmesh2, kmesh3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
        kindex2[ik,jk,kk] = kp
    end


    DosEmin, DosEmax = Dos_Erange
    de = 0.01/eV2Hartree
    Dos_N = floor(Int64, (DosEmax-DosEmin)/de)  # length of Dos

    DosE = zeros(Float64, Dos_N)
    Dos = zeros(Float64, Dos_N, spinsize)
    for ie = 1:Dos_N
        DosE[ie] = DosEmin + (DosEmax-DosEmin)*(ie-1)/(Dos_N-1)
    end



    Dos = Vector{Vector{Vector{Float64}}}(undef, spinsize)
    for spin = 1:spinsize
        Dos[spin] = Vector{Vector{Float64}}(undef, gsize)
        for ist = 1:gsize
            Dos[spin][ist] = zeros(Float64, Dos_N)
        end
    end



    cell_e = zeros(Float64, 8)
    cell_a = zeros(Float64, gsize, 8)
    tetra_e = zeros(Float64, 4)
    tetra_a = zeros(Float64, 4)
    tetra_id = [1 2 3 6; 2 3 4 6; 3 4 6 8; 1 3 5 6; 3 5 6 7; 3 6 7 8]



    # tetrahedron
    for spin = 1:spinsize, ieg = 1:neg
        for ik = 1:Nkpt, ist = 1:gsize
            i = kindex[ik,1]-1
            j = kindex[ik,2]-1
            k = kindex[ik,3]-1

            for i_in = 0:1, j_in = 0:1, k_in = 0:1
                ii = mod(i+i_in, kmesh1)+1
                jj = mod(j+j_in, kmesh2)+1
                kk = mod(k+k_in, kmesh3)+1
                cell_e[4*i_in+2*j_in+k_in+1] = Enk[spin][ii][jj][kk][ieg]
                
                kp = kindex2[ii,jj,kk]
                rval = EVec[spin][kp][ist][ieg]
                cell_a[4*i_in+2*j_in+k_in+1] = rval
            end


            for itetra = 1:6
                for ic = 1:4
                    tetra_e[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_a[ic] = cell_a[tetra_id[itetra,ic]]
                end
                    
    
                OrderE!(tetra_e, tetra_a, 4)
                    
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
                        result = ATM_Spectrum(tetra_e, tetra_a, DosE[ie+1])
                        Dos[spin][ist][ie+1] += result
                    end
                end
            end
        end
    end



    # Normalize
    for spin = 1:spinsize, ist = 1:gsize, ie = 1:Dos_N
        Dos[spin][ist][ie] = Dos[spin][ist][ie]/Nkpt/6/eV2Hartree
    end

    
    Calc_PDos_Orbital_proj(filename, material, Dos_Erange, DosE, Dos)
    # Write_PDos_gnuplot(filename, Dos_Erange)
end