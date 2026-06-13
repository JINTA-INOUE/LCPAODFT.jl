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


    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]
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


function Calc_DosMain(filename::String, material::LCPAO_model, Enk, EVec, neg, kmesh, Dos_Erange)

    println("<DosMain>  Generate Density of State using Tetrahedron method")
    
    Natom = material.Natom
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
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


    DosEmin = Dos_Erange[1]
    DosEmax = Dos_Erange[2]

    de = 0.01/eV2Hartree
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
            wval = 0.0
            for atom = 1:Natom, ist = 1:Total_NumOrbs[atom]
                wval += EVec[spin][ik][atom][ist][ieg]
            end

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
                        Dos[ie,spin] += wval * result
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



function DosMain(filepath::String, kmesh, Erange::Vector{Float64}; mode::String="Dos")

    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)

    if nprocs > 1
        error("please run serial.")
    end

    if length(kmesh) ≠ 3
        error("pleas check kmesh")
    end

    if length(Erange) ≠ 2
        error("pleas check Erange")
    end


    mode = lowercase(mode)
    if mode ∉ ("all", "dos", "pdos")
        error("please check mode.")
    end


    filename, _ = splitext(basename(filepath))
    println("filename = $filename")
    
    material = Load_LCPAODFT_model(filepath)
    Print_LCPAO_model(filepath, material)

    Natom = material.Natom
    FNAN = material.FNAN
    natn = material.natn
    ncn = material.ncn
    atv_ijk = material.atv_ijk
    SpinPol = material.SpinPol
    Total_NumOrbs = material.Total_NumOrbs
    MP = material.MP
    ChemP = material.ChemP
    fsize = sum(Total_NumOrbs)
    Nfsize = ifelse(SpinPol ∈ ("off", "on"), fsize, 2*fsize)
    spinsize = ifelse(SpinPol == "off", 1, 2)

    OLP = material.OLP
    Hks = material.Hks
    iHks = material.iHks
    ChemP = material.ChemP

    Dos_Erange = zeros(Float64, 2)
    @. Dos_Erange = Erange/eV2Hartree


    kmesh1, kmesh2, kmesh3 = kmesh
    Nkpt = prod(kmesh)
    kpts = zeros(Float64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kpts[kp,1] = (ik-1)/kmesh1
        kpts[kp,2] = (jk-1)/kmesh2
        kpts[kp,3] = (kk-1)/kmesh3
    end



    # Find iemin, iemax at Γ point
    # iemin : minimal band index
    # iemax : maximum band index
    iemin = 1
    iemax = 1
    n1min = 1
    kpts_zeros = zeros(Float64, 3)
    EΓ = zeros(Float64, Nfsize)
    S = zeros(ComplexF64, Nfsize, Nfsize)
    H = zeros(ComplexF64, Nfsize, Nfsize)
    tmpH = zeros(ComplexF64, fsize, fsize)

    for spin = 1:spinsize

        if SpinPol ∈ ("off", "on")
            HS_matrix!(S, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            HS_matrix!(H, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        elseif SpinPol == "nc"
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts_zeros)
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            EΓ = eigvals(Hermitian(H), Hermitian(S))
        end

        iemin0 = 1
        n1min = ifelse(n1min<Nfsize, Nfsize, n1min)
        
        for μ = 1:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[1])
                iemin0 = μ - 1
                break
            end
        end

        iemin0 = ifelse(iemin0<1, 1, iemin0)
        iemax0 = Nfsize

        for μ = iemin0:Nfsize
            if EΓ[μ] > (ChemP + Dos_Erange[2])
                iemax0 = μ
                break
            end
        end

        iemax0 = ifelse(iemax0>Nfsize, Nfsize, iemax0)
        iemin = ifelse(iemin>iemin0, iemin0, iemin)
        iemax = ifelse(iemax<iemax0, iemax0, iemax)
    end

    if SpinPol ∈ ("off", "on")
        iemin -= max(div(fsize, 20), 10)
        iemax += max(div(fsize, 20), 10)
    elseif SpinPol == "nc"
        iemin -= max(div(fsize, 10), 10)
        iemax += max(div(fsize, 10), 10)
    end

    iemin = ifelse(iemin<1, 1, iemin)
    iemax = ifelse(iemax>Nfsize, Nfsize, iemax)
    iemax = ifelse(iemax>n1min, n1min, iemax)
    neg = iemax-iemin+1


    # Solve Eigen Problem
    Enk = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, spinsize)
    for spin = 1:spinsize
        Enk[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Enk[spin][ik] = Vector{Vector{Vector{Float64}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Enk[spin][ik][jk] = Vector{Vector{Float64}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Enk[spin][ik][jk][kk] = zeros(Float64, neg)
                end
            end
        end
    end
    Cnk = Vector{Vector{Matrix{ComplexF64}}}(undef, spinsize)
    for spin = 1:spinsize
        Cnk[spin] = Vector{Matrix{ComplexF64}}(undef, Nkpt)
        for ik = 1:Nkpt
            Cnk[spin][ik] = zeros(ComplexF64, Nfsize, Nfsize)
        end
    end


    Enk_tmp = zeros(Float64, Nfsize)
    if SpinPol ∈ ("off", "on")
        for spin = 1:spinsize
            kp = 0
            for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
                kp += 1
                HS_matrix!(S, H, OLP, Hks[spin], Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp,:])

                Enk_tmp, Cnk[spin][kp] = eigen(Hermitian(H), Hermitian(S))
                for μ = 1:neg
                    Enk[spin][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP
                end
            end
        end
    elseif SpinPol == "nc"
        kp = 0
        for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
            kp += 1
            HS_matrix_NC!(tmpH, H, Hks, iHks, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp,:])
            HS_matrix!(tmpH, OLP, Natom, Total_NumOrbs, MP, FNAN, natn, ncn, atv_ijk, kpts[kp,:])
            @. S[1:fsize, 1:fsize] = tmpH
            @. S[fsize+1:end, fsize+1:end] = tmpH
            Enk_tmp, Cnk[1][kp] = eigen(Hermitian(H), Hermitian(S))

            for μ = 1:neg
                Enk[1][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP
                Enk[2][ik][jk][kk][μ] = Enk_tmp[μ+iemin-1] - ChemP      # copy for EVec
            end
        end
    end

    

    EVec = Vector{Vector{Vector{Vector{Vector{Float32}}}}}(undef, spinsize)
    for spin = 1:spinsize
        EVec[spin] = Vector{Vector{Vector{Vector{Float32}}}}(undef, Nkpt)
        for ik = 1:Nkpt
            EVec[spin][ik] = Vector{Vector{Vector{Float32}}}(undef, Natom)
            for atom = 1:Natom
                EVec[spin][ik][atom] = Vector{Vector{Float32}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    EVec[spin][ik][atom][ist] = zeros(Float32, neg)
                end
            end
        end
    end

    
    if SpinPol ∈ ("off", "on")
        Get_EVec_Collinear!(material, Nkpt, kpts, Cnk, EVec, iemin, iemax)
    elseif SpinPol == "nc"
        Get_EVec_NonCollinear!(material, Nkpt, kpts, Cnk, EVec, iemin, iemax)
    end
    

    if mode == "all"
        Calc_DosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
        Calc_PDosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
    elseif mode == "dos"
        Calc_DosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
    elseif mode == "pdos"
        Calc_PDosMain(filename, material, Enk, EVec, neg, kmesh, Dos_Erange)
    end
end
