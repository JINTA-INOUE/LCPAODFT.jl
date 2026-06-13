mutable struct Hubbard_U
    SpinPol::String
    Nspin::Int32
    Natom::Int32
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    RMI1::Vector{Vector{Int32}}
    Total_NumOrbs::Vector{Int32}
    Hub_U_occ::String
    Hub_Type::String
    Atom_MaxL_Basis::Vector{Int32}
    Atom_Num_Basis::Vector{Vector{Int32}}
    Hub_U_Basis::Vector{Vector{Vector{Float64}}}
    Hub_U_orbpol::Vector{Bool}
    trans_index::Vector{Vector{Vector{Vector{Int32}}}}
    DM_onsite::Vector{Vector{Vector{Vector{Float64}}}}
    NC_OcpN::Union{Vector{Vector{Vector{Vector{Vector{Float64}}}}}, Nothing}
    v_eff::Union{Vector{Vector{Vector{Vector{Float64}}}}, Nothing}
    NC_v_eff::Union{Vector{Vector{Vector{Vector{Vector{ComplexF64}}}}}, Nothing}
    H_Hub::Vector{Vector{Vector{Vector{Vector{Float64}}}}}
end


"""
    switch is Hub_U_switch
    Hub_U_occ ∈ (`on_site`, `dual`, `full`)      default `dual`
    Hub_Type ∈ (`Dudarev`, `General`)   defalut `Dudarev`
    dc_type ∈ (`sFLL`, `sAMF`, `cFLL`, `cAMF`) (Hub_Type = `General` case only)
"""
function DFT_Hubbard_U(
    SpinPol::String,
    Hub_U_Atom::Vector{Vector{Float64}},
    Hub_U_orbpol::Vector{Bool},
    Hub_U_occ::String,
    Hub_Type::String,
    dc_Type::String,
    pao::Vector{PAO},
    system_grid::System_Grid)


    Hub_U_occ = lowercase(Hub_U_occ)
    if Hub_U_occ ∉ ("onsite", "on_site", "dual")
        println("Now Hubbard Occupation is $Hub_U_occ")
        error("please check Hubbard Occupation.")
    else
        Hub_U_occ = "dual"
    end

    Hub_Type = lowercase(Hub_Type)
    if Hub_Type ∉ ("general", "dudarev")
        println("No Hubbard Type is $Hub_Type")
        error("please check Hubbard Type.")
    else
        Hub_Type = "dudarev"
    end

    dc_Type = lowercase(dc_Type)
    if dc_Type ∉ ("sfll", "samf", "cfll", "camf")
        println("Now DC Type is $dc_Type")
        error("please check DC Type.")
    else
        dc_Type = "sfull"
    end

    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end


    Natom = system_grid.Natom
    atom2spe = system_grid.atom2spe
    FNAN = system_grid.FNAN
    Total_NumOrbs = system_grid.Total_NumOrbs
    natn = system_grid.natn
    RMI = system_grid.RMI
    Total_Hsize = system_grid.Total_Hsize

    # For Noncollinear case
    RMI1 = Vector{Vector{Int32}}(undef, Natom)
    for atom = 1:Natom
        RMI1[atom] = zeros(Int32, FNAN[atom]+1)
        for Rn = 1:FNAN[atom]+1
            RMI1[atom][Rn] = RMI[atom][Rn][1]
        end
    end


    loop2atom = zeros(Int32, Total_Hsize)
    loop2Rn = zeros(Int32, Total_Hsize)
    loop2 = zeros(Int32, Total_Hsize)
    loop2atom = zeros(Int32, Total_Hsize)
    

    
    Atom_MaxL_Basis = zeros(Int32, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        Atom_MaxL_Basis[atom] = pao[spe].Spe_MaxL_Basis
    end


    Atom_Num_Basis = Vector{Vector{Int32}}(undef, Natom)
    for atom = 1:Natom
        spe = atom2spe[atom]
        Atom_Num_Basis[atom] = zeros(Int32, Atom_MaxL_Basis[atom]+1)
        for l = 1:Atom_MaxL_Basis[atom]+1
            Atom_Num_Basis[atom][l] = pao[spe].Spe_Num_Basis[l]
        end
    end
    

    Hub_U_Basis = Set_Hub_U_Basis(Natom, Atom_MaxL_Basis, Atom_Num_Basis, Hub_U_Atom)
    


    trans_index = Vector{Vector{Vector{Vector{Int32}}}}(undef, Natom)
    for atom = 1:Natom
        trans_index[atom] = Vector{Vector{Vector{Int32}}}(undef, Atom_MaxL_Basis[atom]+1)
        for l = 1:Atom_MaxL_Basis[atom]+1
            trans_index[atom][l] = Vector{Vector{Int32}}(undef, Atom_Num_Basis[atom][l])
            for p = 1:Atom_Num_Basis[atom][l]
                trans_index[atom][l][p] = zeros(Int32, 2*l-1)
            end
        end
    end

    
    for atom = 1:Natom
        to1 = 1
        for l = 1:Atom_MaxL_Basis[atom]+1, p = 1:Atom_Num_Basis[atom][l], m = 1:2*l-1
            trans_index[atom][l][p][m] = to1
            to1 += 1
        end
    end



    DM_onsite = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspin)
    for spin = 1:Nspin
        DM_onsite[spin] = Vector{Vector{Vector{Float64}}}(undef, Natom)
        for atom = 1:Natom
            DM_onsite[spin][atom] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
            for ist = 1:Total_NumOrbs[atom]
                DM_onsite[spin][atom][ist] = zeros(Float64, Total_NumOrbs[atom])
            end
        end
    end



    if SpinPol ≠ "nc"
        NC_OcpN = nothing
    else
        NC_OcpN = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, 2)
        for s1 = 1:2
            NC_OcpN[s1] = Vector{Vector{Vector{Vector{Float64}}}}(undef, 2)
            for s2 = 1:2
                NC_OcpN[s1][s2] = Vector{Vector{Vector{Float64}}}(undef, Natom+1)
                for atom = 1:Natom+1
                    if atom == 1
                        NO0 = 1
                    else
                        NO0 = Total_NumOrbs[atom]
                    end

                    NC_OcpN[s1][s2][atom] = Vector{Vector{Float64}}(undef, NO0)
                    for ist = 1:NO0
                        NC_OcpN[s1][s2][atom][ist] = zeros(Float64, NO0)
                    end
                end
            end
        end
    end


    if SpinPol ∈ ("off", "on")
        NC_v_eff = nothing
        v_eff = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspin)
        for spin = 1:Nspin
            v_eff[spin] = Vector{Vector{Vector{Float64}}}(undef, Natom)
            for atom = 1:Natom
                v_eff[spin][atom] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    v_eff[spin][atom][ist] = zeros(Float64, Total_NumOrbs[atom])
                end
            end
        end
    elseif SpinPol == "nc"
        v_eff = nothing
        NC_v_eff = Vector{Vector{Vector{Vector{Vector{ComplexF64}}}}}(undef, 2)
        for spin = 1:2
            NC_v_eff[spin] = Vector{Vector{Vector{Vector{ComplexF64}}}}(undef, 2)
            for jspin = 1:2
                NC_v_eff[spin][jspin] = Vector{Vector{Vector{ComplexF64}}}(undef, Natom)
                for atom = 1:Natom
                    NC_v_eff[spin][jspin][atom] = Vector{Vector{ComplexF64}}(undef, Total_NumOrbs[atom])
                    for ist = 1:Total_NumOrbs[atom]
                        NC_v_eff[spin][jspin][atom][ist] = zeros(ComplexF64, Total_NumOrbs[atom])
                    end
                end
            end
        end
    end


    H_Hub = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, Nspin)
    for spin = 1:Nspin
        H_Hub[spin] = Vector{Vector{Vector{Vector{Float64}}}}(undef, Natom)
        for atom = 1:Natom
            H_Hub[spin][atom] = Vector{Vector{Vector{Float64}}}(undef, FNAN[atom]+1)
            for Rn = 1:FNAN[atom]+1
                H_Hub[spin][atom][Rn] = Vector{Vector{Float64}}(undef, Total_NumOrbs[atom])
                for ist = 1:Total_NumOrbs[atom]
                    H_Hub[spin][atom][Rn][ist] = zeros(Float64, Total_NumOrbs[natn[atom][Rn]])
                end
            end
        end
    end


    return Hubbard_U( SpinPol, Nspin, Natom, FNAN, natn, RMI1, Total_NumOrbs,
                      Hub_U_occ, Hub_Type, 
                      Atom_MaxL_Basis, Atom_Num_Basis, 
                      Hub_U_Basis, Hub_U_orbpol, 
                      trans_index, 
                      DM_onsite, NC_OcpN, 
                      v_eff, NC_v_eff, H_Hub )
end



function Set_Eff_Hub_Pot!(hubbard_u::Union{Hubbard_U}, OLP)

    Hub_U_occ = hubbard_u.Hub_U_occ

    if Hub_U_occ == "onsite"
        error("Hub_U_occ == on site is not support")
        # H_U_onsite!(hubbard_u, OLP)
    elseif Hub_U_occ == "full"
        error("Hub_U_occ == full site is not support")
        # H_U_full!(hubbard_u, OLP)
    elseif Hub_U_occ == "dual"
        H_U_dual!(hubbard_u, OLP)
    else
        error("please check Hub_U_occ")
    end
end


function Occupation_Number_DFT_U!(SCF_iter, hubbard_u::Hubbard_U, DM, OLP)
    
    Hub_U_occ = hubbard_u.Hub_U_occ
    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Hub_U_orbpol = hubbard_u.Hub_U_orbpol
    
    if Hub_U_occ == "onsite"
        error("Hub_U_occ == on site is not support")
        # occupation_onsite!(hubbard_u, DM)
    elseif Hub_U_occ == "full"
        error("Hub_U_occ == full site is not support")
        # occupation_full!(hubbard_u, DM, OLP)
    elseif Hub_U_occ == "dual"
        occupation_dual!(hubbard_u, DM, OLP)
    else
        error("please check Hub_U_occ")
    end


    
    # Orbital polarization
    if SCF_iter < SCF_Enhance
        for atom = 1:Natom
            if Hub_U_orbpol[atom]
                if SpinPol == "off"
                    Induce_Orbital_Polarization!(atom, DM, DM_onsite)
                elseif SpinPol == "on"
                    Induce_Orbital_Polarization_Together!(hubbard_u, atom, DM)
                elseif SpinPol == "nc"
                    Induce_NC_Orbital_Polarization!()
                end
            end
        end
    end
    

    
    if SpinPol ∈ ("off", "on")
        make_v_eff!(hubbard_u)
    elseif SpinPol == "nc"
        error("not support")
        # make_NC_v_eff!(hubbard_u, DM, OLP)
    end
end




# calc DM_onsite
function occupation_onsite!(hubbard_u::Hubbard_U, DM)

    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    Total_NumOrbs = hubbard_u.Total_NumOrbs

    DM_onsite = hubbard_u.DM_onsite

    if SpinPol ∈ ("off", "on")
        # for spin = 1:Nspin, atom = 1:Natom, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
        #     DM_onsite[spin][atom][ist][jst] = DM[spin][atom][1][ist][jst]
        # end

        for spin = 1:Nspin, loop = 1:Total_Hsize
            atom = loop2atom[loop]
            Rn = loop2Rn[loop]
            if Rn == 1
                NO0 = Total_NumOrbs[atom]
                for ist = 1:NO0, jst = 1:NO0
                    DM_onsite[spin][atom][ist][jst] = DM[spin][loop]
                end
            end
        end
    elseif SpinPol == "nc"
        for atom = 1:Natom, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            NC_OcpN[1][1][atom][ist][jst] = ComplexF64(DM[1][atom][1][ist][jst], iDM[1][atom][1][ist][jst])
            NC_OcpN[2][2][atom][ist][jst] = ComplexF64(DM[2][atom][1][ist][jst], iDM[2][atom][1][ist][jst])
            NC_OcpN[1][2][atom][ist][jst] = ComplexF64(DM[3][atom][1][ist][jst],  DM[4][atom][1][ist][jst])
            NC_OcpN[2][1][atom][ist][jst] = ComplexF64(DM[3][atom][1][jst][ist], -DM[4][atom][1][jst][ist])
        end
    end

    hubbard_u.DM_onsite = DM_onsite
end




function occupation_full!(hubbard_u::Hubbard_U, DM, OLP)

    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    Total_NumOrbs = hubbard_u.Total_NumOrbs
    RMI1 = hubbard_u.RMI1

    DM_onsite = hubbard_u.DM_onsite

    if SpinPol ∈ ("off", "on")
        for spin = 1:Nspin, atom = 1:Natom, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]

            Sum = 0.0
            for Rn = 1:FNAN[atom]+1
                katom = natn[atom][Rn]
                for Rm = 1:FNAN[atom]+1
                    latom = natn[atom][Rm]
                    kl = RMI1[atom][Rn][Rm]
                    if kl >= 0
                        for kst = 1:Total_NumOrbs[katom], lst = 1:Total_NumOrbs[latom]
                            Sum += DM[spin][katom][kl][kst][lst]*OLP[atom][Rn][ist][kst]*OLP[atom][Rm][jst][lst]
                        end
                    end
                end
            end

            DM_onsite[spin][atom][ist][jst] = Sum
        end
    elseif SpinPol == "nc"
        error("not support yet")
    end


    hubbard_u.DM_onsite = DM_onsite
end


function occupation_dual!(hubbard_u::Hubbard_U, DM, OLP)

    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    FNAN = hubbard_u.FNAN
    natn = hubbard_u.natn
    Total_NumOrbs = hubbard_u.Total_NumOrbs


    if SpinPol ∈ ("off", "on")

        DM_onsite = hubbard_u.DM_onsite

        for spin = 1:Nspin, atom = 1:Natom, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            Sum = 0.0
            for Rn = 1:FNAN[atom]+1
                jatom = natn[atom][Rn]
                for kst = 1:Total_NumOrbs[jatom]
                    Sum += DM[spin][atom][Rn][jst][kst]*OLP[atom][Rn][ist][kst] + DM[spin][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst]
                end
            end

            DM_onsite[spin][atom][ist][jst] = 0.5*Sum
        end


        hubbard_u.DM_onsite = DM_onsite

    elseif SpinPol == "nc"

        error("not support occupation_dual for NonCollinear case")

        RMI1 = hubbard_u.RMI1
        NC_OcpN = hubbard_u.NC_OcpN

        for atom = 1:Natom, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[atom]
            ReOcn00 = 0.0
            ReOcn11 = 0.0
            ReOcn01 = 0.0
      
            ImOcn00 = 0.0
            ImOcn11 = 0.0
            ImOcn01 = 0.0
            for Rn = 1:FNAN[atom]+1
                jatom = natn[atom][Rn]
                kl = RMI1[atom][Rn]
                for kst = 1:Total_NumOrbs[jatom]
                    ReOcn00 += 0.5*( DM[1][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] +  DM[1][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                    ReOcn11 += 0.5*( DM[2][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] +  DM[2][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                    ReOcn11 += 0.5*( DM[3][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] +  DM[2][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                    ImOcn00 += 0.5*(iDM[1][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] + iDM[1][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                    ImOcn11 += 0.5*(iDM[2][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] + iDM[2][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                    ImOcn01 += 0.5*( DM[4][atom][Rn][ist][kst]*OLP[atom][Rn][jst][kst] +  DM[4][atom][kl][kst][jst]*OLP[atom][Rn][ist][kst])
                end
            end

            NC_OcpN[1][1][atom][ist][jst] = ComplexF64(ReOcn00, ImOcn00)
            NC_OcpN[2][2][atom][ist][jst] = ComplexF64(ReOcn11, ImOcn11)
            NC_OcpN[1][2][atom][ist][jst] = ComplexF64(ReOcn01, ImOcn01)
            NC_OcpN[2][1][atom][ist][jst] = ComplexF64(ReOcn01, -ImOcn01)
        end

        hubbard_u.NC_OcpN = NC_OcpN
    end
end


function make_v_eff!(hubbard_u::Hubbard_U)

    Atom_MaxL_Basis = hubbard_u.Atom_MaxL_Basis
    Atom_Num_Basis = hubbard_u.Atom_Num_Basis
    Hub_U_Basis = hubbard_u.Hub_U_Basis
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    DM_onsite = hubbard_u.DM_onsite
    v_eff = hubbard_u.v_eff

    # Hub_Type = "Dudarev" case
    for spin = 1:Nspin, atom = 1:Natom
        ist = 1
        for l1 = 0:Atom_MaxL_Basis[atom], p1 = 1:Atom_Num_Basis[atom][l1+1], m1 = 1:2*l1+1
            jst = 1
            for l2 = 0:Atom_MaxL_Basis[atom], p2 = 1:Atom_Num_Basis[atom][l2+1], m2 = 1:2*l2+1
                if l1 == l2 && p1 == p2
                    Uvalues = Hub_U_Basis[atom][l1+1][p1]
                else
                    Uvalues = 0.0
                end

                if ist == jst
                    v_eff[spin][atom][ist][jst] = Uvalues*(0.5 - DM_onsite[spin][atom][ist][jst])
                else
                    v_eff[spin][atom][ist][jst] = Uvalues*(0.0 - DM_onsite[spin][atom][ist][jst])
                end
                jst += 1
            end
            ist += 1
        end
    end

    hubbard_u.v_eff = v_eff
end


function make_v_eff_NC(hubbard_u::Hubbard_U)

    Atom_MaxL_Basis = hubbard_u.Atom_MaxL_Basis
    Atom_Num_Basis = hubbard_u.Atom_Num_Basis
    Hub_U_Basis = hubbard_u.Hub_U_Basis
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    DM_onsite = hubbard_u.DM_onsite
    NC_v_eff = hubbard_u.NC_v_eff

    # Hub_Type = "Dudarev" case
    for spin = 1:Nspin, atom = 1:Natom
        ist = 1
        for l1 = 0:Atom_MaxL_Basis[atom], p1 = 1:Atom_Num_Basis[atom][l1+1], m1 = 1:2*l1+1
            jst = 1
            for l2 = 0:Atom_MaxL_Basis[atom], p2 = 1:Atom_Num_Basis[atom][l2+1], m2 = 1:2*l2+1
                if l1 == l2 && p1 == p2
                    Uvalues = Hub_U_Basis[atom][l1+1][p1]
                else
                    Uvalues = 0.0
                end

                if ist == jst
                    NC_v_eff[spin][atom][ist][jst] = Uvalues*(0.5 - DM_onsite[spin][atom][ist][jst])
                else
                    NC_v_eff[spin][atom][ist][jst] = Uvalues*(0.0 - DM_onsite[spin][atom][ist][jst])
                end
                jst += 1
            end
            ist += 1
        end
    end

    hubbard_u.NC_v_eff = NC_v_eff
end


function Induce_Orbital_Polarization!(atom, DM, DM_onsite)

    Ncut = 0.3
    Ns = 20
    A = zeros(Float64, Ns, Ns)

    mul1 = 1
    for l1 = 2:Spe_MaxL_Basis
        for spin = 1:Nspin
            for m1 = 1:2*l1+1, m2 = 1:2*l1+1
                to1 = trans_index[atom][l1][mul1][m1]
                to2 = trans_index[atom][l1][mul1][m2]

                A[m1+1,m2+1] = DM[spin][atom][1][to1][to2]
            end

            Sum1 = 0.0
            for m1 = 1:2*l1+1
                Sum1 += A[m1+1,m1+1]
            end

            if Sum1 > Ncut

                val, vec = eigen(A)

                if Hub_U_orbpol[atom] == 1

                elseif Hub_U_orbpol[atom] == 2

                end

                for m1 = 1:2*l1+1, m2 = 1:2*l1+1
                    Sum1 = 0.0
                    for m3 = 1:2*l1+1
                        Sum1 += 1
                    end

                    to1 = trans_index[atom][l1][mul1][m1]
                    to2 = trans_index[atom][l1][mul1][m2]

                    DM_onsite[spin][atom][to1][to2] = Sum1
                end
            end
        end
    end                 
end


function Induce_Orbital_Polarization_Together!(hubbard_u::Hubbard_U, atom, DM)

    Atom_MaxL_Basis = hubbard_u.Atom_MaxL_Basis
    trans_index = hubbard_u.trans_index[atom]
    DM_onsite = hubbard_u.DM_onsite

    Ncut = 0.3
    Ns = 4*4*2 + 1
    A = zeros(Float64, Ns, Ns)


    # @show atom

    mul1 = 1
    for l1 = 2:Atom_MaxL_Basis[atom]

        k = 2*l1 + 1

        for m1 = 1:2*l1+1, m2 = 1:2*l1+1
            to1 = trans_index[l1+1][mul1][m1]
            to2 = trans_index[l1+1][mul1][m2]

            A[m1+1,m2+1]     = DM[1][atom][1][to1][to2] + rand()*1e-13
            A[m1+k+1,m2+k+1] = DM[2][atom][1][to1][to2] + rand()*1e-13
            A[m1+1,m2+k+1]   = rand()*1e-13
            A[m1+k+1,m2+1]   = rand()*1e-13
        end

        sum1 = 0.0
        for m1 = 1:2*k
            sum1 += A[m1+1,m1+1]
        end

        # OK
        # @show l1, sum1

        if sum1 > Ncut

            values, vectors = eigen(A[1:2*k+1,1:2*k+1])

            #OK 
            for m1 = 1:2*k
                @show m1, values[m1+1]
            end

            toccpn = 0.0
            for m1 = 1:2*k
                toccpn += values[m1+1]
            end
            fill!(values, 0.0)
            # OK
            # @show toccpn


            for m1 = 1:2*k
                values[m1+1] = 0.0
            end


            m0 = 4*l1 + 2 - trunc(Int64, toccpn)
            # @show m0, trunc(Int64, toccpn)
            if m0 < 0
                m0 = 0
            end

            for m1 = 4*l1+2:-1:m0+1
                values[m1+1] = 1.0
            end

            if 0 <= 4*l1 + 1 - trunc(Int64, toccpn)
                values[4*l1+2-trunc(Int64, toccpn)] = toccpn-trunc(Int64, toccpn)
            end

            for m1 = 1:2*k
                # @printf("Col2 m1=%2d %15.12f\n",m1,values[m1+1])
            end

            for m1 = 1:2*k, m2 = 1:2*k
                # @printf("Col2 m1=%2d m2=%d %15.12f\n",m1-1,m2-1,vectors[m1+1,m2+1]);
            end


            for m1 = 1:2*k, m2 = 1:2*k
                sum1 = 0.0
                for m3 = 1:2*k
                    sum1 += vectors[m1+1,m3+1]*values[m3+1]*vectors[m2+1,m3+1]
                end

                # TODO: check
                # @show m1-1, m2-1, sum1

                to1 = trans_index[l1+1][mul1][(m1-1)%k+1]
                to2 = trans_index[l1+1][mul1][(m2-1)%k+1]

                if div(m1,k)==0 && div(m2,k)==0
                    DM_onsite[1][atom][to1][to2] = sum1
                elseif div(m1,k)==1 && div(m2,k)==1
                    DM_onsite[2][atom][to1][to2] = sum1
                end
            end
        end
    end

    hubbard_u.DM_onsite = DM_onsite
end


# NonCollinear spin polarization case
function Induce_NC_Orbital_Polarization!()

end


function H_U_onsite!(hubbard_u::Hubbard_U)

    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    natn = hubbard_u.natn
    Total_NumOrbs = hubbard_u.Total_NumOrbs
    H_Hub = hubbard_u.H_Hub

    if SpinPol ∈ ("off", "on")

        v_eff = hubbard_u.v_eff
        for spin = 1:Nspin, atom = 1:Natom
            jatom = natn[atom][begin]
            for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]
                H_Hub[spin][atom][1][ist][jst] = v_eff[spin][atom][ist][jst]
            end
        end
    elseif SpinPol == "nc"

        NC_v_eff = hubbard_u.NC_v_eff
        for atom = 1:Natom
            jatom = natn[atom][begin]
            for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]
                H_Hub[1][atom][1][ist][jst] = NC_v_eff[1][1][atom][ist][jst]
                H_Hub[2][atom][1][ist][jst] = NC_v_eff[2][2][atom][ist][jst]
                H_Hub[3][atom][1][ist][jst] = NC_v_eff[1][2][atom][ist][jst]
            end
        end
    end

    hubbard_u.H_Hub = H_Hub
end


function H_U_full!(hubbard_u::Hubbard_U, OLP)
    
    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    FNAN = hubbard_u.FNAN
    natn = hubbard_u.natn
    Total_NumOrbs = hubbard_u.Total_NumOrbs
    RMI1 = hubbard_u.RMI1
    H_Hub = hubbard_u.H_Hub


    for spin = 1:Nspin, atom = 1:Natom, Rn = 1:FNAN[atom]+1, ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]
        H_Hub[spin][atom][Rn][ist][jst] = 0.0
    end


    if SpinPol ∈ ("off", "on")

        v_eff = hubbard_u.v_eff
        for spin = 1:Nspin, atom = 1:Natom, Rn = 1:FNAN[atom]+1, Rm = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            katom = natn[atom][Rm]
            kl = RMI1[atom][Rn][Rm]
            if kl >= 0
                for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]
                    tmp = 0.0
                    for kst = 1:Total_NumOrbs[katom], lst = 1:Total_NumOrbs[katom]
                        tmp += OLP[atom][Rm][ist][kst]*v_eff[spin][katom][kst][lst]*OLP[jatom][kl+1][jst][lst]
                    end

                    H_Hub[spin][atom][Rn][ist][jst] += tmp
                end
            end
        end

    elseif SpinPol == "nc"

        NC_v_eff = hubbard_u.NC_v_eff
        error("not support yet.")
    end

    hubbard_u.H_Hub = H_Hub
end


function H_U_dual!(hubbard_u::Hubbard_U, OLP)

    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Nspin = hubbard_u.Nspin
    FNAN = hubbard_u.FNAN
    natn = hubbard_u.natn
    Total_NumOrbs = hubbard_u.Total_NumOrbs
    H_Hub = hubbard_u.H_Hub

    if SpinPol ∈ ("off", "on")

        v_eff = hubbard_u.v_eff
        for spin = 1:Nspin, atom = 1:Natom, Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]
                tmp = 0.0
                for kst = 1:Total_NumOrbs[atom]
                    tmp += v_eff[spin][atom][ist][kst]*OLP[atom][Rn][kst][jst]
                end

                for kst = 1:Total_NumOrbs[jatom]
                    tmp += v_eff[spin][jatom][kst][jst]*OLP[atom][Rn][ist][kst]
                end
                H_Hub[spin][atom][Rn][ist][jst] = 0.5*tmp
            end
        end

    elseif SpinPol == "nc"

        NC_v_eff = hubbard_u.NC_v_eff
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            
            for ist = 1:Total_NumOrbs[atom], jst = 1:Total_NumOrbs[jatom]

                sum00 = ComplexF64(0.0, 0.0)
                sum11 = ComplexF64(0.0, 0.0)
                sum01 = ComplexF64(0.0, 0.0)

                for kst = 1:Total_NumOrbs[atom]
                    sum00 += NC_v_eff[1][1][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                    sum11 += NC_v_eff[2][2][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                    sum01 += NC_v_eff[1][2][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                end

                for kst = 1:Total_NumOrbs[jatom]
                    sum00 += NC_v_eff[1][1][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                    sum11 += NC_v_eff[2][2][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                    sum01 += NC_v_eff[1][2][atom][ist][kst] * OLP[atom][Rn][kst][jst]
                end
                
                H_Hub[1][atom][Rn][ist][jst] = 0.5*real(sum00)
                H_Hub[2][atom][Rn][ist][jst] = 0.5*real(sum11)
                H_Hub[3][atom][Rn][ist][jst] = 0.5*real(sum01)

                iHNL[1][atom][Rn][ist][jst] = iHNL[1][atom][Rn][ist][jst] + 0.5*imag(sum00)
                iHNL[2][atom][Rn][ist][jst] = iHNL[2][atom][Rn][ist][jst] + 0.5*imag(sum11)
                iHNL[3][atom][Rn][ist][jst] = iHNL[3][atom][Rn][ist][jst] + 0.5*imag(sum01)
            end
        end
    end

    hubbard_u.H_Hub = H_Hub
end



function Calc_EHub(hubbard_u::Hubbard_U)


    SpinPol = hubbard_u.SpinPol
    Natom = hubbard_u.Natom
    Atom_MaxL_Basis = hubbard_u.Atom_MaxL_Basis
    Atom_Num_Basis = hubbard_u.Atom_Num_Basis
    Hub_U_Basis = hubbard_u.Hub_U_Basis
    DM_onsite = hubbard_u.DM_onsite
    spinsize = ifelse(SpinPol=="off", 1, 2)

    EHub = 0.0

    # Dudarev form
    if SpinPol ∈ ("off", "on")
        for atom = 1:Natom, spin = 1:spinsize
            tot1 = 1
            for l1 = 1:Atom_MaxL_Basis[atom]+1, p1 = 1:Atom_Num_Basis[atom][l1], m1 = 1:2*l1-1
                EHub += 0.5 * Hub_U_Basis[atom][l1][p1] * DM_onsite[spin][atom][tot1][tot1]
                tot1 += 1
            end

            tot1 = 1
            for l1 = 1:Atom_MaxL_Basis[atom]+1, p1 = 1:Atom_Num_Basis[atom][l1], m1 = 1:2*l1-1
                tot2 = 1
                Sum = 0.0
                for l2 = 1:Atom_MaxL_Basis[atom]+1, p2 = 1:Atom_Num_Basis[atom][l2], m2 = 1:2*l2-1
                    if l1 == l2 && p1 == p2
                        Sum -= 0.5 * Hub_U_Basis[atom][l1][p1] * DM_onsite[spin][atom][tot1][tot2] * DM_onsite[spin][atom][tot2][tot1]
                    end
                    tot2 += 1
                end
                tot1 += 1
                EHub += Sum
            end
        end
    elseif SpinPol == "nc"
        # TODO
        for atom = 1:Natom, spin = 1:spinsize
            tot1 = 1
            for l1 = 1:Atom_MaxL_Basis[atom]+1, p1 = 1:Atom_Num_Basis[atom][l1], m1 = 1:2*l1-1
                EHub += 1
                tot1 += 1
            end

            tot1 = 1
            for l1 = 1:Atom_MaxL_Basis[atom]+1, p1 = 1:Atom_Num_Basis[atom][l1], m1 = 1:2*l1-1
                tot2 = 1
                Sum = 0.0
                for l2 = 1:Atom_MaxL_Basis[atom]+1, p2 = 1:Atom_Num_Basis[atom][l2], m2 = 1:2*l2-1
                    if l1 == l2 && p1 == p2
                        Sum -= 1
                    end
                    tot2 += 1
                end
                tot1 += 1
                EHub += Sum
            end
        end
    else
        error("please check SpinPolarization.")
    end


    if SpinPol == "off"
        EHub = 2.0*EHub
    end


    return EHub
end
