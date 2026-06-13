function Calc_Vnk!(boltz_setup::Boltz_Setup, kpoints::KPoints, Enk, Cnk)

    material = boltz_setup.material
    Latvecs = material.Latvecs
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    kmesh = boltz_setup.kmesh
    kmesh1, kmesh2, kmesh3 = kmesh
    ChemP = material.ChemP
    NCell = material.NCell
    cell_list_ijk = material.cell_list_ijk
    HmnR = material.HmnR

    if SpinPol ∈ ("off", "nc")
        spinsize = 1
    elseif SpinPol == "on"
        spinsize = 2
    else
        error("please check SpinPol")
    end


    Nkpt = kpoints.Nkpt
    kpts = kpoints.MPI_kpts
    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:kmesh1, jk = 1:kmesh2, kk = 1:kmesh3
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    cartesian_cell = zeros(Float64, 3, NCell)
    for cell = 1:NCell
        cartesian_cell[1,cell] = Latvecs[1,1]*cell_list_ijk[cell][1] + Latvecs[2,1]*cell_list_ijk[cell][2] + Latvecs[3,1]*cell_list_ijk[cell][3]
        cartesian_cell[2,cell] = Latvecs[1,2]*cell_list_ijk[cell][1] + Latvecs[2,2]*cell_list_ijk[cell][2] + Latvecs[3,2]*cell_list_ijk[cell][3]
        cartesian_cell[3,cell] = Latvecs[1,3]*cell_list_ijk[cell][1] + Latvecs[2,3]*cell_list_ijk[cell][2] + Latvecs[3,3]*cell_list_ijk[cell][3]
        cartesian_cell[1,cell] = cartesian_cell[1,cell]/Ang_to_bohr
        cartesian_cell[2,cell] = cartesian_cell[2,cell]/Ang_to_bohr
        cartesian_cell[3,cell] = cartesian_cell[3,cell]/Ang_to_bohr
    end



    for spin = 1:spinsize, cell = 1:NCell, ist = 1:Nwann, jst = 1:Nwann
        HmnR[jst,ist,cell,spin] = HmnR[jst,ist,cell,spin]*eV2Hartree
    end


    for spin = 1:spinsize, cell = 1:NCell, ist = 1:Nwann
        if cell_list_ijk[cell] == [0,0,0]
            HmnR[ist,ist,cell,spin] = HmnR[ist,ist,cell,spin] - ChemP*eV2Hartree
        end
    end
    

    Vnk = Vector{Vector{Vector{Vector{Vector{Vector{Float64}}}}}}(undef, spinsize)
    for spin = 1:spinsize
        Vnk[spin] = Vector{Vector{Vector{Vector{Vector{Float64}}}}}(undef, kmesh1)
        for ik = 1:kmesh1
            Vnk[spin][ik] = Vector{Vector{Vector{Vector{Float64}}}}(undef, kmesh2)
            for jk = 1:kmesh2
                Vnk[spin][ik][jk] = Vector{Vector{Vector{Float64}}}(undef, kmesh3)
                for kk = 1:kmesh3
                    Vnk[spin][ik][jk][kk] = Vector{Vector{Float64}}(undef, 3)
                    for xyz = 1:3
                        Vnk[spin][ik][jk][kk][xyz] = zeros(Float64, Nwann)
                    end
                end
            end
        end
    end

    


    Hmnk_x = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_y = zeros(ComplexF64, Nwann, Nwann)
    Hmnk_z = zeros(ComplexF64, Nwann, Nwann)
    Htemp1 = zeros(ComplexF64, Nwann, Nwann)
    Htemp2 = zeros(ComplexF64, Nwann, Nwann)

    for spin = 1:spinsize, k = 1:Nkpt

        fill!(Hmnk_x, 0.0)
        fill!(Hmnk_y, 0.0)
        fill!(Hmnk_z, 0.0)
            
        for cell = 1:NCell
            kRn = dot(kpts[k], cell_list_ijk[cell])
            phase = cispi(2*kRn)

            tmpx = -im*cartesian_cell[1,cell]
            tmpy = -im*cartesian_cell[2,cell]
            tmpz = -im*cartesian_cell[3,cell]

            @inbounds for ist = 1:Nwann, jst = 1:Nwann
                Hmnk_x[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpx
                Hmnk_y[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpy
                Hmnk_z[ist,jst] += HmnR[jst,ist,cell,spin] * phase * tmpz
            end
        end

        ik, jk, kk = kindex[k,:]

        mul!(Htemp1, Cnk[spin][k]', Hmnk_x)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][1] = real(diag(Htemp2))
        mul!(Htemp1, Cnk[spin][k]', Hmnk_y)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][2] = real(diag(Htemp2))
        mul!(Htemp1, Cnk[spin][k]', Hmnk_z)
        mul!(Htemp2, Htemp1, Cnk[spin][k])
        Vnk[spin][ik][jk][kk][3] = real(diag(Htemp2))
    end


    # change unit (Hartree -> eV)
    for spin = 1:spinsize, k = 1:Nkpt
        ik, jk, kk = kindex[k,:]
        @inbounds for μ = 1:Nwann
            Enk[spin][ik][jk][kk][μ] = Enk[spin][ik][jk][kk][μ]*eV2Hartree - ChemP*eV2Hartree
        end
    end


    return Vnk
end