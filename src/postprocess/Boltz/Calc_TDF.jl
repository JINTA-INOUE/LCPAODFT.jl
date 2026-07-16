@timeit timer "Calc_TDF" function Calc_TDF(boltz_setup::Boltz_Setup, Enk, Vnk)
    
    material = boltz_setup.material
    SpinPol = material.SpinPol
    # spinsize = ifelse(SpinPol=="on", 2, 1)
    # Nwann = material.Ngsize
    filename = boltz_setup.filename
    # filepath = boltz_setup.filepath
    # kmesh = boltz_setup.kmesh
    # decomp = boltz_setup.decomp
    plane_type = boltz_setup.plane_type
    mat_type = boltz_setup.mat_type
    Write_TDF = boltz_setup.Write_TDF
    tau = boltz_setup.tau
    # TDF_Erange = boltz_setup.TDF_Erange
    # TDF_dE = boltz_setup.TDF_dE

    if plane_type
        TDF_Energy, TDF = Calc_TDF_3element(boltz_setup, Enk, Vnk)
    else
        TDF_Energy, TDF = Calc_TDF_6element(boltz_setup, Enk, Vnk)
    end

    if Write_TDF
        if plane_type
            Write_TDF_3elements(filename, mat_type, SpinPol, TDF_Energy, tau, TDF)
        else
            Write_TDF_6elements(filename, mat_type, SpinPol, TDF_Energy, tau, TDF)
        end
    end

    #=
    if Write_TDF
        println("Write $(filename).TDF.jld2")
        jldopen("$(filename).TDF.jld2", "w") do file
            file["Dates"] = now()
            file["filepath"] = filepath
            file["SpinPol"] = SpinPol
            file["spinsize"] = spinsize
            file["Nwann"] = Nwann
            file["kmesh"] = kmesh
            file["decomp"] = decomp
            file["plane_type"] = plane_type
            file["tau"] = tau
            file["TDF_Erange"] = TDF_Erange
            file["TDF_dE"] = TDF_dE
            file["TDF_Energy"] = TDF_Energy
            file["TDF"] = TDF
        end
    end=#


    return TDF_Energy, TDF
end


function Calc_TDF_3element(boltz_setup::Boltz_Setup, Enk, Vnk)

    material = boltz_setup.material
    Latvecs = material.Latvecs
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nstate = Set_Nstate(material)
    kmesh = boltz_setup.kmesh
    knum_i, knum_j, knum_k = kmesh
    Nkpt = prod(kmesh)
    tau = boltz_setup.tau
    TDF_Erange = boltz_setup.TDF_Erange
    TDF_dE = boltz_setup.TDF_dE
    cell_volume_AU = abs(det(Latvecs))
    cell_volume_Ang = cell_volume_AU/Ang_to_bohr^3
    plane_type = boltz_setup.plane_type

    if !plane_type
        error("please check plane_type")
    end



    TDF_Emin, TDF_Emax = TDF_Erange
    TDF_EneNum = floor(Int, (TDF_Emax-TDF_Emin)/TDF_dE)  # length of TDF

    TDF_Energy = zeros(Float64, TDF_EneNum)
    TDF = zeros(Float64, TDF_EneNum, spinsize, 6)
    for ie = 1:TDF_EneNum
        TDF_Energy[ie] = TDF_Emin + (TDF_Emax-TDF_Emin)*(ie-1)/(TDF_EneNum-1)
    end


    Nkpt = prod(kmesh)
    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:knum_i, jk = 1:knum_j, kk = 1:knum_k
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    cell_e = zeros(Float64, 8)
    cell_axx = zeros(Float64, 8)
    cell_axy = zeros(Float64, 8)
    cell_ayy = zeros(Float64, 8)
    tetra_exx = zeros(Float64, 4)
    tetra_exy = zeros(Float64, 4)
    tetra_eyy = zeros(Float64, 4)
    tetra_axx = zeros(Float64, 4)
    tetra_axy = zeros(Float64, 4)
    tetra_ayy = zeros(Float64, 4)
    tetra_id = [1 2 3 6; 2 3 4 6; 3 4 6 8; 1 3 5 6; 3 5 6 7; 3 6 7 8]


    # tetrahedron
    for spin = 1:spinsize, ieg = 1:Nstate
        for ik = 1:Nkpt
            i = kindex[ik,1]-1
            j = kindex[ik,2]-1
            k = kindex[ik,3]-1

            for i_in = 0:1, j_in = 0:1, k_in = 0:1
                ii = mod(i+i_in,knum_i)+1
                jj = mod(j+j_in,knum_j)+1
                kk = mod(k+k_in,knum_k)+1
                cell_e[4*i_in+2*j_in+k_in+1] = Enk[spin][ii][jj][kk][ieg]

                vvalx = Vnk[spin][ii][jj][kk][1][ieg]
                vvaly = Vnk[spin][ii][jj][kk][2][ieg]
                cell_axx[4*i_in+2*j_in+k_in+1] = vvalx^2
                cell_axy[4*i_in+2*j_in+k_in+1] = vvalx*vvaly
                cell_ayy[4*i_in+2*j_in+k_in+1] = vvaly^2
            end

            for itetra = 1:6
                for ic = 1:4
                    tetra_exx[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_exy[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_eyy[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_axx[ic] = cell_axx[tetra_id[itetra,ic]]
                    tetra_axy[ic] = cell_axy[tetra_id[itetra,ic]]
                    tetra_ayy[ic] = cell_ayy[tetra_id[itetra,ic]]
                end
                

                OrderE!(tetra_exx, tetra_axx, 4)
                OrderE!(tetra_exy, tetra_axy, 4)
                OrderE!(tetra_eyy, tetra_ayy, 4)
                
                xx = (tetra_exx[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_xx = trunc(Int, xx)
                xx = (tetra_exx[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_xx = trunc(Int, xx)

                xy = (tetra_exy[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_xy = trunc(Int, xy)
                xy = (tetra_exy[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_xy = trunc(Int, xy)

                yy = (tetra_eyy[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_yy = trunc(Int, yy)
                yy = (tetra_eyy[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_yy = trunc(Int, yy)

                # xx elements
                if iemin_xx < 0
                    iemin_xx = 0
                end
                if iemax_xx >= TDF_EneNum
                    iemax_xx = TDF_EneNum - 1
                end
                if (0 <= iemin_xx < TDF_EneNum) && (0 <= iemax_xx < TDF_EneNum)
                    for ie = iemin_xx:iemax_xx
                        resultxx = ATM_Spectrum(tetra_exx, tetra_axx, TDF_Energy[ie+1])
                        TDF[ie+1,spin,1] += resultxx
                    end
                end

                # xy elements
                if iemin_xy < 0
                    iemin_xy = 0
                end
                if iemax_xy >= TDF_EneNum
                    iemax_xy = TDF_EneNum - 1
                end
                if (0 <= iemin_xy < TDF_EneNum) && (0 <= iemax_xy < TDF_EneNum)
                    for ie = iemin_xy:iemax_xy
                        resultxy = ATM_Spectrum(tetra_exy, tetra_axy, TDF_Energy[ie+1])
                        TDF[ie+1,spin,2] += resultxy
                    end
                end

                # yy elements
                if iemin_yy < 0
                    iemin_yy = 0
                end
                if iemax_yy >= TDF_EneNum
                    iemax_yy = TDF_EneNum - 1
                end
                if (0 <= iemin_yy < TDF_EneNum) && (0 <= iemax_yy < TDF_EneNum)
                    for ie = iemin_yy:iemax_yy
                        resultyy = ATM_Spectrum(tetra_eyy, tetra_ayy, TDF_Energy[ie+1])
                        TDF[ie+1,spin,3] += resultyy
                    end
                end
            end
        end
    end

    
    factor = 1/knum_i/knum_j/knum_k/6

    for spin = 1:spinsize, ie = 1:TDF_EneNum
        TDF[ie,spin,1] = TDF[ie,spin,1]*factor*tau/cell_volume_Ang   # xx elements
        TDF[ie,spin,2] = TDF[ie,spin,2]*factor*tau/cell_volume_Ang   # xy elements
        TDF[ie,spin,3] = TDF[ie,spin,3]*factor*tau/cell_volume_Ang   # yy elements
    end


    return TDF_Energy, TDF
end


function Calc_TDF_6element(boltz_setup::Boltz_Setup, Enk, Vnk)

    material = boltz_setup.material
    Latvecs = material.Latvecs
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    Nstate = Set_Nstate(material)
    kmesh = boltz_setup.kmesh
    knum_i, knum_j, knum_k = kmesh
    Nkpt = prod(kmesh)
    tau = boltz_setup.tau
    TDF_Erange = boltz_setup.TDF_Erange
    TDF_dE = boltz_setup.TDF_dE
    cell_volume_AU = abs(det(Latvecs))
    cell_volume_Ang = cell_volume_AU/Ang_to_bohr^3


    TDF_Emin, TDF_Emax = TDF_Erange
    TDF_EneNum = floor(Int, (TDF_Emax-TDF_Emin)/TDF_dE)  # length of TDF

    TDF_Energy = zeros(Float64, TDF_EneNum)
    TDF = zeros(Float64, TDF_EneNum, spinsize, 6)
    for ie = 1:TDF_EneNum
        TDF_Energy[ie] = TDF_Emin + (TDF_Emax-TDF_Emin)*(ie-1)/(TDF_EneNum-1)
    end


    Nkpt = prod(kmesh)
    kindex = zeros(Int64, Nkpt, 3)
    kp = 0
    for ik = 1:knum_i, jk = 1:knum_j, kk = 1:knum_k
        kp += 1
        kindex[kp,1] = ik
        kindex[kp,2] = jk
        kindex[kp,3] = kk
    end


    cell_e = zeros(Float64, 8)
    cell_axx = zeros(Float64, 8)
    cell_axy = zeros(Float64, 8)
    cell_axz = zeros(Float64, 8)
    cell_ayy = zeros(Float64, 8)
    cell_ayz = zeros(Float64, 8)
    cell_azz = zeros(Float64, 8)
    tetra_exx = zeros(Float64, 4)
    tetra_exy = zeros(Float64, 4)
    tetra_exz = zeros(Float64, 4)
    tetra_eyy = zeros(Float64, 4)
    tetra_eyz = zeros(Float64, 4)
    tetra_ezz = zeros(Float64, 4)
    tetra_axx = zeros(Float64, 4)
    tetra_axy = zeros(Float64, 4)
    tetra_axz = zeros(Float64, 4)
    tetra_ayy = zeros(Float64, 4)
    tetra_ayz = zeros(Float64, 4)
    tetra_azz = zeros(Float64, 4)
    tetra_id = [1 2 3 6; 2 3 4 6; 3 4 6 8; 1 3 5 6; 3 5 6 7; 3 6 7 8]


    # tetrahedron
    for spin = 1:spinsize, ieg = 1:Nstate
        for ik = 1:Nkpt
            i = kindex[ik,1]-1
            j = kindex[ik,2]-1
            k = kindex[ik,3]-1

            for i_in = 0:1, j_in = 0:1, k_in = 0:1
                ii = mod(i+i_in,knum_i)+1
                jj = mod(j+j_in,knum_j)+1
                kk = mod(k+k_in,knum_k)+1
                cell_e[4*i_in+2*j_in+k_in+1] = Enk[spin][ii][jj][kk][ieg]

                vvalx = Vnk[spin][ii][jj][kk][1][ieg]
                vvaly = Vnk[spin][ii][jj][kk][2][ieg]
                vvalz = Vnk[spin][ii][jj][kk][3][ieg]
                cell_axx[4*i_in+2*j_in+k_in+1] = vvalx^2
                cell_axy[4*i_in+2*j_in+k_in+1] = vvalx*vvaly
                cell_axz[4*i_in+2*j_in+k_in+1] = vvalx*vvalz
                cell_ayy[4*i_in+2*j_in+k_in+1] = vvaly^2
                cell_ayz[4*i_in+2*j_in+k_in+1] = vvaly*vvalz
                cell_azz[4*i_in+2*j_in+k_in+1] = vvalz^2
            end

            for itetra = 1:6
                for ic = 1:4
                    tetra_exx[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_exy[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_exz[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_eyy[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_eyz[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_ezz[ic] = cell_e[tetra_id[itetra,ic]]
                    tetra_axx[ic] = cell_axx[tetra_id[itetra,ic]]
                    tetra_axy[ic] = cell_axy[tetra_id[itetra,ic]]
                    tetra_axz[ic] = cell_axz[tetra_id[itetra,ic]]
                    tetra_ayy[ic] = cell_ayy[tetra_id[itetra,ic]]
                    tetra_ayz[ic] = cell_ayz[tetra_id[itetra,ic]]
                    tetra_azz[ic] = cell_azz[tetra_id[itetra,ic]]
                end
                

                OrderE!(tetra_exx, tetra_axx, 4)
                OrderE!(tetra_exy, tetra_axy, 4)
                OrderE!(tetra_exz, tetra_axz, 4)
                OrderE!(tetra_eyy, tetra_ayy, 4)
                OrderE!(tetra_eyz, tetra_ayz, 4)
                OrderE!(tetra_ezz, tetra_azz, 4)
                
                xx = (tetra_exx[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_xx = trunc(Int, xx)
                xx = (tetra_exx[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_xx = trunc(Int, xx)

                xy = (tetra_exy[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_xy = trunc(Int, xy)
                xy = (tetra_exy[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_xy = trunc(Int, xy)

                xz = (tetra_exz[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_xz = trunc(Int, xz)
                xz = (tetra_exz[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_xz = trunc(Int, xz)

                yy = (tetra_eyy[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_yy = trunc(Int, yy)
                yy = (tetra_eyy[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_yy = trunc(Int, yy)

                yz = (tetra_eyz[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_yz = trunc(Int, yz)
                yz = (tetra_eyz[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_yz = trunc(Int, yz)

                zz = (tetra_ezz[1]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)-1
                iemin_zz = trunc(Int, zz)
                zz = (tetra_ezz[4]-TDF_Emin)/(TDF_Emax-TDF_Emin)*(TDF_EneNum-1)+1
                iemax_zz = trunc(Int, zz)

                # xx elements
                if iemin_xx < 0
                    iemin_xx = 0
                end
                if iemax_xx >= TDF_EneNum
                    iemax_xx = TDF_EneNum - 1
                end
                if (0 <= iemin_xx < TDF_EneNum) && (0 <= iemax_xx < TDF_EneNum)
                    for ie = iemin_xx:iemax_xx
                        resultxx = ATM_Spectrum(tetra_exx, tetra_axx, TDF_Energy[ie+1])
                        TDF[ie+1,spin,1] += resultxx
                    end
                end

                # xy elements
                if iemin_xy < 0
                    iemin_xy = 0
                end
                if iemax_xy >= TDF_EneNum
                    iemax_xy = TDF_EneNum - 1
                end
                if (0 <= iemin_xy < TDF_EneNum) && (0 <= iemax_xy < TDF_EneNum)
                    for ie = iemin_xy:iemax_xy
                        resultxy = ATM_Spectrum(tetra_exy, tetra_axy, TDF_Energy[ie+1])
                        TDF[ie+1,spin,2] += resultxy
                    end
                end

                # xz elements
                if iemin_xz < 0
                    iemin_xz = 0
                end
                if iemax_xz >= TDF_EneNum
                    iemax_xz = TDF_EneNum - 1
                end
                if (0 <= iemin_xz < TDF_EneNum) && (0 <= iemax_xz < TDF_EneNum)
                    for ie = iemin_xz:iemax_xz
                        resultxz = ATM_Spectrum(tetra_exz, tetra_axz, TDF_Energy[ie+1])
                        TDF[ie+1,spin,4] += resultxz
                    end
                end

                # yy elements
                if iemin_yy < 0
                    iemin_yy = 0
                end
                if iemax_yy >= TDF_EneNum
                    iemax_yy = TDF_EneNum - 1
                end
                if (0 <= iemin_yy < TDF_EneNum) && (0 <= iemax_yy < TDF_EneNum)
                    for ie = iemin_yy:iemax_yy
                        resultyy = ATM_Spectrum(tetra_eyy, tetra_ayy, TDF_Energy[ie+1])
                        TDF[ie+1,spin,3] += resultyy
                    end
                end

                # yz elements
                if iemin_yz < 0
                    iemin_yz = 0
                end
                if iemax_yz >= TDF_EneNum
                    iemax_yz = TDF_EneNum - 1
                end
                if (0 <= iemin_yz < TDF_EneNum) && (0 <= iemax_yz < TDF_EneNum)
                    for ie = iemin_yz:iemax_yz
                        resultyz = ATM_Spectrum(tetra_eyz, tetra_ayz, TDF_Energy[ie+1])
                        TDF[ie+1,spin,5] += resultyz
                    end
                end

                # zz elements
                if iemin_zz < 0
                    iemin_zz = 0
                end
                if iemax_zz >= TDF_EneNum
                    iemax_zz = TDF_EneNum - 1
                end
                if (0 <= iemin_zz < TDF_EneNum) && (0 <= iemax_zz < TDF_EneNum)
                    for ie = iemin_zz:iemax_zz
                        resultzz = ATM_Spectrum(tetra_ezz, tetra_azz, TDF_Energy[ie+1])
                        TDF[ie+1,spin,6] += resultzz
                    end
                end
            end
        end
    end

    

    factor = 1/knum_i/knum_j/knum_k/6
    for spin = 1:spinsize, ie = 1:TDF_EneNum
        TDF[ie,spin,1] = TDF[ie,spin,1]*factor*tau/cell_volume_Ang   # xx elements
        TDF[ie,spin,2] = TDF[ie,spin,2]*factor*tau/cell_volume_Ang   # xy elements
        TDF[ie,spin,3] = TDF[ie,spin,3]*factor*tau/cell_volume_Ang   # yy elements
        TDF[ie,spin,4] = TDF[ie,spin,4]*factor*tau/cell_volume_Ang   # xz elements
        TDF[ie,spin,5] = TDF[ie,spin,5]*factor*tau/cell_volume_Ang   # yz elements
        TDF[ie,spin,6] = TDF[ie,spin,6]*factor*tau/cell_volume_Ang   # zz elements
    end



    return TDF_Energy, TDF
end