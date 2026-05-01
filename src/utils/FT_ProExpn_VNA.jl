function FT_ProExpn_VNA!(pao::PAO, pspot::Pspot, maxL, Projector_VNA, Spe_VNA_Bessel)

    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Num_Mesh_VPS = Int64(pspot.Spe_Num_Mesh_VPS)
    Spe_VPS_RV = pspot.Spe_VPS_RV

    GL_Abscissae = zeros(Float64, GL_Mesh)
    GL_Weight = zeros(Float64, GL_Mesh)
    GL_Abscissae, GL_Weight = Gauss_Legendre(GL_Mesh)

    kmin = Radial_kmin
    kmax = Nkmax
    Sk = kmax + kmin
    Dk = kmax - kmin

    rmin = Spe_VPS_RV[begin]
    rmax = Spe_Atom_Cut1 + 0.5
    Sr = rmax + rmin
    Dr = rmax - rmin

    SphB = zeros(Float64, GL_Mesh)
    RGL = zeros(Float64, GL_Mesh)
    for i = 1:GL_Mesh
        RGL[i] = 0.5*(Dr*GL_Abscissae[i] + Sr)
    end


    TmpPhiF = Vector{Vector{Vector{Float64}}}(undef, maxL+1)
    for l = 0:maxL
        TmpPhiF[l+1] = Vector{Vector{Float64}}(undef, maxM)
        for m = 1:maxM
            TmpPhiF[l+1][m] = zeros(Float64, GL_Mesh)
        end
    end

    for l = 0:maxL, m = 1:maxM
        for i = 1:GL_Mesh
            r = RGL[i]
            TmpPhiF[l+1][m][i] = PhiF(Spe_Num_Mesh_VPS, r, Spe_VPS_RV, Projector_VNA[l+1][m])
        end
    end



    for ik = 1:GL_Mesh
        k = 0.5*(Dk*GL_Abscissae[ik] + Sk)
        for l = 0:maxL
            @. SphB = SphericalBesselj(l, k*RGL)
            @. SphB = RGL*RGL*SphB*GL_Weight
            for m = 1:maxM
                Sum = dot(TmpPhiF[l+1][m], SphB)
                Spe_VNA_Bessel[l+1][m][ik] = 0.5*Dr*Sum
            end
        end
    end
end


function Calc_ProExpn_VNA!(pao::PAO, pspot::Pspot, maxL, Projector_VNA, VNA_proj_ene)
    
    Spe_PAO_LMAX = pao.Spe_PAO_Lmax
    Spe_PAO_Mul = pao.Spe_PAO_Mul
    Spe_Num_Mesh_PAO = pao.Spe_Num_Mesh_PAO
    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_PAO_RWF = pao.Spe_PAO_RWF
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV


    GL_Abscissae = zeros(Float64, GL_Mesh)
    GL_Weight = zeros(Float64, GL_Mesh)
    GL_Abscissae, GL_Weight = Gauss_Legendre(GL_Mesh)

    Spe_Vna = Calc_Spe_Vna(pao, pspot)

    phi = Vector{Vector{Float64}}(undef, maxM)
    for m = 1:maxM
        phi[m] = zeros(Float64, Spe_Num_Mesh_VPS)
    end
    phi2 = zeros(Float64, Spe_Num_Mesh_VPS)
    pe = zeros(Float64, maxM)
    


    for L = 0:maxL
        for m = 1:maxM
            mm = m - 1
            if L <= Spe_PAO_LMAX && mm < Spe_PAO_Mul
                for i = 1:Spe_Num_Mesh_VPS
                    r = Spe_VPS_RV[i]
                    phi[m][i] = RadialF(L, Spe_Num_Mesh_PAO, r, Spe_PAO_RV, Spe_PAO_RWF[L+1][m])
                end
            elseif L <= Spe_PAO_LMAX
                for i = 1:Spe_Num_Mesh_VPS
                    phi[m][i] = (0.1*Spe_Vna[i]+1.0e-13)^mm * phi[1][i]
                end
            elseif mm < Spe_PAO_Mul
                for i = 1:Spe_Num_Mesh_VPS
                    r = Spe_VPS_RV[i]
                    temp = RadialF(Spe_PAO_LMAX, Spe_Num_Mesh_PAO, r, Spe_PAO_RV, Spe_PAO_RWF[Spe_PAO_LMAX+1][m])
                    phi[m][i] = temp * r^(L-Spe_PAO_LMAX)
                end
            else
                for i = 1:Spe_Num_Mesh_VPS
                    r = Spe_VPS_RV[i]
                    temp = RadialF(Spe_PAO_LMAX, Spe_Num_Mesh_PAO, r, Spe_PAO_RV, Spe_PAO_RWF[Spe_PAO_LMAX+1][Spe_PAO_Mul])
                    phi[m][i] = temp * (0.1*Spe_Vna[i]+1.0e-13)^(mm-Spe_PAO_Mul+1)
                end
            end
        end

        
        for m = 1:maxM
            dum0 = Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, phi[m], phi[m])
            dum0 = ifelse(dum0>1.0e-17, 1/sqrt(dum0), 0.0)
            for i = 1:Spe_Num_Mesh_VPS
                phi[m][i] = dum0*phi[m][i]
            end
        end

        for i = 1:Spe_Num_Mesh_VPS
            Projector_VNA[L+1][1][i] = phi[1][i]
            phi2[i] = phi[1][i]*Spe_Vna[i]
        end

        dum0 = Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, Projector_VNA[L+1][1], phi2)
        pe[1] = ifelse(abs(dum0)<1.0e-15, 0.0, 1/dum0)


        for m = 2:maxM
            fill!(Projector_VNA[L+1][m], 0.0)
            for n = 1:m-1
                for i = 1:Spe_Num_Mesh_VPS
                    phi2[i] = phi[m][i]*Spe_Vna[i]
                end
                dum0 = Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, Projector_VNA[L+1][n], phi2)
                for i = 1:Spe_Num_Mesh_VPS
                    Projector_VNA[L+1][m][i] += Projector_VNA[L+1][n][i]*pe[n]*dum0
                end
            end

            for i = 1:Spe_Num_Mesh_VPS
                Projector_VNA[L+1][m][i] = phi[m][i] - Projector_VNA[L+1][m][i]
                phi2[i] = Projector_VNA[L+1][m][i]*Spe_Vna[i]
            end

            dum0 = Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, Projector_VNA[L+1][m], phi2)
            pe[m] = ifelse(abs(dum0)<1.0e-15, 0.0, 1/dum0)
        end


        for m = 1:maxM
            dum0 = Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, Projector_VNA[L+1][m], Projector_VNA[L+1][m])
            for i = 1:Spe_Num_Mesh_VPS
                Projector_VNA[L+1][m][i] = dum0*Projector_VNA[L+1][m][i]
            end

            VNA_proj_ene[L+1][m] = ifelse(abs(dum0)<1.0e-15, 0.0, pe[m]/(dum0*dum0))
        end


        for m = 1:maxM
            for i = 1:Spe_Num_Mesh_VPS
                Projector_VNA[L+1][m][i] = Spe_Vna[i]*Projector_VNA[L+1][m][i]
            end
        end
    end
end