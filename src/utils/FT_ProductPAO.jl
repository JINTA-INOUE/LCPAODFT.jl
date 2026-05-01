function FT_ProductPAO!(pao::PAO, Spe_VPS_RV, Spe_ProductRF_Bessel)

    Spe_MaxL_Basis = pao.Spe_MaxL_Basis
    Spe_Num_Basis = pao.Spe_Num_Basis
    Spe_Num_Mesh_PAO = pao.Spe_Num_Mesh_PAO
    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_PAO_RWF = pao.Spe_PAO_RWF
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    
    kmin = Radial_kmin
    kmax = Nkmax
    Sk = kmax + kmin
    Dk = kmax - kmin

    rmin = Spe_VPS_RV
    rmax = Spe_Atom_Cut1 + 0.5
    Sr = rmax + rmin
    Dr = rmax - rmin

    GL_Abscissae = zeros(Float64, GL_Mesh)
    GL_Weight = zeros(Float64, GL_Mesh)
    GL_Abscissae, GL_Weight = Gauss_Legendre(GL_Mesh)


    RGL = zeros(Float64, GL_Mesh)
    RGL2 = zeros(Float64, GL_Mesh)
    for i = 1:GL_Mesh
        RGL[i] = 0.5*(Dr*GL_Abscissae[i] + Sr)
        RGL2[i] = RGL[i]*RGL[i]*GL_Weight[i]
    end

    
    GL_PAO = Vector{Vector{Vector{Float64}}}(undef, Spe_MaxL_Basis+1)
    for l = 0:Spe_MaxL_Basis
        GL_PAO[l+1] = Vector{Vector{Float64}}(undef, Spe_Num_Basis[l+1])
        for p = 1:Spe_Num_Basis[l+1]
            GL_PAO[l+1][p] = zeros(Float64, GL_Mesh)
        end
    end


    for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
        for i = 1:GL_Mesh
	        GL_PAO[l+1][p][i] = RadialF(l, Spe_Num_Mesh_PAO, RGL[i], Spe_PAO_RV, Spe_PAO_RWF[l+1][p])
        end
    end

    Lmax = 2*Spe_MaxL_Basis+3
    SphB = Vector{Vector{Float64}}(undef, Lmax)
    for l = 1:Lmax
        SphB[l] = zeros(Float64, GL_Mesh)
    end


    for ik = 1:GL_Mesh

        k = 0.50*(Dk*GL_Abscissae[ik] + Sk)
        for l = 0:2*Spe_MaxL_Basis
            @. SphB[l+1] = SphericalBesselj(l, k*RGL)
        end

        for l1 = 0:Spe_MaxL_Basis, p1 = 1:Spe_Num_Basis[l1+1]
            for l2 = 0:Spe_MaxL_Basis
                if l1 <= l2
                    Lmax = 2*l2
                    for p2 = 1:Spe_Num_Basis[l2+1]
                        for l = 0:Lmax
                            Sum = 0.0
                            for i = 1:GL_Mesh
							    sj = SphB[l+1][i]
							    Sum += RGL2[i]*sj*GL_PAO[l1+1][p1][i]*GL_PAO[l2+1][p2][i]
                            end

                            Spe_ProductRF_Bessel[l1+1][p1][l2+1][p2][l+1][ik] = 0.5*Dr*Sum
                        end
                    end
                end
            end
        end
    end
end