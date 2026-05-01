function FT_PAO!(Spe_Atom_Cut1, Spe_MaxL_Basis, Spe_Num_Basis, Spe_PAO_RV, Spe_PAO_RWF, Spe_RF_Bessel)

    Spe_Num_Mesh_PAO = length(Spe_PAO_RV)

    rmin = Spe_PAO_RV[begin]
    rmax = Spe_Atom_Cut1 + 0.5
    dk = Nkmax/NkGrid
    h = (rmax - rmin)/OneD_Grid
    r1 = zeros(Float64, OneD_Grid+1)
    r2 = zeros(Float64, OneD_Grid+1)
    for ir = 1:OneD_Grid+1
        r1[ir] = (ir-1)*h + rmin
    end
    @. r2 = r1*r1
    r2[begin] = 0.5*r2[begin]
    r2[end] = 0.5*r2[end]

    TmpRF = Vector{Vector{Vector{Float64}}}(undef, Spe_MaxL_Basis+1)
    for l = 0:Spe_MaxL_Basis
        TmpRF[l+1] = Vector{Vector{Float64}}(undef, Spe_Num_Basis[l+1])
        for p = 1:Spe_Num_Basis[l+1]
            TmpRF[l+1][p] = zeros(Float64, OneD_Grid+1)
        end
    end

    for l = 0:Spe_MaxL_Basis, p = 1:Spe_Num_Basis[l+1]
        for i = 1:OneD_Grid+1
            r = r1[i]
            TmpRF[l+1][p][i] = RadialF(l, Spe_Num_Mesh_PAO, r, Spe_PAO_RV, Spe_PAO_RWF[l+1][p])
        end
    end



    SphB = zeros(Float64, OneD_Grid+1)
    SphB2 = zeros(Float64, OneD_Grid+1)

    for ik = 1:NkGrid
        k = (ik-1)*dk
        for l = 0:Spe_MaxL_Basis

            @. SphB = SphericalBesselj(l, k*r1)
            @. SphB2 = SphB*r2

            for p = 1:Spe_Num_Basis[l+1]
                Sum = dot(TmpRF[l+1][p], SphB2)
                Spe_RF_Bessel[l+1][p][ik] = Sum*h
            end
        end
    end
end