function FT_NLP(pao::Vector{PAO}, pspot::Vector{Pspot})
    Nspecies = length(pao)
    NLRF_Bessel = Vector{Vector{Vector{Vector{Float64}}}}(undef, Nspecies)
    for spe = 1:Nspecies
        Spe_Num_RVPS = pspot[spe].Spe_Num_RVPS
        VPS_j_dependency = pspot[spe].VPS_j_dependency
        
        NLRF_Bessel[spe] = Vector{Vector{Vector{Float64}}}(undef, VPS_j_dependency+1)
        for so = 1:VPS_j_dependency+1
            NLRF_Bessel[spe][so] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
            for l = 1:Spe_Num_RVPS
                NLRF_Bessel[spe][so][l] = zeros(Float64, NkGrid+1)
            end
        end
    end

    for spe = 1:Nspecies
        Spe_Atom_Cut1 = pao[spe].Spe_Atom_Cut1
        FT_NLP!(Spe_Atom_Cut1, NLRF_Bessel[spe], pspot[spe])
    end

    return NLRF_Bessel
end


function FT_NLP!(Spe_Atom_Cut1, NLRF_Bessel, pspot::Pspot)

    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Num_RVPS = pspot.Spe_Num_RVPS
    Spe_VPS_List = pspot.Spe_VPS_List
    VPS_j_dependency = pspot.VPS_j_dependency
    Spe_VNL = pspot.Spe_VNL

    dk = Nkmax/NkGrid

    rmin = Spe_VPS_RV[begin]
    rmax = Spe_Atom_Cut1 + 0.5
    h = (rmax - rmin)/OneD_Grid
    r1 = zeros(Float64, OneD_Grid+1)
    r2 = zeros(Float64, OneD_Grid+1)
    for ir = 1:OneD_Grid+1
        r1[ir] = (ir-1)*h + rmin
    end
    @. r2 = r1^2



    TmpNRF = Vector{Vector{Vector{Float64}}}(undef, VPS_j_dependency+1)
    for so = 1:VPS_j_dependency+1
        TmpNRF[so] = Vector{Vector{Float64}}(undef, Spe_Num_RVPS)
        for L = 1:Spe_Num_RVPS
            TmpNRF[so][L] = zeros(Float64, OneD_Grid+1)
        end
    end

    for so = 1:VPS_j_dependency+1, L = 1:Spe_Num_RVPS
        GL = Spe_VPS_List[L]
        for i = 1:OneD_Grid+1
            r = r1[i]
            TmpNRF[so][L][i] = Nonlocal_RadialF(GL, Spe_Num_Mesh_VPS, r, Spe_VPS_RV, Spe_VNL[so][L])
        end
    end


    Lmax = maximum(Spe_VPS_List)
    SphB = zeros(Float64, OneD_Grid+1)
    SphB2 = Vector{Vector{Float64}}(undef, Lmax+1)
    for L = 1:Lmax+1
        SphB2[L] = zeros(Float64, OneD_Grid+1)
    end


    Sum = zeros(Float64, VPS_j_dependency+1)
    for ik = 1:NkGrid
        
        k = (ik-1)*dk
        for GL = 0:Lmax
            @. SphB = SphericalBesselj(GL, k*r1)
            @. SphB2[GL+1] = SphB * r2
        end

        for GL = 1:Lmax+1
            SphB2[GL][begin] = 0.5*SphB2[GL][begin]
            SphB2[GL][end] = 0.5*SphB2[GL][end]
        end
            
        for so = 1:VPS_j_dependency+1, L = 1:Spe_Num_RVPS
            GL = Spe_VPS_List[L]
            Sum = dot(TmpNRF[so][L], SphB2[GL+1])

            NLRF_Bessel[so][L][ik] = Sum*h
        end
    end
end
