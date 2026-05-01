function FT_VNA!(pao::Vector{PAO}, pspot::Vector{Pspot}, Spe_CrudeVNA_Bessel)
    if length(pao) ≠ length(pspot)
        error("please check FT_VNA")
    end

    Nspecies = length(pao)
    for spe = 1:Nspecies
        FT_VNA!(pao[spe], pspot[spe], Spe_CrudeVNA_Bessel[spe])
    end
end


function FT_VNA!(pao::PAO, pspot::Pspot, Spe_CrudeVNA_Bessel)

    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV

    kmin = Radial_kmin
    kmax = Nkmax
    Sk = kmax + kmin
    Dk = kmax - kmin
    xmin = sqrt(Spe_PAO_RV[begin])
    xmax = sqrt(Spe_Atom_Cut1 + 0.5)
    h = (xmax - xmin)/OneD_Grid

    GL_Abscissae = Gauss_Legendre_x(GL_Mesh)
    Spe_Vna = Calc_Spe_Vna(pao, pspot)


    TmpVNAF = zeros(Float64, OneD_Grid+1)
    for i = 1:OneD_Grid+1
        x = xmin + (i-1)*h
        r = x*x
        TmpVNAF[i] = VNAF(Spe_Num_Mesh_VPS, r, Spe_Atom_Cut1, Spe_VPS_RV, Spe_Vna)
    end


    for ik = 1:GL_Mesh
        k = 0.5*(Dk*GL_Abscissae[ik] + Sk)

        r = xmin*xmin
        SphB = SphericalBesselj(0, k*r)
        Sum = TmpVNAF[1]*r*r*xmin*SphB

        r = xmax*xmax
        SphB = SphericalBesselj(0, k*r)
        Sum += TmpVNAF[OneD_Grid+1]*r*r*xmax*SphB

        for i = 2:OneD_Grid
            x = xmin + (i-1)*h
            r = x*x
            SphB = SphericalBesselj(0, k*r)
            Sum += 2*TmpVNAF[i]*r*r*x*SphB
        end

        Spe_CrudeVNA_Bessel[ik] = Sum*h
    end
end