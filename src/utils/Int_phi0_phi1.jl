function Int_phi0_phi1(GL_Abscissae, GL_Weight, Spe_Num_Mesh_VPS, Spe_VPS_RV, phi0, phi1)

    rmin = Spe_VPS_RV[begin]
    rmax = Spe_VPS_RV[end]
    Sr = rmax + rmin
    Dr = rmax - rmin
    Sum = 0.0
    for i = 1:GL_Mesh
        r = 0.5*(Dr*GL_Abscissae[i] + Sr)
        tmp0 = PhiF(Spe_Num_Mesh_VPS, r, Spe_VPS_RV, phi0)
        tmp1 = PhiF(Spe_Num_Mesh_VPS, r, Spe_VPS_RV, phi1)
        Sum += r*r*GL_Weight[i]*tmp0*tmp1
    end

    return 0.5*Dr*Sum
end