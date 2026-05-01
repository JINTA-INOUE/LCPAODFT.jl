function Spe_VHart_Atom!(Spe_VPS_RV, Spe_PAO_RV, Spe_VH_Atom, Spe_Atom_Den)
    
    FineGL_x, FineGL_Weight = Gauss_Legendre(FineGL_Mesh)

    Num_Mesh_VPS = length(Spe_VPS_RV)
    Num_Mesh_PAO = length(Spe_PAO_RV)

    for i = 1:Num_Mesh_VPS

        # inside
        R = Spe_VPS_RV[i]
        xmin = log(Spe_PAO_RV[begin])
        xmax = log(R)
        Sx = xmax + xmin
        Dx = xmax - xmin
        Inside = 0.0

        for j = 1:FineGL_Mesh
            x = 0.5*(Dx*FineGL_x[j] + Sx)
            rp = exp(x)
            temp = KumoF(Num_Mesh_PAO, x, Spe_PAO_RV, Spe_Atom_Den)
            Inside += temp*FineGL_Weight[j]*rp^3
        end
        Inside = 0.5*Dx*Inside
        Inside = 4*pi*Inside/R


        # outside
        xmin = log(R)
        xmax = log(Spe_PAO_RV[end])
        Sx = xmax + xmin
        Dx = xmax - xmin

        Outside = 0.0
        for j = 1:FineGL_Mesh
            x = 0.5*(Dx*FineGL_x[j] + Sx)
            rp = exp(x)
            temp = KumoF(Num_Mesh_PAO, x, Spe_PAO_RV, Spe_Atom_Den)
            Outside += temp*FineGL_Weight[j]*rp^2
        end
        Outside = 2*Dx*pi*Outside

        Spe_VH_Atom[i+1] = Inside + Outside
    end
end


function Calc_Spe_VH_Atom(pao::PAO, pspot::Pspot)

    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VH_Atom = zeros(Float64, Spe_Num_Mesh_VPS+2)
    Calc_Spe_VH_Atom!(pao, pspot, Spe_VH_Atom)

    return Spe_VH_Atom
end


function Calc_Spe_VH_Atom!(pao::PAO, pspot::Pspot, Spe_VH_Atom)

    Spe_PAO_RV = pao.Spe_PAO_RV
    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Atomic_Den = pao.Spe_Atomic_Den
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Vcore = pspot.Spe_Vcore
    Spe_Core_Charge = pspot.Spe_Core_Charge

    Spe_VHart_Atom!( Spe_VPS_RV, Spe_PAO_RV, Spe_VH_Atom, Spe_Atomic_Den ) 
    Spe_VH_Atom[1] = 2*Spe_VH_Atom[2] - Spe_VH_Atom[3]
    Spe_VH_Atom[end] = 2*Spe_VH_Atom[end-1] - Spe_VH_Atom[end-2]


    # the nearest point to Spe_Atom_Cut1
    dum = 1000.0
    ii = 0
    for i = 1:Spe_Num_Mesh_VPS

        r = Spe_VPS_RV[i]
        dum1 = abs(r - Spe_Atom_Cut1)

        if dum1 < dum
            dum = dum1
            ii = i
        end
    end

    # correct the asymptotic behaviour of Spe_VH_Atom
    if Spe_Core_Charge > 1.0e-15
        dum = -Spe_Vcore[ii]/Spe_VH_Atom[ii+1]
        for i = 1:Spe_Num_Mesh_VPS+2
            Spe_VH_Atom[i] = dum*Spe_VH_Atom[i]
        end
    end
end


function Calc_Spe_Vna(pao::PAO, pspot::Pspot)

    Spe_Atom_Cut1 = pao.Spe_Atom_Cut1
    Spe_Num_Mesh_VPS = pspot.Spe_Num_Mesh_VPS
    Spe_VPS_RV = pspot.Spe_VPS_RV
    Spe_Vcore = pspot.Spe_Vcore

    Spe_VH_Atom = Calc_Spe_VH_Atom(pao, pspot)

    Spe_Vna = zeros(Float64, Spe_Num_Mesh_VPS)
    for i = 1:Spe_Num_Mesh_VPS
        r = Spe_VPS_RV[i]
        temp = 1/(1+exp(20*(r-Spe_Atom_Cut1)))
        Spe_Vna[i] = temp*(Spe_Vcore[i] + Spe_VH_Atom[i+1])
    end


    return Spe_Vna
end