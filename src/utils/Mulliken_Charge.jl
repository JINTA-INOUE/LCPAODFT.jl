mutable struct Mulliken_Charge
    Natom::Int32
    Nspin::Int32
    SpinPol::String
    # Atom_symbol::Vector{String}
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    Total_NumOrbs::Vector{Int32}
    Atoms_Core_Charge::Vector{Float64}
    MulP::Vector{Float64}
    DecMulP::Vector{Vector{Vector{Float64}}}
    InitN_USpin::Vector{Float64}
    InitN_DSpin::Vector{Float64}
    Angle_Spin::Vector{Vector{Float64}}
    Total_SpinS::Float64
    Total_SpinAngle0::Float64
    Total_SpinAngle1::Float64
    TotalZ::Int32
end


function Mulliken_Charge(SpinPol::String, system_grid::System_Grid, Atoms_Core_Charge)

    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs

    TotalZ = sum(Atoms_Core_Charge)

    if SpinPol == "off"
        Nspin = 1
    elseif SpinPol == "on"
        Nspin = 2
    elseif SpinPol == "nc"
        Nspin = 4
    end

    MulP = zeros(Float64, 4)
    DecMulP = Vector{Vector{Vector{Float64}}}(undef, Nspin)
    for spin = 1:Nspin
        DecMulP[spin] = Vector{Vector{Float64}}(undef, Natom)
        for atom = 1:Natom
            DecMulP[spin][atom] = zeros(Float64, maximum(Total_NumOrbs))
        end
    end

    InitN_USpin = zeros(Float64, Natom)
    InitN_DSpin = zeros(Float64, Natom)
    
    Angle_Spin = Vector{Vector{Float64}}(undef, Natom)
    for atom = 1:Natom
        Angle_Spin[atom] = zeros(Float64, 2)
    end


    return Mulliken_Charge(
        Natom, Nspin, SpinPol, FNAN, natn, Total_NumOrbs, Atoms_Core_Charge,
        MulP, DecMulP, InitN_USpin, InitN_DSpin, Angle_Spin, 0.0, 0.0, 0.0, TotalZ)
end



function Mulliken_Charge!(mulliken_charge::Mulliken_Charge, DM, OLP)

    Natom = mulliken_charge.Natom
    Nspin = mulliken_charge.Nspin
    SpinPol = mulliken_charge.SpinPol
    FNAN = mulliken_charge.FNAN
    natn = mulliken_charge.natn
    Total_NumOrbs = mulliken_charge.Total_NumOrbs
    
    DecMulP = mulliken_charge.DecMulP
    InitN_USpin = mulliken_charge.InitN_USpin
    InitN_DSpin = mulliken_charge.InitN_DSpin
    Angle_Spin = mulliken_charge.Angle_Spin

    MulP = zeros(Float64, Nspin, Natom)

    for spin = 1:Nspin, atom = 1:Natom
        fill!(DecMulP[spin][atom], 0.0)
    end


    Total_SpinSx = 0.0
    Total_SpinSy = 0.0
    Total_SpinSz = 0.0
    Total_SpinS = 0.0


    for spin = 1:Nspin
        hst = 0
        for atom = 1:Natom, Rn = 1:FNAN[atom]+1
            jatom = natn[atom][Rn]
            for ist = 1:Total_NumOrbs[atom]
                tmp0 = 0.0
                for jst = 1:Total_NumOrbs[jatom]
                    hst += 1
                    tmp0 += DM[spin][hst] * OLP[hst]
                end


                if spin == 4
                    DecMulP[spin][atom][ist] -= tmp0
                    MulP[spin,atom] -= tmp0
                else
                    DecMulP[spin][atom][ist] += tmp0
                    MulP[spin,atom] += tmp0
                end
            end
        end
    end



    if SpinPol == "off"

        for atom = 1:Natom
            InitN_USpin[atom] = MulP[1,atom]
            InitN_DSpin[atom] = MulP[1,atom]

            Total_SpinS += 0.5*(InitN_USpin[atom] - InitN_DSpin[atom])
        end

        mulliken_charge.Total_SpinS = Total_SpinS

    elseif SpinPol == "on"

        for atom = 1:Natom
            InitN_USpin[atom] = MulP[1,atom]
            InitN_DSpin[atom] = MulP[2,atom]

            Total_SpinS += 0.5*(InitN_USpin[atom] - InitN_DSpin[atom])
        end

        mulliken_charge.Total_SpinS = Total_SpinS

    elseif SpinPol == "nc"

        for atom = 1:Natom

            Nup, Ndown, theta, phi = EulerAngle_Spin(MulP[1,atom], MulP[2,atom], MulP[3,atom], MulP[4,atom])

            MulP[1,atom] = Nup
            MulP[2,atom] = Ndown
            MulP[3,atom] = theta
            MulP[4,atom] = phi
                
            for ist = 1:Total_NumOrbs[atom]

                Nup, Ndown, theta, phi = EulerAngle_Spin(DecMulP[1][atom][ist], DecMulP[2][atom][ist], DecMulP[3][atom][ist], DecMulP[4][atom][ist])

                DecMulP[1][atom][ist] = Nup
                DecMulP[2][atom][ist] = Ndown
                DecMulP[3][atom][ist] = theta
                DecMulP[4][atom][ist] = phi
            end

            InitN_USpin[atom] = MulP[1,atom]
            InitN_DSpin[atom] = MulP[2,atom]
            Angle_Spin[atom][1] = MulP[3,atom]
            Angle_Spin[atom][2] = MulP[4,atom]

            theta = Angle_Spin[atom][1]
            phi = Angle_Spin[atom][2]
            sden = 0.5*(InitN_USpin[atom] - InitN_DSpin[atom])
            Total_SpinSx += sden*sin(theta)*cos(phi)
            Total_SpinSy += sden*sin(theta)*sin(phi)
            Total_SpinSz += sden*cos(theta)
        end

        Total_SpinS, Total_SpinAngle0, Total_SpinAngle1 = xyz_to_spherical(Total_SpinSx, Total_SpinSy, Total_SpinSz)
        mulliken_charge.Total_SpinS = Total_SpinS
        mulliken_charge.Total_SpinAngle0 = Total_SpinAngle0
        mulliken_charge.Total_SpinAngle1 = Total_SpinAngle1
    else
        error("please check SpinPol")
    end
end