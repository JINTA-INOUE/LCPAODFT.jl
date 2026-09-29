function Calc_ChemP(electron::ClusterBloch, kpoints::KPoints)
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    ChemP = Calc_ChemP(electron.spinsize, electron.Spindeg, electron.Nfsize, electron.TotalZ, Beta, electron.Enk)
    return ChemP
end


function Calc_ChemP(spinsize, Spindeg, Nfsize, TZ, Beta, Enk; loopmax=2000)
    
    Beta_trial = 1.0/kb/1200.0*eV2Hartree
    ChemP = 0.0
    ChemP_max = 30.0
    ChemP_min = -30.0
    loopN = 0

    while loopN < loopmax
        
        loopN += 1
        ChemP = 0.5*(ChemP_min + ChemP_max)
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, μ = 1:Nfsize
                
            x = (Enk[μ,spin] - ChemP)*Beta_trial
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF
        end

        Num_state = Spindeg*Num_state
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-14
            break
        end
    end


    ChemP_max = 30.0
    ChemP_min = -30.0
    loopN = 0
    
    while loopN < loopmax
        
        loopN += 1
        if loopN ≠ 1
            ChemP = 0.5*(ChemP_min + ChemP_max)
        end
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, μ = 1:Nfsize
                
            x = (Enk[μ,spin] - ChemP)*Beta
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF
        end

        Num_state = Spindeg*Num_state
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-14
            break
        end
    end

    
    return ChemP
end


function Calc_ChemP(electron::CrystalBloch, kpoints::KPoints)
    Beta = 1/electron.E_Temp/kb*eV2Hartree
    ChemP = Calc_ChemP(electron.spinsize, electron.Spindeg, electron.Nfsize, kpoints.Nkpt, electron.TotalZ, kpoints.All_kweight, Beta, electron.Enk)
    return ChemP
end


function Calc_ChemP(spinsize, Spindeg, Nfsize, Nkpt, TZ, kweight, Beta, Enk; loopmax=2000)
    
    All_Nkpt = sum(kweight)
    Beta_trial = 1.0/kb/3000.0*eV2Hartree

    ChemP = 0.0
    ChemP_max = 20.0
    ChemP_min = -20.0
    loopN = 0

    while loopN < loopmax
        
        loopN += 1
        ChemP = 0.5*(ChemP_min + ChemP_max)
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:Nfsize
                
            x = (Enk[μ,ik,spin] - ChemP)*Beta_trial
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF*kweight[ik]
        end

        Num_state = Spindeg*Num_state/All_Nkpt
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-12
            break
        end
    end


    ChemP_max = 20.0
    ChemP_min = -20.0
    loopN = 0
    
    while loopN < loopmax
        
        loopN += 1
        if loopN ≠ 1
            ChemP = 0.5*(ChemP_min + ChemP_max)
        end
        Num_state = 0.0

        @inbounds for spin = 1:spinsize, ik = 1:Nkpt, μ = 1:Nfsize
                
            x = (Enk[μ,ik,spin] - ChemP)*Beta
            if x <= -max_x
                x = -max_x
            end

            if x >= max_x
                x = max_x
            end

            FermiF = 1/(1 + exp(x))
            Num_state += FermiF*kweight[ik]
        end

        Num_state = Spindeg*Num_state/All_Nkpt
        Dnum = TZ - Num_state

        if Dnum >= 0.0
            ChemP_min = ChemP
        else
            ChemP_max = ChemP
        end

        if abs(Dnum) < 1e-12
            break
        end
    end

    
    return ChemP
end