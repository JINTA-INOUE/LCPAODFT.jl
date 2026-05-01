function Generate_CWF(cwf_setup::CWF_Setup)


    material = cwf_setup.material
    AtomOrbital = material.AtomOrbital
    xc_type = material.xc_type
    
    GNatom = cwf_setup.GNatom
    GNspecies = cwf_setup.GNspecies
    Guide_Orbs_scale = cwf_setup.Guide_Orbs_scale
    gtype = cwf_setup.gtype
    CWF_HmnR = cwf_setup.CWF_HmnR
    CWF_Wannier = cwf_setup.CWF_Wannier
    filename = cwf_setup.filename
    verbose = cwf_setup.verbose
    
        
    @show AtomOrbital
    @show xc_type

    
    pao, _ = PAO_Pspot(AtomOrbital, xc_type, false)


    if gtype == "AO"

        Spe_Symbol, Spe_cutoff, Spe_orb, Spe_extra, Spe_scale = Get_Guide_Atoms_data(unique(Guide_Orbs_scale))
        proj_pao = Vector{proj_PAO}(undef, GNspecies)
        for spe = 1:GNspecies
            @show Spe_Symbol[spe], Spe_cutoff[spe], join(Spe_orb[spe]), Spe_extra[spe], Spe_scale[spe]
            proj_pao[spe] = Read_ProjPAO(Spe_Symbol[spe], Spe_cutoff[spe], join(Spe_orb[spe]), Spe_extra[spe], Spe_scale[spe])
        end

        projOLP = Set_projOLP(pao, proj_pao, cwf_setup)

    elseif gtype == "HO"
        projOLP = Set_projOLP(cwf_setup)
    else
        error("please check gtype")
    end
    



    Enk, Cnk = Calc_Enk_Cnk(cwf_setup)
    
    
    Amnk = Calc_Amnk(Enk, Cnk, projOLP, cwf_setup)
    # Rotate_Amnk!(Amnk, cwf_setup)
    Σmk, Umnk = Calc_Σmk_Umnk(Amnk, cwf_setup)
    


    # Calculate DM Function
    DM = Calc_DM(Σmk, cwf_setup)
    

    #=
    if CWF_HmnR
        HmnR = Calc_HmnR(Enk, Umnk, cwf_setup)
        # Write_CWF_HmnR(HmnR, DM, cwf_setup)
        Write_CWF_binary(HmnR, DM, cwf_setup)
    end



    
    if CWF_Wannier

        Total_NumOrbs = material.Total_NumOrbs
        
        ucell = Set_UCell(cwf_setup)
        Orbs_Grid = Set_Orbitals_Grid(Total_NumOrbs, pao, ucell)

        CWF_ExpnCoef = Set_CWF_ExpnCoef(Cnk, Umnk, cwf_setup)
            
        Set_CWF_Grid(CWF_ExpnCoef, Orbs_Grid, ucell, cwf_setup)
    end
    =#
end