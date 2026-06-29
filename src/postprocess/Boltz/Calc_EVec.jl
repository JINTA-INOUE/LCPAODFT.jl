function Calc_EVec(boltz_setup::Boltz_Setup, Cnk)

    mat_type = boltz_setup.mat_type
    if mat_type == "LCPAO"
        # return Calc_EVec_LCPAO(boltz_setup, Cnk)
    elseif mat_type == "CWF"
        return Calc_EVec_Wannier(boltz_setup, Cnk)
    elseif mat_type == "HWF"
        # return Calc_EVec_HWFs(boltz_setup)
    else
        error("please check WF_mode")
    end
end


function Calc_EVec_Wannier(boltz_setup::Boltz_Setup, Cnk)

    material = boltz_setup.material
    Nwann = material.Ngsize
    SpinPol = material.SpinPol
    spinsize = ifelse(SpinPol=="on", 2, 1)
    kmesh = boltz_setup.kmesh
    plane_type = boltz_setup.plane_type

    Nkpt = prod(kmesh)
    if plane_type
        if kmesh[3] ≠ 1
            error("please check plane_type, kmesh.")
        end
    end


    EVec = Vector{Vector{Vector{Vector{Float32}}}}(undef, spinsize)
    for spin = 1:spinsize
        EVec[spin] = Vector{Vector{Vector{Float32}}}(undef, Nkpt)
        for ik = 1:Nkpt
            EVec[spin][ik] = Vector{Vector{Float32}}(undef, Nwann)
            for ist = 1:Nwann
                EVec[spin][ik][ist] = zeros(Float32, Nwann)
            end
        end
    end

    @inbounds for spin = 1:spinsize, ik = 1:Nkpt, ist = 1:Nwann, μ = 1:Nwann
        tmp = conj(Cnk[spin][ik][ist,μ]) * Cnk[spin][ik][ist,μ]
        EVec[spin][ik][ist][μ] = real(tmp)
    end
    

    return EVec
end