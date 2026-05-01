function G_to_R!(planbw, fG::Vector{ComplexF64})
    planbw*fG
    return
end


function R_to_G!(planfw, fR::Vector{ComplexF64})
    planfw*fR
    return
end