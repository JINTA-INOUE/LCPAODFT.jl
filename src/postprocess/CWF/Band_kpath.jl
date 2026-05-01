function Calc_Enkpath!(Enk, spin, Nkpath, kpath_Nk, kpath_start, kpath_end, NCell, cell_list_ijk, Nwann, HmnR)

    H = zeros(ComplexF64, Nwann, Nwann)

    for ik = 1:Nkpath, ipath = 1:kpath_Nk[ik]

        fill!(H, 0.0)

        k1 = kpath_start[ik][1] + (kpath_end[ik][1]-kpath_start[ik][1])*(ipath-1)/(kpath_Nk[ik]-1)
        k2 = kpath_start[ik][2] + (kpath_end[ik][2]-kpath_start[ik][2])*(ipath-1)/(kpath_Nk[ik]-1)
        k3 = kpath_start[ik][3] + (kpath_end[ik][3]-kpath_start[ik][3])*(ipath-1)/(kpath_Nk[ik]-1)

        for cell = 1:NCell
            kRn = k1*cell_list_ijk[cell][1] + k2*cell_list_ijk[cell][2] + k3*cell_list_ijk[cell][3]
            ex = cispi(2*kRn)

            for ist = 1:Nwann, jst = 1:Nwann
                H[ist,jst] += HmnR[spin,cell,ist,jst] * ex
            end
        end

        Enk[ik][ipath] = eigvals(Hermitian(H))*Hartree2eV
    end
end


function Band_kpath(filepath::String, kpath::Vector{Vector{Float64}}, kname::Vector{String}; filename = splitext(basename(filepath))[1], PAO_file=Nothing, Nk=30, verbose=true)
    
    spinsize, Recvecs, Erange, ChemP, NCell, cell_list_ijk, NNcell, NNcell_list, Nwann, HmnR, CWF_SOC_method, SO_real, SO_imag = Load_CWF_jld2(filepath)
    println("filename = $filename")


    if verbose
        @show Recvecs
        @show Erange
        @show ChemP
        @show kpath
        @show kname
        @show NCell
        @show NNcell
        @show NNcell_list
        @show Nwann
        @show HmnR[1,1,1,1]
        @show CWF_SOC_method
    end


    @show Nkpath = length(kpath)-1
    @show kpath_Nk = ones(Int64, Nkpath)*Nk


    kpath_start = Vector{Vector{Float64}}(undef, Nkpath)
    kpath_end = Vector{Vector{Float64}}(undef, Nkpath)
    for ik = 1:Nkpath
        kpath_start[ik] = kpath[ik]
        kpath_end[ik] = kpath[ik+1]
    end

    
    Enk = Vector{Vector{Vector{Float64}}}(undef, Nkpath)
    for ik = 1:Nkpath
        Enk[ik] = Vector{Vector{Float64}}(undef, kpath_Nk[ik])
        for ipath = 1:kpath_Nk[ik]
            Enk[ik][ipath] = zeros(Float64, Nwann)
        end
    end

    for spin = 1:spinsize
        Calc_Enkpath!(Enk, spin, Nkpath, kpath_Nk, kpath_start, kpath_end, NCell, cell_list_ijk, Nwann, HmnR)
        write_BANDDAT(filename, spin, kpath_start, kpath_end, kpath_Nk, Nkpath, Nwann, Enk, ChemP, Recvecs)
    end


    write_GNUBAND(filename, spinsize, Erange, Nkpath, kpath, kname, Recvecs; PAO_file)
end