function FFT_Grid!(Ngrid, Gridxyz::AbstractArray{<:ComplexF64,3})

    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    in1 = zeros(ComplexF64, Ngrid1)
    in2 = zeros(ComplexF64, Ngrid2)
    in3 = zeros(ComplexF64, Ngrid3)

    planfw = plan_fft!(in3, flags=FFTW.ESTIMATE)
    @inbounds for j = 1:Ngrid2, i = 1:Ngrid1
        @. @views in3 = Density_xyz[i,j,:]
        R_to_G!(planfw, in3)
        @. @views Gridxyz[i,j,:] = in3
    end

    planfw = plan_fft!(in2, flags=FFTW.ESTIMATE)
    @inbounds for i = 1:Ngrid1, k = 1:Ngrid3
        @. @views in2 = Gridxyz[i,:,k]
        R_to_G!(planfw, in2)
        @. @views Gridxyz[i,:,k] = in2
    end

    planfw = plan_fft!(in1, flags=FFTW.ESTIMATE)
    @inbounds for k = 1:Ngrid3, j = 1:Ngrid2
        @. @views in1 = Gridxyz[:,j,k]
        R_to_G!(planfw, in1)
        @. @views Gridxyz[:,j,k] = in1
    end
end


function iFFT_Grid!(Ngrid, Gridxyz::AbstractArray{<:ComplexF64,3})

    Ngrid1, Ngrid2, Ngrid3 = Ngrid

    in1 = zeros(ComplexF64, Ngrid1)
    in2 = zeros(ComplexF64, Ngrid2)
    in3 = zeros(ComplexF64, Ngrid3)

    planbw = plan_ifft!(in1, flags=FFTW.ESTIMATE)
    @inbounds for k = 1:Ngrid3, j = 1:Ngrid2
        @. @views in1 = Gridxyz[:,j,k]
        G_to_R!(planbw, in1)
        @. @views Gridxyz[:,j,k] = in1
    end

    planbw = plan_ifft!(in2, flags=FFTW.ESTIMATE)
    @inbounds for i = 1:Ngrid1, k = 1:Ngrid3
        @. @views in2 = Gridxyz[i,:,k]
        G_to_R!(planbw, in2)
        @. @views Gridxyz[i,:,k] = in2
    end

    planbw = plan_ifft!(in3, flags=FFTW.ESTIMATE)
    @inbounds for j = 1:Ngrid2, i = 1:Ngrid1
        @. @views in3 = Gridxyz[i,j,:]
        G_to_R!(planbw, in3)
        @. @views Gridxyz[i,j,:] = in3
    end
end