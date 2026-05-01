@timeit timer "Poisson" function Solve_Poisson!(dft_mixing::Rhor_Mixing, dDensity_Grid, dVHart_Grid)
    Ngrid = dft_mixing.Ngrid
    Recvecs = dft_mixing.Recvecs
    Density_xyz = dft_mixing.Density_xyz
    Solve_Poisson!(Ngrid, Recvecs, Density_xyz, dDensity_Grid, dVHart_Grid)
end


function Solve_Poisson!(Ngrid, Recvecs, Density_xyz, dDensity_Grid, dVHart_Grid)
    
    Ngrid1, Ngrid2, Ngrid3 = Ngrid
    NN = prod(Ngrid)


    counts = 1
    for i = 1:Ngrid1, j = 1:Ngrid2, k = 1:Ngrid3
        Density_xyz[i,j,k] = dDensity_Grid[counts]
        counts += 1
    end


    in1 = zeros(ComplexF64, Ngrid1)
    in2 = zeros(ComplexF64, Ngrid2)
    in3 = zeros(ComplexF64, Ngrid3)

    
    planfw = plan_fft!( in3, flags=FFTW.MEASURE )
    for j = 1:Ngrid2, i = 1:Ngrid1
        @. @views in3 = Density_xyz[i,j,:]
        R_to_G!(planfw, in3)
        @. @views Density_xyz[i,j,:] = in3
    end

    planfw = plan_fft!( in1, flags=FFTW.MEASURE )
    for k = 1:Ngrid3, j = 1:Ngrid2
        @. @views in1 = Density_xyz[:,j,k]
        R_to_G!(planfw, in1)
        @. @views Density_xyz[:,j,k] = in1
    end

    planfw = plan_fft!( in2, flags=FFTW.MEASURE )
    for i = 1:Ngrid1, k = 1:Ngrid3
        @. @views in2 = Density_xyz[i,:,k]
        R_to_G!(planfw, in2)
        @. @views Density_xyz[i,:,k] = in2
    end

    
    Ngrid11 = div(Ngrid1, 2)
    Ngrid22 = div(Ngrid2, 2)
    Ngrid33 = div(Ngrid3, 2)
    for k = 0:Ngrid3-1, j = 0:Ngrid2-1, i = 0:Ngrid1-1

        sk1 = ifelse(i < Ngrid11, i, i-Ngrid1)
        sk2 = ifelse(j < Ngrid22, j, j-Ngrid2)
        sk3 = ifelse(k < Ngrid33, k, k-Ngrid3)
                
        Gx = sk1*Recvecs[1,1] + sk2*Recvecs[2,1] + sk3*Recvecs[3,1]
        Gy = sk1*Recvecs[1,2] + sk2*Recvecs[2,2] + sk3*Recvecs[3,2] 
        Gz = sk1*Recvecs[1,3] + sk2*Recvecs[2,3] + sk3*Recvecs[3,3]

        Density_xyz[i+1,j+1,k+1] = 4*pi/NN*Density_xyz[i+1,j+1,k+1]/(Gx^2 + Gy^2 + Gz^2)
    end
    Density_xyz[1,1,1] = 0.0



    planbw = plan_ifft!( in3, flags=FFTW.MEASURE )
    for j = 1:Ngrid2, i = 1:Ngrid1
        @. @views in3 = Density_xyz[i,j,:]
        G_to_R!(planbw, in3)
        @. @views Density_xyz[i,j,:] = in3
    end

    planbw = plan_ifft!( in2, flags=FFTW.MEASURE )
    for i = 1:Ngrid1, k = 1:Ngrid3
        @. @views in2 = Density_xyz[i,:,k]
        G_to_R!(planbw, in2)
        @. @views Density_xyz[i,:,k] = in2
    end

    planbw = plan_ifft!( in1, flags=FFTW.MEASURE )
    for k = 1:Ngrid3, j = 1:Ngrid2
        @. @views in1 = Density_xyz[:,j,k]
        G_to_R!(planbw, in1)
        @. @views Density_xyz[:,j,k] = in1
    end



    counts = 1
    for i = 1:Ngrid1, j = 1:Ngrid2, k = 1:Ngrid3
        dVHart_Grid[counts] = real(Density_xyz[i,j,k]) * NN
        counts += 1
    end
end