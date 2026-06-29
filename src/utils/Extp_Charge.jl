function symMatrix_inv!(N, IA)
    val, vec = eigen(Symmetric(IA))
    for i = 1:N
        if abs(val[i]) < 1.0e-15
            val[i] = 0.0
        else
            val[i] = 1/val[i]
        end
    end

    for i = 1:N, j = 1:N
        Sum = 0.0
        for k = 1:N
            Sum += vec[i,k]*val[k]*vec[j,k]
        end
        IA[i,j] = Sum
    end
end


function Extp_Charge(Opt_iter, Extra_CHistory, Natom, His_Gxyz, Gxyz)

    extpln_coes = zeros(Float64, Extra_CHistory+2)
    extpln_coes[1] = 1.0

    if 1 < Extra_CHistory < Opt_iter
        A = zeros(Float64, Extra_CHistory+2, Extra_CHistory+2)
        B = zeros(Float64, Extra_CHistory+2)

        NumHis = Opt_iter - 1
        if Extra_CHistory < NumHis
            NumHis = Extra_CHistory
        end

        for i = 1:NumHis, j = i:NumHis
            Sum = 1.0e-16*(rand()-0.5)
            for k = 1:3*Natom
                Sum += His_Gxyz[i][k]*His_Gxyz[j][k]
            end

            A[i,j] = Sum
            A[j,i] = Sum
        end


        for i = 1:NumHis
            Sum = 0.0
            k = 0
            for j = 1:Natom
                Sum += His_Gxyz[i][k+1]*Gxyz[j][1]
                Sum += His_Gxyz[i][k+2]*Gxyz[j][2]
                Sum += His_Gxyz[i][k+3]*Gxyz[j][3]
                k += 3
            end

            B[i] = Sum
        end

        IA = zeros(Float64, NumHis, NumHis)
        for i = 1:NumHis, j = 1:NumHis
            IA[i,j] = A[i,j]
        end

        symMatrix_inv!(NumHis, IA)


        for i = 1:NumHis
            Sum = 0.0
            for j = 1:NumHis
                Sum += IA[i,j]*B[j]
            end
            
            extpln_coes[i] = Sum
        end


        flag_nan = false
        for i = 1:NumHis
            flag_nan = flag_nan || isnan(extpln_coes[i]) || isinf(extpln_coes[i])
        end

        if flag_nan
            for i = 1:Extra_CHistory
                extpln_coes[i] = 0.0
            end
            extpln_coes[begin] = 1.0
        end

        flag_improper = 0
        for i = 1:NumHis
            if abs(extpln_coes[i]) > 10.0
                flag_improper = 1
            end
        end

        Sum = 0.0
        for i = 1:NumHis
            Sum += extpln_coes[i]
        end

        flag_improper = false
        if abs(Sum-1.0) > 0.001
            flag_improper = true
        end

        if flag_improper
            for i = 1:Extra_CHistory
                extpln_coes[i] = 0.0
            end
            extpln_coes[begin] = 1.0
        end
    end


    for i = Extra_CHistory-1:-1:1
        for k = 1:3*Natom
            His_Gxyz[i+1][k] = His_Gxyz[i][k]
        end
    end

    k = 0
    for i = 1:Natom
        His_Gxyz[1][k+1] = Gxyz[i][1]
        His_Gxyz[1][k+2] = Gxyz[i][2]
        His_Gxyz[1][k+3] = Gxyz[i][3]
        k += 3
    end


    return extpln_coes
end
