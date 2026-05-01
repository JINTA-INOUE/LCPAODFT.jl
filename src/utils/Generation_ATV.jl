function Generation_ATV!(N::Integer, Latvecs, atv, atv_ijk, ratv)

    Rn = 1
    di = -N-1

    for i = -N:N
        di = di + 1
        dj = -N-1
        for j = -N:N
            dj = dj + 1
            dk = -N-1
            for k = -N:N
                dk = dk + 1

                if i == 0 && j == 0 && k == 0
                    atv[1][1] = 0.0
                    atv[1][2] = 0.0
                    atv[1][3] = 0.0
                    atv_ijk[1][1] = 0
                    atv_ijk[1][2] = 0
                    atv_ijk[1][3] = 0
                    
                    ratv[i+N+1,j+N+1,k+N+1] = 0
                else
                    atv[Rn+1][1] = di*Latvecs[1,1] + dj*Latvecs[2,1] + dk*Latvecs[3,1]
                    atv[Rn+1][2] = di*Latvecs[1,2] + dj*Latvecs[2,2] + dk*Latvecs[3,2]
                    atv[Rn+1][3] = di*Latvecs[1,3] + dj*Latvecs[2,3] + dk*Latvecs[3,3]
                    atv_ijk[Rn+1][1] = i
                    atv_ijk[Rn+1][2] = j
                    atv_ijk[Rn+1][3] = k
                    ratv[i+N+1,j+N+1,k+N+1] = Rn

                    Rn = Rn + 1
                end
            end
        end
    end
end


function Generation_ATV!(N::Integer, Latvecs, atv)

    Rn = 1
    di = -N-1

    for i = -N:N
        di = di + 1
        dj = -N-1
        for j = -N:N
            dj = dj + 1
            dk = -N-1
            for k = -N:N
                dk = dk + 1

                if i == 0 && j == 0 && k == 0
                    atv[1][1] = 0.0
                    atv[1][2] = 0.0
                    atv[1][3] = 0.0
                else
                    atv[Rn+1][1] = di*Latvecs[1,1] + dj*Latvecs[2,1] + dk*Latvecs[3,1]
                    atv[Rn+1][2] = di*Latvecs[1,2] + dj*Latvecs[2,2] + dk*Latvecs[3,2]
                    atv[Rn+1][3] = di*Latvecs[1,3] + dj*Latvecs[2,3] + dk*Latvecs[3,3]

                    Rn = Rn + 1
                end
            end
        end
    end
end


function Generation_ATV_ijk!(N::Integer, atv_ijk)

    Rn = 1
    di = -N-1

    for i = -N:N
        di = di + 1
        dj = -N-1
        for j = -N:N
            dj = dj + 1
            dk = -N-1
            for k = -N:N
                dk = dk + 1

                if i == 0 && j == 0 && k == 0
                    atv_ijk[1][1] = 0
                    atv_ijk[1][2] = 0
                    atv_ijk[1][3] = 0
                else
                    atv_ijk[Rn+1][1] = i
                    atv_ijk[Rn+1][2] = j
                    atv_ijk[Rn+1][3] = k

                    Rn = Rn + 1
                end
            end
        end
    end
end




function Generation_RATV!(N::Integer, ratv)

    Rn = 1
    di = -N-1

    for i = -N:N
        di = di + 1
        dj = -N-1
        for j = -N:N
            dj = dj + 1
            dk = -N-1
            for k = -N:N
                dk = dk + 1

                if i == 0 && j == 0 && k == 0
                    continue
                else
                    ratv[i+N+1,j+N+1,k+N+1] = Rn
                    Rn = Rn + 1
                end
            end
        end
    end
end

