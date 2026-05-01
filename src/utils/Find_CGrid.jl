function Find_CGrids!(NOC, Cxyz, CpyCell, Ns, n1, n2, n3, atv, ratv, gtv, Grid_Origin)

    l1, nn1 = Find_ln(Ns[1], n1)
    l2, nn2 = Find_ln(Ns[2], n2)
    l3, nn3 = Find_ln(Ns[3], n3)

    l1 = l1 + 1
    l2 = l2 + 1
    l3 = l3 + 1
    
    if  (CpyCell < abs(l1)) || (CpyCell < abs(l2)) || (CpyCell < abs(l3))
        NOC[1] = 0
        NOC[2] = nn1
        NOC[3] = nn2
        NOC[4] = nn3
    else
        Rn = ratv[CpyCell+l1, CpyCell+l2, CpyCell+l3]
        
        Cxyz[1] = atv[Rn+1][1] + nn1*gtv[1,1] + nn2*gtv[2,1] + nn3*gtv[3,1] + Grid_Origin[1]
        Cxyz[2] = atv[Rn+1][2] + nn1*gtv[1,2] + nn2*gtv[2,2] + nn3*gtv[3,2] + Grid_Origin[2]
        Cxyz[3] = atv[Rn+1][3] + nn1*gtv[1,3] + nn2*gtv[2,3] + nn3*gtv[3,3] + Grid_Origin[3]

        NOC[1] = Rn
        NOC[2] = nn1
        NOC[3] = nn2
        NOC[4] = nn3
    end
end



function Find_ln(N, n0)

    l0 = 0
    n1 = 0

    if n0 < 0
        l0 = -trunc(Int64, (abs(n0)-1)/N) - 1
        n1 = N - (abs(n0) - N*(abs(l0) - 1))
    elseif N <= n0
        l0 = trunc(Int64, n0/N)
        n1 = n0 - N*l0
    else
        l0 = 0
        n1 = n0
    end 
    
   return l0, n1
end

