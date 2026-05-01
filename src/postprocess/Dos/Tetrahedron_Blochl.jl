"""
/*****************************************************************
   Improved analytical tetrahedron method 
   Idea: PRB 49, 16223 (1994), by Bloechl, Jepsen and Andersen 
*******************************************************************/

ref: OpenMX/Tetrahedron_Blochl.c
"""


function OrderE0!(e, n)
    for i = 1:n, j = i:n
        if e[j] < e[i]
            t = e[j]
            e[j] = e[i]
            e[i] = t
        end
    end
end


function OrderE!(e, a, n)
    for i = 1:n, j = i:n
        if e[j] < e[i]
            t = e[j]
            e[j] = e[i]
            e[i] = t

            t = a[j]
            a[j] = a[i]
            a[i] = t
        end
    end
end


function ATM_Dos(et, e)

    e21 = et[2] - et[1]
    e31 = et[3] - et[1]
    e32 = et[3] - et[2]
    e41 = et[4] - et[1]
    e42 = et[4] - et[2]
    e43 = et[4] - et[3]
    if e < et[1]
        dos = 0.0
    elseif e > et[1] && e < et[2]
        e1 = e - et[1]
        dos =  3.0 * e1*e1 / (e21 * e31 * e41)
    elseif e > et[2] && e < et[3]
        e2 = e - et[2]
        dos = (e21*3.0 + e2*6.0 - (e31 + e42)*3.0*e2*e2 / (e32*e42)) / (e31*e41)
    elseif e > et[3] && e < et[4]
        e4 = et[4] - e
        dos =  3.0* e4*e4 / (e41 * e42 * e43)
    elseif e > et[4]
        dos = 0.0
    end

    return dos
end


function ATM_Spectrum(et, at, e)

    dos = ATM_Dos(et, e)

    e21 = et[2] - et[1]
    e31 = et[3] - et[1]
    e32 = et[3] - et[2]
    e41 = et[4] - et[1]
    e42 = et[4] - et[2]
    e43 = et[4] - et[3]
    a21 = at[2] - at[1]
    a31 = at[3] - at[1]
    a32 = at[3] - at[2]
    a41 = at[4] - at[1]
    a42 = at[4] - at[2]
    a43 = at[4] - at[3]

    spectrum = 0.0

    if e < et[1]
        spectrum = 0.0
    elseif e > et[1] && e < et[2]
        spectrum = dos * (at[1] + (e - et[1]) * 1/3 * (a21 / e21 + a31 / e31 + a41 / e41))
    elseif e > et[2] && e < et[3]
        spectrum = dos * (at[1] + (a21 + e21 * a31 / e31 + e21 * a41 / e41) *
                 1/3 * (et[3] - e) / e32 + (at[4] - at[1] - (
                e43 * a41 / e41 + e43 * a42 / e42 + a43) * 1/3)
                 * (e - et[2]) / e32)
    elseif e > et[3] && e < et[4]
        spectrum = dos * (at[4] + (e - et[4]) * 1/3 * (a41 / e41 + a42 / e42 + a43 / e43))
    elseif e > et[4]
        spectrum = 0.0
    end

    return spectrum
end