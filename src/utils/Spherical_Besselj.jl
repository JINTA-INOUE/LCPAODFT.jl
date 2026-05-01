# For Set_OLP_Kin.jl, Set_Nonlocal.jl, Set_ProExpn_VNA.jl
function SphericalBesselj(l::Integer, x::T)::T where {T}

    iszero(l) && iszero(x) && return one(T)
    iszero(x) && return zero(T)
    
    return sqrt(0.5*pi/x)*besselj(l+0.5,x)
end


@inline function Calc_SphericalBesselj!(lmax::Integer, x, tsb, SphB, dSphB)

    xmin = 0.0
    nmax = lmax + trunc(Int64, 1.5*x) + 10
    if nmax < 30
        nmax = 30
    end

    if x > xmin
        invx = 1/x
        vsbi = 0.0
        vsb0 = 0.0
        vsb1 = 1.0e-14

        for n = nmax-1:-1:lmax+3
            vsb2 = (2.0*n + 1.0)*invx*vsb1 - vsb0
            if 0.0<(vsb2-1.0e+250)
                tmp = 1.0/vsb2
                vsb1 *= tmp
                vsb2 = 1.0
            end
            vsbi = vsb0
            vsb0 = vsb1
            vsb1 = vsb2
        end

        n = lmax + 3
        tsb[n] = vsb1
        tsb[n+1] = vsb0
        tsb[n+2] = vsbi

        tmp = tsb[n]    
        tsb[n] /= tmp
        tsb[n+1] /= tmp

        for n = lmax+2:-1:1
            tsb[n] = (2*n + 1)*invx*tsb[n+1] - tsb[n+2]
            if 1.0e+250<tsb[n]
                tmp = tsb[n]
                for m = n:lmax+1
                    tsb[m] /= tmp
                end
            end
        end

        si = sin(x)
        co = cos(x)
        ix = 1/x
        j0 = si*ix
        j1 = si*ix*ix - co*ix

        sf = ifelse(abs(tsb[2])<abs(tsb[1]), j0/tsb[1], j1/tsb[2])

        for n = 1:lmax+2
            SphB[n] = tsb[n]*sf
        end

        dSphB[1] = co*ix - si*ix*ix
        for n = 1:lmax
            dSphB[n+1] = (n*SphB[n] - (n+1)*SphB[n+2])/(2*n+1)
        end
    else
        for n = 2:lmax+1
            SphB[n] = 0.0
        end
        SphB[begin] = 1.0

        dSphB[1] = 0.0
        for n = 1:lmax
            dSphB[n+1] = (n*SphB[n] - (n+1)*SphB[n+2])/(2*n+1)
        end

        return nothing
    end
end

# For Set_ProExpn_VNA.jl
@inline function Calc_SphericalBesselj2!(lmax::Integer, x, tsb, SphB)

    xmin = 0.0
    nmax = lmax + trunc(Int64, 3*x) + 20
    if nmax < 100
        nmax = 100
    end

    if x > xmin
        invx = 1/x
        vsbi = 0.0
        vsb0 = 0.0
        vsb1 = 1.0e-14

        for n = nmax-1:-1:lmax+3
            vsb2 = (2.0*n + 1.0)*invx*vsb1 - vsb0
            if 0.0<(vsb2-1.0e+250)
                tmp = 1.0/vsb2
                vsb2 *= tmp
                vsb1 *= tmp
            end
            vsbi = vsb0
            vsb0 = vsb1
            vsb1 = vsb2
        end

        n = lmax + 3
        tsb[n] = vsb1
        tsb[n+1] = vsb0
        tsb[n+2] = vsbi

        tmp = tsb[n]    
        tsb[n] /= tmp
        tsb[n+1] /= tmp

        for n = lmax+2:-1:1
            tsb[n] = (2*n + 1)*invx*tsb[n+1] - tsb[n+2]
            if 0.0<(tsb[n]-1.0e+250)
                tmp = tsb[n]
                for m = n-1:lmax+1
                    tsb[m+1] /= tmp
                end
            end
        end

        si = sin(x)
        co = cos(x)
        ix = 1/x
        j0 = si*ix
        j1 = si*ix*ix - co*ix

        sf = ifelse(abs(tsb[2])<abs(tsb[1]), j0/tsb[1], j1/tsb[2])

        for n = 1:lmax+2
            SphB[n] = tsb[n]*sf
        end
    else
        for n = 2:lmax+1
            SphB[n] = 0.0
        end
        SphB[begin] = 1.0

        return nothing
    end
end


# For Set_ProExpn_VNAforce.jl
@inline function Calc_SphericalBesselj2!(lmax::Integer, x, tsb, SphB, dSphB)

    xmin = 0.0
    nmax = lmax + trunc(Int64, 3*x) + 20
    if nmax < 100
        nmax = 100
    end

    if x > xmin
        invx = 1/x
        vsbi = 0.0
        vsb0 = 0.0
        vsb1 = 1.0e-14

        for n = nmax-1:-1:lmax+3
            vsb2 = (2.0*n + 1.0)*invx*vsb1 - vsb0
            if 0.0<(vsb2-1.0e+250)
                tmp = 1.0/vsb2
                vsb2 *= tmp
                vsb1 *= tmp
            end
            vsbi = vsb0
            vsb0 = vsb1
            vsb1 = vsb2
        end

        n = lmax + 3
        tsb[n] = vsb1
        tsb[n+1] = vsb0
        tsb[n+2] = vsbi

        tmp = tsb[n]    
        tsb[n] /= tmp
        tsb[n+1] /= tmp

        for n = lmax+2:-1:1
            tsb[n] = (2*n + 1)*invx*tsb[n+1] - tsb[n+2]
            if 0.0<(tsb[n]-1.0e+250)
                tmp = tsb[n]
                for m = n-1:lmax+1
                    tsb[m+1] /= tmp
                end
            end
        end

        si = sin(x)
        co = cos(x)
        ix = 1/x
        j0 = si*ix
        j1 = si*ix*ix - co*ix

        sf = ifelse(abs(tsb[2])<abs(tsb[1]), j0/tsb[1], j1/tsb[2])

        for n = 1:lmax+2
            SphB[n] = tsb[n]*sf
        end

        dSphB[1] = co*ix - si*ix*ix
        for n = 1:lmax
            dSphB[n+1] = (n*SphB[n] - (n+1)*SphB[n+2])/(2*n+1)
        end
    else
        for n = 2:lmax+1
            SphB[n] = 0.0
        end
        SphB[begin] = 1.0

        dSphB[1] = 0.0
        for n = 1:lmax
            dSphB[n+1] = (n*SphB[n] - (n+1)*SphB[n+2])/(2*n+1)
        end

        return nothing
    end
end