"""
    AbstractXC - LDA --  XC_LDA
               |
               - LSDA --  XC_LSDA
               |
               - GGA -- XC_GGA_PBE
"""


abstract type AbstractXC end
abstract type LDA  <: AbstractXC end
abstract type LSDA <: AbstractXC end
abstract type GGA  <: AbstractXC end


mutable struct XC_LDA <: LDA
    xc_type::String
    SpinPol::String
    Nspin::Int32
    Ngrid::Tuple{Int32,Int32,Int32}
    Vxc_Grid::Vector{Vector{Float64}}
end


mutable struct XC_LSDA <: LSDA
    xc_type::String
    SpinPol::String
    Nspin::Int32
    Ngrid::Tuple{Int32,Int32,Int32}
    Vxc_Grid::Vector{Vector{Float64}}
end


mutable struct XC_GGA_PBE <: GGA
    xc_type::String
    SpinPol::String
    Nspin::Int32
    Ngrid::Tuple{Int32,Int32,Int32}
    Diff_Coef::Vector{Vector{Float64}}
    Vxc_Grid::Vector{Vector{Float64}}
    dDensity_Grid::Vector{Vector{Vector{Float64}}}
    dEXC_dGD::Vector{Vector{Vector{Float64}}}
end



function XC_Func(xc_type::String, SpinPol::String, Nspin, Ngrid, gLatvecs)

    xc_type = uppercase(xc_type)
	if xc_type ∉ ("LDA", "LSDA", "GGA_PBE")
		println("xc_type is $(xc_type) ")
		error("please check xc_type name")
	end

    SpinPol = lowercase(SpinPol)
    if SpinPol ∉ ("off", "on", "nc")
		println("SpinPol is $(SpinPol) ")
		error("please check SpinPol name")
	end


    if xc_type == "LDA" && SpinPol ∈ ("on", "nc")
        println("xc_type is $(xc_type) ")
        println("SpinPol is $(SpinPol) ")
        error("please check xc_type and SpinPol")
    elseif xc_type == "LSDA" && SpinPol == "off"
        println("xc_type is $(xc_type) ")
        println("SpinPol is $(SpinPol) ")
        error("please check xc_type and SpinPol")
    end

    NN = prod(Ngrid)
    Vxc_Grid = Vector{Vector{Float64}}(undef, Nspin)
	for spin = 1:Nspin
		Vxc_Grid[spin] = zeros(Float64, NN)
	end


    if xc_type == "LDA"

        return XC_LDA(xc_type, SpinPol, Nspin, Ngrid, Vxc_Grid)

    elseif xc_type == "LSDA"

        return XC_LSDA(xc_type, SpinPol, Nspin, Ngrid, Vxc_Grid)

    elseif xc_type == "GGA_PBE"

        # For dDen_Grid, dEXC_dGD calculation
        Diff_Coef = Calc_Diff_Coef(gLatvecs)

        dDen_Grid = Vector{Vector{Vector{Float64}}}(undef, Nspin)
        dEXC_dGD = Vector{Vector{Vector{Float64}}}(undef, Nspin)
        for spin = 1:Nspin
            dDen_Grid[spin] = Vector{Vector{Float64}}(undef, 3)
            dEXC_dGD[spin] = Vector{Vector{Float64}}(undef, 3)
            for xyz = 1:3
                dDen_Grid[spin][xyz] = zeros(Float64, NN)
                dEXC_dGD[spin][xyz] = zeros(Float64, NN)
            end
        end

        return XC_GGA_PBE(xc_type, SpinPol, Nspin, Ngrid, Diff_Coef, Vxc_Grid, dDen_Grid, dEXC_dGD)
    else
        println("now xc type is $xc_type")
        error("please check xc_type")
    end
end


function Set_XC_Grid!(xc_func::XC_LDA, PCCDensity_Grid, Density_Grid)
    @. xc_func.Vxc_Grid[1] = LDA_CA(2*(Density_Grid[1]+PCCDensity_Grid), 1)
end


# fix XC_P_switch = 1
function Set_XC_Grid!(xc_func::XC_LSDA, PCCDensity_Grid, Density_Grid)

    SpinPol = xc_func.SpinPol
    Vxc_Grid = xc_func.Vxc_Grid
    
    NN = prod(xc_func.Ngrid)
    for i = 1:NN
        Vxc_up, Vxc_down = LSDA_CA(Density_Grid[1][i]+PCCDensity_Grid[i], Density_Grid[2][i]+PCCDensity_Grid[i], 1)
        Vxc_Grid[1][i] = Vxc_up
        Vxc_Grid[2][i] = Vxc_down
    end

    
    if SpinPol == "nc"
        for i = 1:NN
            Vxc_up = Vxc_Grid[1][i]
            Vxc_down = Vxc_Grid[2][i]
            tmp0 = 0.5*(Vxc_up + Vxc_down)
            tmp1 = 0.5*(Vxc_up - Vxc_down)
            
            Vxc_Grid[4][i] = -tmp1*sin(Density_Grid[3][i])*sin(Density_Grid[4][i])
            Vxc_Grid[3][i] =  tmp1*sin(Density_Grid[3][i])*cos(Density_Grid[4][i])
            Vxc_Grid[2][i] =  tmp0 - tmp1*cos(Density_Grid[3][i])
            Vxc_Grid[1][i] =  tmp0 + tmp1*cos(Density_Grid[3][i])
        end
    end

    xc_func.Vxc_Grid = Vxc_Grid
end


function Set_XC_Grid!(xc_func::XC_GGA_PBE, PCCDensity_Grid, Density_Grid)

    Ngrid = xc_func.Ngrid
    SpinPol = xc_func.SpinPol
    Diff_Coef = xc_func.Diff_Coef
    dDensity_Grid = xc_func.dDensity_Grid
    dEXC_dGD = xc_func.dEXC_dGD
    Vxc_Grid = xc_func.Vxc_Grid

    Set_dDensity_Grid!(Ngrid, SpinPol, Density_Grid, PCCDensity_Grid, Diff_Coef, dDensity_Grid)
    Set_Vxc_GGA_PBE!(Ngrid, SpinPol, Density_Grid, PCCDensity_Grid, dDensity_Grid, Diff_Coef, dEXC_dGD, Vxc_Grid)
    xc_func.Vxc_Grid = Vxc_Grid
end


# For Energy.jl
function Calc_Vxc_Grid(xc_type::String, SpinPol::String, gtv, Ngrid, PCCDensity_Grid, Density_Grid)

    NN = prod(Ngrid)
    spinmax = ifelse(SpinPol=="off", 1, 2)
    Vxc_Grid = Vector{Vector{Float64}}(undef, spinmax)
	for spin = 1:spinmax
		Vxc_Grid[spin] = zeros(Float64, NN)
	end

    
    if xc_type == "LDA"
        @. Vxc_Grid[1] = LDA_CA(2*(Density_Grid[1]+PCCDensity_Grid), 0)
    elseif xc_type == "LSDA"
        for i = 1:NN
            Vxc_up, Vxc_down = LSDA_CA(Density_Grid[1][i]+PCCDensity_Grid[i], Density_Grid[2][i]+PCCDensity_Grid[i], 0)
            Vxc_Grid[1][i] = Vxc_up
            Vxc_Grid[2][i] = Vxc_down
        end
    elseif xc_type == "GGA_PBE"
        Diff_Coef = Calc_Diff_Coef(gtv)
        _Calc_Vxc_Grid_GGA!(SpinPol, Ngrid, PCCDensity_Grid, Density_Grid, Diff_Coef, Vxc_Grid)
    end

    
    return Vxc_Grid
end



function _Calc_Vxc_Grid_GGA!(SpinPol::String, Ngrid, PCCDensity_Grid, Density_Grid, Diff_Coef, Vxc_Grid)

    NN = prod(Ngrid)
    spinmax = ifelse(SpinPol=="off", 1, 2)

    dDen_Grid = Vector{Vector{Vector{Float64}}}(undef, spinmax)
    for spin = 1:spinmax
        dDen_Grid[spin] = Vector{Vector{Float64}}(undef, 3)
        for xyz = 1:3
            dDen_Grid[spin][xyz] = zeros(Float64, NN)
        end
    end
    Set_dDensity_Grid!(Ngrid, SpinPol, Density_Grid, PCCDensity_Grid, Diff_Coef, dDen_Grid)

    
    if SpinPol == "off"
        _Set_Vxc_GGA_PBE_woSpin!(Ngrid, Density_Grid, PCCDensity_Grid, dDen_Grid, Vxc_Grid)
    elseif SpinPol ∈ ("on", "nc")
        _Set_Vxc_GGA_PBE_wSpin!(Ngrid, Density_Grid, PCCDensity_Grid, dDen_Grid, Vxc_Grid)
    else
        println("SpinPol = $(SpinPol) in XC/XC_Func.jl")
        error("please check SpinPol")
    end
end