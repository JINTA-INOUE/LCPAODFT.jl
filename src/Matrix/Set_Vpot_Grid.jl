function Set_Vpot_Grid!(xc_func::AbstractXC, dVHart_Grid, Vpot_Grid)

    SpinPol = xc_func.SpinPol
    Vxc_Grid = xc_func.Vxc_Grid

    if SpinPol == "off"
        @. Vpot_Grid[begin] = dVHart_Grid + Vxc_Grid[begin]
    elseif SpinPol == "on"
        @. Vpot_Grid[begin] = dVHart_Grid + Vxc_Grid[begin]
        @. Vpot_Grid[end] = dVHart_Grid + Vxc_Grid[end]
    elseif SpinPol == "nc"
        @. Vpot_Grid[1] = dVHart_Grid + Vxc_Grid[1]
        @. Vpot_Grid[2] = dVHart_Grid + Vxc_Grid[2]
        @. Vpot_Grid[3] = deepcopy(Vxc_Grid[3])
        @. Vpot_Grid[4] = deepcopy(Vxc_Grid[4])
    else
        error("please check SpinPol")
    end
end

