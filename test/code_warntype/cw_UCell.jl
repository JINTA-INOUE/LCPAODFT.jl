include("../../src/LCPAODFT.jl")
using .LCPAODFT

# Si case

function check_Set_Periodic_code_warntype()
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    CpyCell = 3
    @code_warntype Set_Periodic(Latvecs, CpyCell)
end


function check_Estimate_Trn_System_code_warntype()
    Natom = 2
    CpyCell = 3
    TCpyCell = (2*CpyCell + 1)^3 - 1
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)
    
    @code_warntype Estimate_Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)
end


function check_Trn_System_code_warntype()
    Natom = 2
    CpyCell = 3
    TCpyCell = (2*CpyCell + 1)^3 - 1
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    atv = Vector{Vector{Float64}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv[cell] = zeros(Float64, 3)
    end
    Generation_ATV!(CpyCell, Latvecs, atv)
    
    @code_warntype Trn_System(Natom, Gxyz, Atom_Cut1, atv, TCpyCell)
end


function check_Get_FNAN_code_warntype()
    Natom = 2
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    
    @code_warntype Get_FNAN(Latvecs, Natom, Gxyz, Atom_Cut1)
end


function check_Get_RMI_code_warntype()
    Natom = 2
    CpyCell = 3
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    CpyCell, FNAN, natn, ncn, _ = Get_FNAN(Latvecs, Natom, Gxyz, Atom_Cut1)

    @code_warntype Get_RMI(Natom, CpyCell, FNAN, natn, ncn)
end


function check_Check_system_code_warntype()
    Natom = 2
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    CpyCell, FNAN, natn, ncn, _ = Get_FNAN(Latvecs, Natom, Gxyz, Atom_Cut1)
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv_ijk[cell] = zeros(Int32, 3)
    end
    Generation_ATV_ijk!(CpyCell, atv_ijk)
    @code_warntype Check_system(FNAN, ncn, atv_ijk)
end


function check_Calc_AtomsGrid_code_warntype()
    Natom = 2
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    CpyCell = 3
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    
    @code_warntype LCPAODFT.Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atom_Cut1, Ngrid, Grid_Origin)
    # @code_warntype LCPAODFT.Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atom_Cut1, Int32.(Ngrid), Grid_Origin)
end


function check_Calc_AtomOverlap_Grid_code_warntype()
    Natom = 2
    Gxyz = [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    Latvecs = [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    Atom_Cut1 = [7.0, 7.0]
    CpyCell = 3
    Ngrid = (3,3,3)
    Grid_Origin = [0.0, 0.0, 0.0]
    Total_NumOrbs = [13, 13]
    atv_ijk = Vector{Vector{Int32}}(undef, (2*CpyCell+1)^3)
    for cell = 1:(2*CpyCell+1)^3
        atv_ijk[cell] = zeros(Int32, 3)
    end
    Generation_ATV_ijk!(CpyCell, atv_ijk)

    CpyCell, FNAN, natn, ncn, Dis = Get_FNAN(Latvecs, Natom, Gxyz, Atom_Cut1)
    RMI = Get_RMI(Natom, CpyCell, FNAN, natn, ncn)
    Total_Hsize, MPI_Hsize, MPHks, Nloop, MPI_size, MPI_atom, MPI_FNAN, MPI_natn, MPI_ncn, MPI_Dis, MPI_RMI = split_system_grid(Natom, FNAN, natn, ncn, Dis, RMI, Total_NumOrbs)
    GridN_Atom, GridListAtom, CellListAtom = Calc_AtomsGrid(Latvecs, Natom, CpyCell, Gxyz, Atom_Cut1, Ngrid, Grid_Origin)
    

    @code_warntype LCPAODFT.Calc_AtomOverlap_Grid(CpyCell, MPI_size, MPI_atom, MPI_natn, MPI_ncn, GridN_Atom, GridListAtom, CellListAtom, atv_ijk)
end


# check_Set_Periodic_code_warntype()
# check_Estimate_Trn_System_code_warntype()
# check_Trn_System_code_warntype()
# check_Get_FNAN_code_warntype()
# check_Get_RMI_code_warntype()
# check_Check_system_code_warntype()
# check_Calc_AtomsGrid_code_warntype()