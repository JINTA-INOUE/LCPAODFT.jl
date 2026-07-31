module LCPAODFT

using LinearAlgebra
using Bessels: besselj
using FFTW
using JLD2
using TimerOutputs
using Spglib
using StaticArrays
using JSON
using SMTPClient
using Printf
using Dates
using MPI
using PrecompileTools


const Ang_to_bohr = 1.8897259886
const kb = 0.00008617251324000000
const eV2Hartree = 27.2113845
export Ang_to_bohr
export kb
export eV2Hartree


# For Gaunt
const S3J_MAX_FACT = 40


# For Gauss Legendre method
const CoarseGL_Mesh = 150
const GL_Mesh = 300
const FineGL_Mesh = 1500
export CoarseGL_Mesh
export GL_Mesh
export FineGL_Mesh


# For Fourier Transform of PAO and VPS
# NkGrid    : number of kmesh for radial Fourier Transform (Ngrid_NormK)
# Nkmax     : max number k for radial Fourier Transform (PAO_Nkmax)
# OneD_Grid : number of rmesh for radial Fourier Transform
const NkGrid = 900
const Nkmax = 60.0
const OneD_Grid = 900
const Radial_kmin = 1e-6
export NkGrid
export Nkmax
export OneD_Grid
export Radial_kmin


# For VNA
const BufferL_ProVNA = 6
const maxM = 12
export BufferL_ProVNA
export maxM


# For ChemP
const max_x = 30.0
export max_x


# For GGA_PBE
const den_min = 1.0e-14
export den_min


# For Set_Lebedev_Grid.jl
const Num_Leb_Grid = 590
export Num_Leb_Grid


# For Boltz
const k_B_SI = 1.380649e-23
const elem_charge_SI = 1.602176565e-19
const hbar_SI = 1.054571726e-34


# For Read_PAO.jl / Read_VPS.jl
const PACKAGE_ROOT = pkgdir(@__MODULE__)
const PAO_File_path = joinpath(PACKAGE_ROOT, "DFT_DATA19", "PAO", "")
const VPS_File_path = joinpath(PACKAGE_ROOT, "DFT_DATA19", "VPS", "")
export PAO_File_path
export VPS_File_path


# For sending_emails
const EmailPass_File_path = joinpath(PACKAGE_ROOT, "smtp.json")
export EmailPass_File_path


const timer = TimerOutput()


include("utils/split_evenly.jl")
export split_evenly


include("DFT_Options.jl")
export default_DFT_Options
export DFT_Options


include("Symmetry/Calc_Symmetry.jl")
export Get_Symmetry_Spglib


include("Atoms/Atom_Data.jl")
include("Atoms/Atom_Lattice_Pos.jl")
include("Atoms/Atom_utils.jl")
include("Atoms/Read_PAO.jl")
include("Atoms/Read_VPS.jl")
export Atompos
export Lattice
export Atom_Core_Charge
export Atom_Znumber
export Get_Atoms_data
export get_ialpha_index
export PAO
export Read_PAO
export Pspot
export Read_VPS


include("DFT_Setup.jl")
export DFT_Setup


include("UCell.jl")
include("KPoints.jl")
include("Electron.jl")
include("Hamiltonian.jl")
export Set_Periodic
export Get_FNAN
export Estimate_Trn_System
export Trn_System
export Get_RMI
export Check_system
export UCell
export System_Grid
export KPoints
export Gen_KPoints
export Electron
export Hamiltonian
export Set_Hamiltonian!
export Calc_MatrixElements_dVH_Vxc_off!
export Calc_MatrixElements_dVH_Vxc_on!
export Calc_MatrixElements_dVH_Vxc_nc!


# exchange calculater
include("XC/XC_Func.jl")
include("XC/LDA_CA.jl")
include("XC/LSDA_CA.jl")
include("XC/XC_PBE.jl")
include("XC/GGA_PBE.jl")
# include("XC/DFT_U_utils.jl")
# include("XC/DFT_U.jl")
export XC_Func
export Calc_Diff_Coef
export LDA_CA
export LSDA_CA
export Set_ex2primitive!
export Set_dDensity_Grid!
export Set_XC_Grid!
export Calc_Vxc_Grid


include("Matrix/Set_OLPpos.jl")
include("Matrix/Set_OLP_Kin.jl")
include("Matrix/Set_Nonlocal.jl")
include("Matrix/Set_ProExpn_VNA.jl")
export Set_OLPpos
export Set_OLP_Kin!
export Set_Nonlocal!
export Set_ProExpn!
export Set_VNA2
export Set_VNA2!
export Set_ProExpn_VNA!


include("Matrix/Set_AdenPCC_Grid.jl")
include("Matrix/Set_Orbitals_Grid.jl")
include("Matrix/Set_Density_Grid.jl")
include("Matrix/Set_Vpot_Grid.jl")
export Set_AdenPCC_Grid
export Set_Orbitals_Grid
export Set_Orbitals_Grid!
export Set_Density_Grid_nonpol!
export Set_Density_Grid_pol!
export Set_Density_Grid_nc!
export Set_Vpot_Grid!


include("Mixing/DFT_Mixing.jl")
include("Mixing/Mixing_H.jl")
export update_NormRD!
export update_Eele!
export update_NormRD_Eele!
export get_Mixing_weight!
export Mixing_H!
export Simple_Mixing_H!
export Pulay_Mixing_H!
export get_metric!
export H_norm
export DFT_Mixing


include("Poisson.jl")
export Solve_Poisson!


include("Cluster_DFT.jl")
include("Crystal_DFT.jl")
export Cluster_DFT!
export Cluster_DFT_Collinear_nonpol!
export Cluster_DFT_Collinear_pol!
export Cluster_DFT_NonCollinear!
export Calc_DM_Cluster_Collinear_nopol!
export Calc_DM_Cluster_Collinear_pol!
export Calc_DM_Cluster_NonCollinear!
export Calc_iDM_Cluster_NonCollinear!
export Crystal_DFT!
export Crystal_DFT_Collinear_nonpol!
export Crystal_DFT_Collinear_pol!
export Crystal_DFT_NonCollinear!
export Calc_fnkCnk!
export Calc_DM_Crystal_Collinear_nopol!
export Calc_DM_Crystal_Collinear_pol!
export Calc_DM_Crystal_NonCollinear!
export Calc_iDM_Crystal_NonCollinear!
export Calc_EDM
export Calc_EDM_Collinear!
export Calc_EDM_NonCollinear!


include("utils/File_utils.jl")
include("utils/Associated_Legendre.jl")
include("utils/Gauss_Legendre.jl")
include("utils/PhiF.jl")
include("utils/KumoF.jl")
include("utils/Dr_KumoF.jl")
include("utils/RadialF.jl")
include("utils/Nonlocal_RadialF.jl")
include("utils/VH_AtomF.jl")
include("utils/Dr_VH_AtomF.jl")
include("utils/VNAF.jl")
include("utils/Int_phi0_phi1.jl")
include("utils/RF_BesselF.jl")
include("utils/Calc_Bessel_Pro00.jl")
include("utils/FT_PAO.jl")
include("utils/FT_NLP.jl")
include("utils/FT_VNA.jl")
include("utils/FT_ProductPAO.jl")
include("utils/FT_ProExpn_VNA.jl")
include("utils/Calc_VNA_Bessel.jl")
include("utils/xyz_to_spherical.jl")
include("utils/EulerAngle_Spin.jl")
include("utils/Set_Comp2Real.jl")
include("utils/Ylm_complex.jl")
include("utils/Ylm_real.jl")
include("utils/Gaunt.jl")
include("utils/Spherical_Besselj.jl")
include("utils/dampingF.jl")
include("utils/Find_CGrid.jl")
include("utils/Generation_ATV.jl")
include("utils/Calc_Spe_Vna.jl")
include("utils/Calc_Gxyz_frac.jl")
include("utils/Calc_Ngrid.jl")
include("utils/Calc_Grid_Origin.jl")
include("utils/HS_matrix.jl")
include("utils/Calc_ChemP.jl")
include("utils/Calc_Band_Energy.jl")
include("utils/Convert_Matrix.jl")
include("utils/wrappers_fft.jl")
include("utils/Calc_dipole_moment.jl")
include("utils/Mulliken_Charge.jl")
include("utils/Set_Lebedev_Grid.jl")
include("utils/Read_restartFile.jl")
export Associated_Legendre
export Associated_Legendre2
export Gauss_Legendre
export Gauss_Legendre_x
export PhiF
export KumoF
export Dr_KumoF
export RadialF
export Nonlocal_RadialF
export VH_AtomF
export Dr_VH_AtomF
export VNAF
export Int_phi0_phi1
export RF_BesselF
export Calc_Bessel_Pro00!
export Calc_ProExpn_VNA!
export FT_PAO!
export FT_NLP
export FT_NLP!
export FT_VNA!
export FT_ProductPAO!
export FT_ProExpn_VNA!
export Calc_VNA_Bessel!
export xyz_to_spherical
export EulerAngle_Spin
export Set_Comp2Real
export Set_Comp2Real!
export Set_NLComp2Real!
export Set_VNAComp2Real
export Ylm_table
export Ylm_complex
export Ylm_real
export calc_Ylm!
export s3j
export Clebsch_Gordan
export Gaunt
export SphericalBesselj
export Calc_SphericalBesselj!
export Calc_SphericalBesselj2!
export dampingF
export deri_dampingF
export Find_ln
export Generation_ATV!
export Generation_ATV_ijk!
export Generation_RATV!
export Spe_VHart_Atom!
export Calc_Spe_VH_Atom
export Calc_Spe_VH_Atom!
export Calc_Spe_Vna
export Calc_Gxyz_frac
export Calc_Ngrid
export Calc_Grid_Origin
export HS_matrix!
export HS_matrix_NC!
export Calc_ChemP
export Calc_Band_Energy!
export Set_DM_Vec2DM
export Set_DM2DM_Vec
export Set_HNL2HNL_Vec
export Set_HVNA2HVNA_Vec
export G_to_R!
export R_to_G!
export Calc_dipole_moment
export Mulliken_Charge
export Mulliken_Charge!
export Set_Lebedev_Grid
export Read_restartFile!


include("Matrix/Set_dOrbitals_Grid.jl")
include("Matrix/Set_OLP_Kinforce.jl")
include("Force.jl")
export Set_dOrbitals_Grid
export Set_dOrbitals_Grid!
export Set_OLP_Kinforce
export Set_OLP_Kinforce!
export PCC_Force
export Kinetic_Force
export Force3_nospin
export Force3_spin
export Force3_nc
export HVNA_Force
export OLP_Force
export HNL_Force
export Core_Force
export EH0_Force
export Exc_Force
export Init_Force
export Force
export Force!


include("Energy.jl")
export Calc_Atomic_Den2
export Init_Energy
export Energy
export Total_Energy!
export Calc_Ecore
export Calc_EH0
export Calc_Ekin
export Calc_Ena
export Calc_Enl
export Calc_EH1
export Calc_EXC1
export Calc_EXC2


include("Optim/GeoOpt_Setup.jl")
include("utils/Extp_Charge.jl")
include("Optim/Geo_Optim.jl")
export Extp_Charge
export GeoOpt_Setup


include("LCPAO_model.jl")
include("utils/Load_LCPAODFT_model.jl")
include("utils/OutData.jl")
export LCPAO_model
export Print_LCPAO_model
export Load_LCPAODFT_model
export WriteFile


include("KSsolve_SCF.jl")
include("DFT.jl")
export KSsolve_SCF
export DFT


# For Cube
include("postprocess/Cube/Print_Density.jl")
include("postprocess/Cube/Calc_psi.jl")
include("postprocess/Cube/Print_psi.jl")
include("postprocess/Cube/Print_CubeData.jl")
include("postprocess/Cube/Print_Cube.jl")
export Print_Cube


# For Closest Wannier Functions
include("postprocess/CWF/Set_MLWF_kgrid.jl")
include("postprocess/CWF/CWF_Setup.jl")
include("postprocess/CWF/CWF_model.jl")
include("postprocess/CWF/Load_CWF_model.jl")
include("postprocess/CWF/Set_OLPexp.jl")
include("postprocess/CWF/CWF_utils.jl")
include("postprocess/CWF/Set_CWF_Guiding_MOs.jl")
include("postprocess/CWF/Calc_Amnk.jl")
include("postprocess/CWF/Calc_Sigmamk_Umnk.jl")
include("postprocess/CWF/Calc_DMfunc.jl")
include("postprocess/CWF/Calc_HmnR.jl")
include("postprocess/CWF/Set_CWF_ExpnCoef.jl")
include("postprocess/CWF/Write_CWF_Cube.jl")
include("postprocess/CWF/Set_CWF_Grid.jl")
include("postprocess/CWF/Calc_WannierCenter.jl")
include("postprocess/CWF/Write_CWF_HmnR.jl")
include("postprocess/CWF/Set_Cnk_work.jl")
include("postprocess/CWF/CWF2Wannier90_utils.jl")
include("postprocess/CWF/Generate_Amnk.jl")
include("postprocess/CWF/Generate_Mmnkb.jl")
include("postprocess/CWF/Generate_eig.jl")
include("postprocess/CWF/Write_win.jl")
include("postprocess/CWF/Write_HmnR_vs_R.jl")
include("postprocess/CWF/CWF2Wannier90.jl")
include("postprocess/CWF/Generate_CWF.jl")
export CWF_Setup
export Calc_WannierCenter
export Write_HmnR_vs_R
export Generate_CWF
export CWF_model


# For Maximally localized Wannier functions
include("postprocess/MLWF/MLWF_Setup.jl")
include("postprocess/MLWF/MLWF_utils.jl")
include("postprocess/MLWF/Set_OLP_WP.jl")
include("postprocess/MLWF/Generate_MLWF.jl")
export MLWF_Setup
export Generate_MLWF


# For Hybrid Wannier functions


# For Boltz
include("postprocess/Boltz/Boltz_Setup.jl")
include("postprocess/Boltz/Boltz_utils.jl")
include("postprocess/Boltz/Calc_Vnk.jl")
include("postprocess/Boltz/Calc_EVec.jl")
include("postprocess/Boltz/Calc_TDF.jl")
include("postprocess/Boltz/Calc_TDF_decomp.jl")
include("postprocess/Boltz/Calc_Sigma.jl")
include("postprocess/Boltz/Calc_Sigma_decomp.jl")
include("postprocess/Boltz/Calc_SigmaS.jl")
include("postprocess/Boltz/Calc_Seebeck.jl")
include("postprocess/Boltz/Calc_Seebeck_decomp.jl")
include("postprocess/Boltz/Calc_Boltz.jl")
include("postprocess/Boltz/Write_TDF.jl")
include("postprocess/Boltz/Write_Sigma.jl")
include("postprocess/Boltz/Write_SigmaS.jl")
include("postprocess/Boltz/Write_Seebeck.jl")
export Boltz_Setup
export Calc_Boltz


# For Band Dispersion
include("postprocess/Band/Write_Band.jl")
include("postprocess/Band/Band_kpath.jl")
export Band_kpath


# For Density of State
include("postprocess/Dos/Tetrahedron_Blochl.jl")
include("postprocess/Dos/Dos_utils.jl")
include("postprocess/Dos/Calc_Enk_Cnk_Dos.jl")
include("postprocess/Dos/Calc_EVec.jl")
include("postprocess/Dos/Calc_DosMain.jl")
include("postprocess/Dos/Calc_PDosMain.jl")
include("postprocess/Dos/Dos.jl")
include("postprocess/Dos/Write_Dos.jl")
export DosMain


# For postprocess
include("postprocess/common/select_model.jl")
include("postprocess/common/Calc_Enk_Cnk.jl")
export Calc_Enk_Cnk
export Calc_Enk_Cnk!


# For result email
include("utils/make_attachment.jl")
include("utils/sending_mail.jl")






# Precompilation block with a basic workflow
@setup_workload begin

    println("Now Precompilation using PrecompileTools.jl")

    # For DFT
    Latvecs = [1.78  1.78  0.00;
               1.78  0.00  1.78;
               0.00  1.78  1.78]*"Ang"
    Atoms_orb = ["C5.0-s2p2d1"]
    Atoms_symbol = ["C", "C"]
    Atoms_pos = [[0.0,0.0,0.0], [0.25,0.25,0.25]]*"Frac"
    system = "Crystal"
    Ecut = 150.0
    SCF_max = 3
    xc_type = "GGA-PBE"
    DFT_kmesh = (3,3,3)
    scf_filename = "Cdia_precompile"
    fileout = true
    verbosity = 1


    # For Band
    filepath = "$scf_filename.jld2"
    kpath = [[0.5,0.75,0.25], [0.5,0.5,0.5], [0.0,0.0,0.0], [0.5,0.5,0.0], [0.5,0.75,0.25], [0.375,0.75,0.375]]
    kname = ["W", "L", "G", "X", "W", "K"]

    # For Dos
    kmesh = (1,1,1)
    Erange = [-25.0,0.0]
    mode = "all"

    # For CWF
    Guide_index = [[1,3,4,5], [1,3,4,5]]
    Dis_Energy = [-26.0, -25.0, 0.0, 16.0]
    kmesh = (3,3,3)
    weight_type = "Poly"
    cwf_filename = "Cdia_AO_Poly"
    CWF_HmnR = true
    CWF_Wannier = true
    CWF2MLWF = true
    CWF_Plot_Cube = [1]
    CWF_Plot_SuperCells = [1,1,1]

    # For Boltz
    TDF_Erange = [-25.0, 0.0]      # eV unit
    kmesh = (1,1,1)
    Temp = 300.0
    decomp = false


    @compile_workload begin

        println("Precompile DFT ...")
        dft_setup = DFT_Setup(Latvecs, Atoms_orb, Atoms_symbol, Atoms_pos, system; 
                              Ecut, SCF_max, xc_type, kmesh=DFT_kmesh, filename=scf_filename, fileout, verbosity)
        DFT(dft_setup)

        println("Precompile Band_kpath ...")
        Band_kpath(filepath, kpath, kname)

        println("Precompile Dos ...")
        DosMain(filepath, kmesh, Erange; mode)

        println("Precompile CWF ...")
        cwf_setup = CWF_Setup(filepath, Guide_index, Dis_Energy; CWF_HmnR, CWF_Wannier, CWF2MLWF, CWF_Plot_Cube, CWF_Plot_SuperCells, filename=cwf_filename, weight_type, Ecut, kmesh)
        Generate_CWF(cwf_setup)
        
        println("Precompile Boltz ...")
        boltz_setup = Boltz_Setup(filepath, kmesh, TDF_Erange, Temp; decomp)
        Calc_Boltz(boltz_setup)
    end
end


end

