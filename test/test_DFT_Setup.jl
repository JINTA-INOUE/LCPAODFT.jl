include("../src/LCPAODFT.jl")
using .LCPAODFT
using Test


@testset "DFT_Setup Si LDA" begin

    Latvecs = [ 5.10   0.00   5.10;
                0.00   5.10   5.10;
                5.10   5.10   0.00]*"AU"
    atomorb = ["Si7.0-s2p2d1"]
    atomsymbol = ["Si", "Si"]
    atompos = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]*"Frac"
    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 50
    xc_type = "LDA"
    kmesh = (7,7,7)
    filename = "Si"
    verbose = false

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, filename, verbose)
    pao_filebase = basename(dft_setup.pao_file[1])
    pspot_filebase = basename(dft_setup.pspot_file[1])
    
    @test dft_setup.Natom == 2
    @test dft_setup.Nspecies == 1
    @test dft_setup.Nspin == 1
    @test dft_setup.spinsize == 1
    @test dft_setup.atom2spe == [1,1]
    @test dft_setup.Latvecs == [5.10 0.00 5.10; 0.00 5.10 5.10; 5.10 5.10 0.00]
    @test dft_setup.Recvecs ≈ [0.6159985595274104 -0.6159985595274104 0.6159985595274104; -0.6159985595274104 0.6159985595274104 0.6159985595274104; 0.6159985595274104 0.6159985595274104 -0.6159985595274104]
    @test dft_setup.gLatvecs ≈ [0.2125 0.0 0.2125; 0.0 0.2125 0.2125; 0.2125 0.2125 0.0]
    @test dft_setup.gRecvecs ≈ [14.78396542865785 -14.78396542865785 14.78396542865785; -14.78396542865785 14.78396542865785 14.78396542865785; 14.78396542865785 14.78396542865785 -14.78396542865785]
    @test dft_setup.Gxyz ≈ [[0.0,0.0,0.0], [2.55,2.55,2.55]]
    @test dft_setup.Gxyz_frac == [[0.0,0.0,0.0], [0.25,0.25,0.25]]
    @test dft_setup.GridVol ≈ 0.019191406249999998
    @test dft_setup.Grid_Origin == [0.0,0.0,0.0]
    @test dft_setup.Atoms_symbol == ["Si", "Si"]
    @test dft_setup.Atoms_cutoff == [7.0, 7.0]
    @test dft_setup.Atoms_pao == ["Si7.0-s2p2d1"]
    @test dft_setup.system == "Crystal"
    @test dft_setup.Init_Atoms_Nspin == [[2.0,2.0], [2.0,2.0]]
    @test dft_setup.Init_Atoms_Angle == [[0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0]]
    @test dft_setup.SpinPol == "off"
    @test dft_setup.SO_switch == false
    @test dft_setup.xc_type == "LDA"
    @test pao_filebase == "Si7.0.jld2"
    @test pspot_filebase == "Si_CA19.jld2"
    @test dft_setup.Ngrid == (24,24,24)
    @test dft_setup.Mixing_method == "RMM-DIISH"
    @test dft_setup.SCF_criterion == 1e-8
    @test dft_setup.SCF_max == 50
    @test dft_setup.Init_Mixing_weight == 0.3
    @test dft_setup.Min_Mixing_weight == 0.001
    @test dft_setup.Max_Mixing_weight == 0.4
    @test dft_setup.Num_Mixing_Pulay == 5
    @test dft_setup.Start_Pulay_SCF == 6
    @test dft_setup.E_Temp == 300.0
    @test dft_setup.kmesh == (7,7,7)
    @test dft_setup.time_rev == true
    # @test dft_setup.verbose == false
    @test dft_setup.fileout == true
    @test dft_setup.filename == "Si"
end


@testset "DFT_Setup Cdia GGA-PBE" begin

    Latvecs = [ 1.78  1.78  0.00;
                1.78  0.00  1.78;
                0.00  1.78  1.78]*"Ang"
    atomorb = ["C5.0-s2p2d1"]
    atomsymbol = ["C", "C"]
    atompos = [[0.0, 0.0, 0.0], [0.89,  0.89,  0.89]]*"Ang"
    Ecut = 150.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (7,7,7)
    filename = "Cdia"
    verbose = false

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, filename, verbose)
    pao_filebase = basename(dft_setup.pao_file[1])
    pspot_filebase = basename(dft_setup.pspot_file[1])
    
    @test dft_setup.Natom == 2
    @test dft_setup.Nspecies == 1
    @test dft_setup.Nspin == 1
    @test dft_setup.spinsize == 1
    @test dft_setup.atom2spe == [1,1]
    @test dft_setup.Latvecs ≈ [3.363712259708 3.363712259708 0.0; 3.363712259708 0.0 3.363712259708; 0.0 3.363712259708 3.363712259708]
    @test dft_setup.Recvecs ≈ [0.9339659313970309 0.9339659313970304 -0.9339659313970307; 0.9339659313970307 -0.9339659313970307 0.9339659313970309; -0.9339659313970304 0.9339659313970309 0.9339659313970307]
    @test dft_setup.gLatvecs ≈ [0.22424748398053335 0.22424748398053335 0.0; 0.22424748398053335 0.0 0.22424748398053335; 0.0 0.22424748398053335 0.22424748398053335]
    @test dft_setup.gRecvecs ≈ [14.00948897095546 14.00948897095546 -14.00948897095546; 14.00948897095546 -14.00948897095546 14.00948897095546; -14.00948897095546 14.00948897095546 14.00948897095546]
    @test dft_setup.Gxyz ≈ [[0.0,0.0,0.0], [1.681856129854,1.681856129854,1.681856129854]]
    @test dft_setup.Gxyz_frac ≈ [[0.0,0.0,0.0], [0.25,0.25,0.25]]
    @test dft_setup.GridVol ≈ 0.022553436885302317
    @test dft_setup.Grid_Origin ≈ [-2.298536710800467, -2.298536710800467, -2.298536710800467]
    @test dft_setup.Atoms_symbol == ["C", "C"]
    @test dft_setup.Atoms_cutoff == [5.0, 5.0]
    @test dft_setup.Atoms_pao == ["C5.0-s2p2d1"]
    @test dft_setup.system == "Crystal"
    @test dft_setup.Init_Atoms_Nspin == [[2.0,2.0], [2.0,2.0]]
    @test dft_setup.Init_Atoms_Angle == [[0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0]]
    @test dft_setup.SpinPol == "off"
    @test dft_setup.SO_switch == false
    @test dft_setup.xc_type == "GGA-PBE"
    @test pao_filebase == "C5.0.jld2"
    @test pspot_filebase == "C_PBE19.jld2"
    @test dft_setup.Ngrid == (15,15,15)
    @test dft_setup.Mixing_method == "RMM-DIISH"
    @test dft_setup.SCF_criterion == 1e-8
    @test dft_setup.SCF_max == 100
    @test dft_setup.Init_Mixing_weight == 0.3
    @test dft_setup.Min_Mixing_weight == 0.001
    @test dft_setup.Max_Mixing_weight == 0.4
    @test dft_setup.Num_Mixing_Pulay == 5
    @test dft_setup.Start_Pulay_SCF == 6
    @test dft_setup.E_Temp == 300.0
    @test dft_setup.kmesh == (7,7,7)
    @test dft_setup.time_rev == true
    # @test dft_setup.verbose == false
    @test dft_setup.fileout == true
    @test dft_setup.filename == "Cdia"
end


@testset "DFT_Setup BCC Fe GGA-PBE" begin

    Latvecs = [-2.71176  2.71176  2.71176;
                2.71176 -2.71176  2.71176;
                2.71176  2.71176 -2.71176]*"AU"
    atomorb = ["Fe6.0S-s2p2d1"]
    atomsymbol = ["Fe"]
    atompos = [[0.0, 0.0, 0.0]]*"AU"
    Atoms_Nspin = [[8.0, 6.0]]
    Ecut = 200.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (11,11,11)
    SpinPol = "on"
    filename = "BCC_Fe"
    verbose = false

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; SpinPol, Ecut, Atoms_Nspin, SCF_max, SCF_criterion, xc_type, kmesh, filename, verbose)
    pao_filebase = basename(dft_setup.pao_file[1])
    pspot_filebase = basename(dft_setup.pspot_file[1])
    
    @test dft_setup.Natom == 1
    @test dft_setup.Nspecies == 1
    @test dft_setup.Nspin == 2
    @test dft_setup.spinsize == 2
    @test dft_setup.atom2spe == [1]
    @test dft_setup.Latvecs ≈ [-2.71176 2.71176 2.71176; 2.71176 -2.71176 2.71176; 2.71176 2.71176 -2.71176]
    @test dft_setup.Recvecs ≈ [-3.487868498008632e-16 1.158506893526637 1.158506893526637; 1.158506893526637 9.486097791694229e-17 1.158506893526637; 1.158506893526637 1.158506893526637 9.486097791694229e-17]
    @test dft_setup.gLatvecs ≈ [-0.15065333333333333 0.15065333333333333 0.15065333333333333; 0.15065333333333333 -0.15065333333333333 0.15065333333333333; 0.15065333333333333 0.15065333333333333 -0.15065333333333333]
    @test dft_setup.gRecvecs ≈ [0.0 20.853124083479464 20.853124083479464; 20.853124083479464 -0.0 20.853124083479464; 20.853124083479464 20.853124083479464 0.0]
    @test dft_setup.Gxyz == [[0.0,0.0,0.0]]
    @test dft_setup.Gxyz_frac == [[0.0,0.0,0.0]]
    @test dft_setup.GridVol ≈ 0.013677169435486814
    @test dft_setup.Grid_Origin ≈ [-1.2805533333333334, -1.2805533333333334, -1.2805533333333334]
    @test dft_setup.Atoms_symbol == ["Fe"]
    @test dft_setup.Atoms_cutoff == [6.0]
    @test dft_setup.Atoms_pao == ["Fe6.0S-s2p2d1"]
    @test dft_setup.system == "Crystal"
    @test dft_setup.Init_Atoms_Nspin == [[8.0,6.0]]
    @test dft_setup.Init_Atoms_Angle == [[0.0,0.0,0.0,0.0]]
    @test dft_setup.SpinPol == "on"
    @test dft_setup.SO_switch == false
    @test dft_setup.xc_type == "GGA-PBE"
    @test pao_filebase == "Fe6.0S.jld2"
    @test pspot_filebase == "Fe_PBE19S.jld2"
    @test dft_setup.Ngrid == (18,18,18)
    @test dft_setup.Mixing_method == "RMM-DIISH"
    @test dft_setup.SCF_criterion == 1e-8
    @test dft_setup.SCF_max == 100
    @test dft_setup.Init_Mixing_weight == 0.3
    @test dft_setup.Min_Mixing_weight == 0.001
    @test dft_setup.Max_Mixing_weight == 0.4
    @test dft_setup.Num_Mixing_Pulay == 5
    @test dft_setup.Start_Pulay_SCF == 6
    @test dft_setup.E_Temp == 300.0
    @test dft_setup.kmesh == (11,11,11)
    @test dft_setup.time_rev == true
    # @test dft_setup.verbose == false
    @test dft_setup.fileout == true
    @test dft_setup.filename == "BCC_Fe"
end



@testset "DFT_Setup Bi2Se3 with SOC GGA-PBE" begin

    Latvecs = [2.0715000     1.1959811     9.5453333;
              -2.0715000     1.1959811     9.5453333;
               0.0000000    -2.3919622     9.5453333]*"Ang"
    atomorb = ["Bi8.0-s3p2d2f1", "Se7.0-s3p2d2"]
    atomsymbol = ["Bi", "Bi", "Se", "Se", "Se"]
    atompos = [[0.4008, 0.4008, 0.4008],
               [0.5992, 0.5992, 0.5992],
               [0.0000, 0.0000, 0.0000],
               [0.2117, 0.2117, 0.2117],
               [0.7883, 0.7883, 0.7883]]*"FRAC"
    Ecut = 200.0
    system = "Crystal"
    SCF_criterion = 1e-8
    SCF_max = 100
    xc_type = "GGA-PBE"
    kmesh = (7,7,7)
    SpinPol = "nc"
    SO_switch = true
    filename = "Bi2Se3_SOC"
    verbose = false

    dft_setup = DFT_Setup(Latvecs, atomorb, atomsymbol, atompos, system; Ecut, SCF_max, SCF_criterion, xc_type, kmesh, SpinPol, SO_switch, filename, verbose)
    pao_filebase = basename.(dft_setup.pao_file)
    pspot_filebase = basename.(dft_setup.pspot_file)
    
    @test dft_setup.Natom == 5
    @test dft_setup.Nspecies == 2
    @test dft_setup.Nspin == 4
    @test dft_setup.spinsize == 1
    @test dft_setup.atom2spe == [1, 1, 2, 2, 2]
    @test dft_setup.Latvecs ≈ [3.9145673853849 2.2600765665444156 18.038064406859; -3.9145673853849 2.2600765665444156 18.038064406859; 0.0 -4.520153133088831 18.038064406859]
    @test dft_setup.Recvecs ≈ [0.8025389128103861 0.4633460506153244 0.116109747429264; -0.802538912810386 0.4633460506153245 0.11610974742926401; -0.0 -0.9266921012306492 0.11610974742926401]
    @test dft_setup.gLatvecs ≈ [0.13048557951283 0.07533588555148052 0.6012688135619667; -0.13048557951283 0.07533588555148052 0.6012688135619667; 0.0 -0.15067177110296104 0.6012688135619667]
    @test dft_setup.gRecvecs ≈ [24.07616738431159 13.900381518459735 3.4832924228779203; -24.07616738431158 13.900381518459735 3.4832924228779203; -0.0 -27.80076303691947 3.4832924228779203]
    @test dft_setup.Gxyz ≈ [[0.0, 0.0, 21.68896864280726], [0.0, 0.0, 32.42522457776974], [0.0, 0.0, 0.0], [0.0, 0.0, 11.455974704796152], [0.0, 0.0, 42.65821851578085]]
    @test dft_setup.Gxyz_frac == [[0.4008, 0.4008, 0.4008], [0.5992, 0.5992, 0.5992], [0.0, 0.0, 0.0], [0.2117, 0.2117, 0.2117], [0.7883, 0.7883, 0.7883]]
    @test dft_setup.GridVol ≈ 0.0354637245653329
    @test dft_setup.Grid_Origin == [0.0, 0.0, 0.0]
    @test dft_setup.Atoms_symbol == ["Bi", "Bi", "Se", "Se", "Se"]
    @test dft_setup.Atoms_cutoff == [8.0, 8.0, 7.0, 7.0, 7.0]
    @test dft_setup.Atoms_pao == ["Bi8.0-s3p2d2f1", "Se7.0-s3p2d2"]
    @test dft_setup.system == "Crystal"
    @test dft_setup.Init_Atoms_Nspin == [[7.5, 7.5], [7.5, 7.5], [3.0, 3.0], [3.0, 3.0], [3.0, 3.0]]
    @test dft_setup.Init_Atoms_Angle == [[0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0]]
    @test dft_setup.SpinPol == "nc"
    @test dft_setup.SO_switch == true
    @test dft_setup.xc_type == "GGA-PBE"
    @test pao_filebase == ["Bi8.0.jld2", "Se7.0.jld2"]
    @test pspot_filebase == ["Bi_PBE19.jld2", "Se_PBE19.jld2"]
    @test dft_setup.Ngrid == (30, 30, 30)
    @test dft_setup.Mixing_method == "RMM-DIISH"
    @test dft_setup.SCF_criterion == 1e-8
    @test dft_setup.SCF_max == 100
    @test dft_setup.Init_Mixing_weight == 0.3
    @test dft_setup.Min_Mixing_weight == 0.001
    @test dft_setup.Max_Mixing_weight == 0.4
    @test dft_setup.Num_Mixing_Pulay == 5
    @test dft_setup.Start_Pulay_SCF == 6
    @test dft_setup.E_Temp == 300.0
    @test dft_setup.kmesh == (7,7,7)
    @test dft_setup.time_rev == false
    # @test dft_setup.verbose == false
    @test dft_setup.fileout == true
    @test dft_setup.filename == "Bi2Se3_SOC"
end

