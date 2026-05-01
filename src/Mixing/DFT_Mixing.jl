abstract type Mixing end
abstract type Rhor_Mixing <: Mixing end
abstract type Rhok_Mixing <: Mixing end
abstract type DM_Mixing <: Rhor_Mixing end
abstract type Ham_Mixing <: Rhor_Mixing end


mutable struct RMM_DIISH_Mixing <: Ham_Mixing
    Mixing_method::String
    SCF_max::Int32
    is_convergence::Bool
    Conv_iter::Int32
    Latvecs::Matrix{Float64}
    Recvecs::Matrix{Float64}
    Natom::Int32
    Nspin::Int32
    FNAN::Vector{Int32}
    natn::Vector{Vector{Int32}}
    Total_NumOrbs::Vector{Int32}
    Total_Hsize::Int32
    ChemP::Float64
    NormRD::Vector{Float64}
    HisEele::Vector{Float64}
    Density_xyz::Array{ComplexF64,3}
    HisH::Vector{Array{Float64,2}}
    ResH::Vector{Array{Float64,2}}
    Ngrid::Tuple{Int32, Int32, Int32}
end


"""
    DFT_Mixing(...)

Define Mixing method for Self Consistent.

Mandatory arguments:

- `Nspin`: number of spin
- `dft_options`: an instance of `DFT_Options`
- `system_grid`: an instance of `System_Grid`
"""
function DFT_Mixing(Nspin, dft_options::DFT_Options, system_grid::System_Grid)
    
    SCF_max = dft_options.SCF_max
    Mixing_method = dft_options.Mixing_method
    Num_Mixing_Pulay = dft_options.Num_Mixing_Pulay

    Latvecs = system_grid.Latvecs
    Natom = system_grid.Natom
    FNAN = system_grid.FNAN
    natn = system_grid.natn
    Total_NumOrbs = system_grid.Total_NumOrbs
    Ngrid = system_grid.Ngrid
    Total_Hsize = system_grid.Total_Hsize



    NormRD = zeros(Float64, SCF_max+2)
    NormRD[1] = 100.0
    HisEele = zeros(Float64, SCF_max+2)
    
    
    if Mixing_method ∈ ("Simple", "Kerker", "RMM-DIISK")
        error("not support $Mixing_method")
    elseif Mixing_method == "RMM-DIISH"

        Density_xyz = zeros(ComplexF64, Ngrid[1], Ngrid[2], Ngrid[3])
        
        # If it is not DFT+U or constrained DFT for noncollinear spin orientations, 
        # the imaginary part of the Hamiltonian matrix does not change during the SCF steps.
        HisH = Vector{Array{Float64,2}}(undef, Num_Mixing_Pulay+1)
        ResH = Vector{Array{Float64,2}}(undef, Num_Mixing_Pulay+1)
        for his = 1:Num_Mixing_Pulay+1
            HisH[his] = zeros(Float64, Total_Hsize, Nspin)
            ResH[his] = zeros(Float64, Total_Hsize, Nspin)
        end

        return RMM_DIISH_Mixing(Mixing_method, SCF_max, false, 0, Latvecs, 2*pi*inv(Latvecs'),
                                Natom, Nspin, FNAN, natn, Total_NumOrbs, Total_Hsize,
                                0.0,
                                NormRD, HisEele, 
                                Density_xyz, 
                                HisH, ResH,
                                Ngrid)
    end
end


function update_NormRD!( Norm, dft_mixing::Mixing )
    
    SCF_max = dft_mixing.SCF_max
    for i = SCF_max:-1:2
        dft_mixing.NormRD[i] = dft_mixing.NormRD[i-1]
    end

    dft_mixing.NormRD[1] = Norm
end


function update_Eele!(Eele, dft_mixing::Mixing)
    
    SCF_max = dft_mixing.SCF_max
    for i = SCF_max:-1:2
        dft_mixing.HisEele[i] = dft_mixing.HisEele[i-1]
    end

    dft_mixing.HisEele[1] = Eele
end


function update_NormRD_Eele!(Norm, Eele, dft_mixing::Mixing)
    update_NormRD!(Norm, dft_mixing)
    update_Eele!(Eele, dft_mixing)
end


function get_Mixing_weight!(dft_options::DFT_Options, dft_mixing::Mixing)
    
    HisEele = dft_mixing.HisEele
    NormRD = dft_mixing.NormRD
    
    Mixing_weight = dft_options.Mixing_weight
    Min_Weight = dft_options.Min_Mixing_weight
    Max_Mixing_weight = dft_options.Max_Mixing_weight
    SCF_RENZOKU = dft_options.SCF_RENZOKU


    if SCF_RENZOKU == -1
        Max_Weight = Max_Mixing_weight
        dft_options.Max_Mixing_weight2 = Max_Mixing_weight
    else
        Max_Weight = dft_options.Max_Mixing_weight2
    end


    if dft_mixing.SCF_max < 2
        println("SCF max = ", dft_mixing.SCF_max)
        error("please check SCF_max\n")
    end


    
    if (sign(HisEele[1] - HisEele[2]) == sign(HisEele[2] - HisEele[3])
        && NormRD[1] < NormRD[2])

        temp = NormRD[2]/max(NormRD[2]-NormRD[1], 1e-10)*Mixing_weight

        if temp < Max_Weight
            Mixing_weight = ifelse(Min_Weight<temp, temp, Min_Weight)
        else
            Mixing_weight = Max_Weight
            dft_options.SCF_RENZOKU += 1
        end

    elseif (sign(HisEele[1] - HisEele[2]) == sign(HisEele[2] - HisEele[3])
        && NormRD[2] < NormRD[1])

        temp = NormRD[2]/max(NormRD[2]+NormRD[1], 1e-10)*Mixing_weight

        if temp < Max_Weight
            Mixing_weight = ifelse(Min_Weight<temp, temp, Min_Weight)
        else
            Mixing_weight = Max_Weight
        end
        dft_options.SCF_RENZOKU = -1

    elseif (sign(HisEele[1] - HisEele[2]) !== sign(HisEele[2] - HisEele[3])
        && NormRD[1] < NormRD[2])

        temp = NormRD[2]/max(NormRD[2]-NormRD[1], 1e-10)*Mixing_weight

        if temp < Max_Weight
            Mixing_weight = ifelse(Min_Weight<temp, temp, Min_Weight)
        else
            Mixing_weight = Max_Weight
            dft_options.SCF_RENZOKU += 1
        end

    elseif (sign(HisEele[1] - HisEele[2]) != sign(HisEele[2] - HisEele[3])
        && NormRD[1] > NormRD[2])

        temp = NormRD[2]/max(NormRD[2]+NormRD[1], 1e-10)*Mixing_weight
        
        if temp < Max_Weight
            Mixing_weight = ifelse(Min_Weight<temp, temp, Min_Weight)
        else
            Mixing_weight = Max_Weight
        end
        dft_options.SCF_RENZOKU = -1
    end


    dft_options.Mixing_weight = Mixing_weight
end
