mutable struct DFT_Options
	Mixing_method::String
    SCF_criterion::Float64
    SCF_max::Int64
	Mixing_weight::Float64
    Min_Mixing_weight::Float64
    Max_Mixing_weight::Float64
    Max_Mixing_weight2::Float64
    Num_Mixing_Pulay::Int64
    SCF_RENZOKU::Int64
    Start_Pulay_SCF::Int64
    Extra_CHistory::Int64
    His_Gxyz::Vector{Vector{Float64}}
    crystal_sym::Bool
    time_rev::Bool
end


"""
    dft_options = DFT_Options(...)

Create an instance of `DFT_Options`.

The following is the most commonly used optional arguments:
- `Mixing_method`: which use Mixing method (`Simple`, `Kerker`, `RMM-DIISK`, `Mixing-H`)  
			  	default `Simple`
- `Mixing_weight`     :  scf.Init.Mixing.Weight (0.3)
- `Min_Mixing_weight` :  scf.Min.Mixing.Weight (0.001)
- `Max_Mixing_weight` :  scf.Max.Mixing.Weight (0.4)
- `Num_Mixing_Pulay`  :  scf.Mixing.History (5)
- `Kerker_factor`     :  scf.Kerker.factor
- `Start_Pulay_SCF`   :  scf.Mixing.StartPulay (6)
- `crystal_sym`           :  whether symmetry is used or not.
- `time_rev`          :  time reversel Symmetry
"""
function default_DFT_Options(;
    Mixing_method = "RMM-DIISH",
    SCF_criterion = 1e-6,
    SCF_max = 10,
    Init_Mixing_weight = 0.3,
    Min_Mixing_weight = 0.001,
    Max_Mixing_weight = 0.4,
    Num_Mixing_Pulay = 5,
    Start_Pulay_SCF = 6,
    Extra_CHistory = 3,
    crystal_sym = false,
    time_rev = true)

    His_Gxyz = Vector{Vector{Float64}}(undef, Extra_CHistory)
    for i = 1:Extra_CHistory
        His_Gxyz[i] = zeros(Float64, 3*1)
    end

    return DFT_Options(
        Mixing_method,
        SCF_criterion,
        SCF_max,
        Init_Mixing_weight, Min_Mixing_weight, Max_Mixing_weight, Max_Mixing_weight,
        Num_Mixing_Pulay,
        -1,
        Start_Pulay_SCF,
        3, His_Gxyz, crystal_sym, time_rev)
end

