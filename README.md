# LCPAODFT.jl
This is a Julia code and modified OpenMX for electronic structure calculations
The original code of OpenMX is available at [https://www.openmx-square.org](https://www.openmx-square.org)

Note that at least Julia 1.10 is required.

## Features
- total energy by band
- local density approximation (LDA, LSDA) and generalized gradient approximation (GGA) to the exchange-correlation potential
- norm-conserving pseudopotentials
- variationally optimized pseudo-atomic basis functions
- fully and scalar relativistic treatment within pseudopotential scheme
non-collinear DFT
- RMM-DIIS for Hamiltonian matrix charge mixing schemes
- dispersion analysis by the band calculation
- density of states (DOS) and projected DOS
- parallel execution by Message Passing Interface (MPI)

## Install
1. LCPAODFT.jl install from https://github.com/JINTA-INOUE/LCPAODFT.jl
2. PAO, VPS file path setting in src/LCPAODFT.jl
    ```
    $ cd LCPAODFT.jl/src
    $ nano LCPAODFT.jl
    ```
    default
    ```
    PAO_File_path = "/Users/user1/.julia/dev/LCPAODFT/DFT_DATA19/PAO/"
    VPS_File_path = "/Users/user1/.julia/dev/LCPAODFT/DFT_DATA19/VPS/"
    ```
    PAO_File_path and VPS_File_path must be specified as absolute paths within Read_PAO.jl and Read_VPS.jl, located under src/Atoms.
3. Please install following packages
    ```
    LinearAlgebra, Bessels, FFTW, JLD2, TimerOutputs, Printf, Dates, MPI
    ```

    MPI.jl requires MPI to be installed. 
4. If you want to use LCPAODFT.jl, you must add include("LCPAODFT.jl") and using .LCPAODFT .

## Usage  
Please check examples.

## Precompile

## Contact
Jin Inoue (Kanazawa University)  
jinoue__at__stu.kanazawa-u.ac.jp
Please replace `__at__` by @.