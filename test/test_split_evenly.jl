using Test
using MPI
include("../src/common/split_evenly.jl")


@testset "split_evenly" begin
    MPI.Init()
    comm = MPI.COMM_WORLD
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)

    N = 100
    my_range = split_evenly(1:N, nprocs)

    @test length(my_range) == nprocs
end
