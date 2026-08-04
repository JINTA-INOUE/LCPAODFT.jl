"""
Print the complete timer tree for MPI rank 0 only.  The unlimited IO context
prevents PrettyTables from cropping rows to the terminal height.
"""
function Print_TimerOutput(timer::TimerOutput, comm=MPI.COMM_WORLD)
    nprocs = MPI.Comm_size(comm)
    myrank = MPI.Comm_rank(comm)
    MPI.Barrier(comm)
    if myrank == 0
        # PrettyTables crops vertically according to `displaysize(io)` even
        # when `:limit` is false, so provide an effectively unbounded height.
        output = IOContext(stdout, :limit => false, :displaysize => (1_000_000, max(displaysize(stdout)[2], 120)))
        TimerOutputs.print_timer(output, timer; compact=false)
        flush(output)
    end
    MPI.Barrier(comm)
end
