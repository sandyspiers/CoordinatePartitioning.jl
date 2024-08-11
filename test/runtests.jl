using CoordinatePartitioning: CoordinatePartitioning
using CoordinatePartitioning:
    STRATEGIES_ALL,
    build_edm,
    build_edms,
    construct,
    euclid_embed,
    grammian,
    isedm,
    ispartition,
    evenish_partition,
    partition,
    rand_edm,
    solve

using GLPK
using Random: shuffle
using JuMP: Model as JumpModel
using JuMP: @variable, @constraint, @objective
using JuMP: objective_value, optimize!

using Aqua
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(#
        CoordinatePartitioning;
        ambiguities=(broken=true,),
    )
end

include("test_edms.jl")
include("test_partitioner.jl")
include("test_modeler.jl")
