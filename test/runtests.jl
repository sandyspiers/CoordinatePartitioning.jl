using CoordinatePartitioning:
    STRATEGIES,
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

using Test

include("test_edms.jl")
include("test_partitioner.jl")
include("test_modeler.jl")
