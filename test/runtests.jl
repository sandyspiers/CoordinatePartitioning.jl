using CoordinatePartitioning:
    FLOAT_TOL,
    STRATEGIES,
    build_edm,
    build_edms,
    euclid_embed,
    grammian,
    isedm,
    ispartition,
    partition
using Test

include("test_edms.jl")
include("test_partitioner.jl")
include("test_modeler.jl")
