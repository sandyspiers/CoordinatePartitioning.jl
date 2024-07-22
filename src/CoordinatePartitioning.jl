module CoordinatePartitioning

# exports
# the main export is the solve routine, which solves
# a completly defined euclidean max sum diversity problem
export solve
# we also export all avaliable partition strategies
export STRATEGIES
# for those wanting to solve the problem on their own,
# the following exports will allow this, and are written
# in the steps that should be followed.
export build_edm, euclid_embed, partition, build_edms, construct

# imports
using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: eigen, eigvals
using LinearAlgebra: Diagonal

using Random: shuffle

using JuMP: Model as JumpModel
using JuMP: MOI
using JuMP: @variable, @constraint, @build_constraint, @objective
using JuMP: VariableRef
using JuMP: callback_value, callback_node_status, set_attribute, optimize!, objective_value

# Constants
const FLOAT_TOL = 1e-10
const STRATEGIES = ["random", "greedy", "stratified", "stepped", "total", "none"]

# # generic utility functions
include("utils.jl")

# # euclidean distance matrix tools
include("edms.jl")

# # partitioner functions and strategies
include("partitioner.jl")

# # main solver rountine
include("modeler.jl")

"""
    solve(
    edm::Matrix{T} where {T<:Real},
    cardinality::Integer,
    strategy::String,
    num::Integer,
    optimizer::Module,
)

The main public solve routine.
Takes a Euclidean distance matrics, given cardinality,
partition strategy and desired number of partitions
and solves the problem using the given optimizer package
"""
function solve(
    edm::Matrix{T} where {T<:Real},
    cardinality::Integer,
    strategy::String,
    num::Integer,
    optimizer::Module,
)
    new_loc, evals = euclid_embed(edm; centered=true)
    par = partition(evals, num, strategy)
    edms = build_edms(new_loc, par)
    mdl = construct(edms, cardinality, optimizer)
    optimize!(mdl)
    return objective_value(mdl)
end

end
