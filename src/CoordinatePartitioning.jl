module CoordinatePartitioning

# 
using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: eigen, eigvals
using LinearAlgebra: Diagonal

using JuMP: Model as JumpModel
using JuMP: MOI
using JuMP: @variable, @constraint, @build_constraint, @objective
using JuMP: VariableRef
using JuMP: callback_value, callback_node_status, set_attribute, optimize!, objective_value

# Constants
const FLOAT_TOL = 1e-10
const STRATEGIES = ["random"]

# # generic utility functions
include("utils.jl")

# # euclidean distance matrix tools
include("edms.jl")

# # partitioner functions and strategies
include("partitioner.jl")

# # main solver rountine
include("modeler.jl")

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
