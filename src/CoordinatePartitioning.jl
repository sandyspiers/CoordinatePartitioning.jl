module CoordinatePartitioning

# 
using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: eigen, eigvals
using LinearAlgebra: Diagonal

# Constants
const FLOAT_TOL = 1e-10
const STRATEGIES = ["random"]

# # generic utility functions
include("utils.jl")

# # euclidean distance matrix tools
include("edms.jl")

# # partitioner functions and strategies
include("partitioner.jl")

end
