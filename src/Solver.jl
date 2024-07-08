module Solver

# should be given an optimiser by the user, but it should use solver independent callbacks
# i.e. one of GLPK, Gurobi, CPLEX
# For test cases, use GLPK
# For implementation, expect Gurobi

function solve(edms::Vector{Matrix{T}}, p::Int) where {T<:Real}
    return nothing
end

end
