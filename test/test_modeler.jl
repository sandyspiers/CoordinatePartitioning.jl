function glover_solve(
    edm::Matrix{T} where {T<:Real}, cardinality::Integer, optimizer::Module
)
    # Instantiate model
    n = size(edm)[1]
    model = JumpModel(optimizer.Optimizer)
    # add variables and cardinality constraint
    @variable(model, 0 <= x[1:n] <= 1, Bin)
    @constraint(model, sum(x) == cardinality)
    # introduce glover variables
    w = @variable(model, w[1:(n - 1)] >= 0)
    # formulate glover linearisation
    @constraint(model, [i = 1:(n - 1)], w[i] <= x[i] * sum(edm[(i + 1):end, i]))
    @constraint(model, [i = 1:(n - 1)], w[i] <= edm[(i + 1):end, i]' * x[(i + 1):end])
    # add glover objective
    @objective(model, Max, sum(w))
    # Solve
    optimize!(model)
    return objective_value(model)
end

@testset "modeler.jl" begin
    # simple test of convergence
    n, s, p = 10, 2, 4
    edm = rand_edm(n, s)
    obj_val = glover_solve(edm, p, GLPK)
    for strat in STRATEGIES_ALL
        for par in 1:n
            val, ncut = solve(edm, p, strat, par, GLPK)
            @test val ≈ obj_val
            @test ncut > 0
        end
    end
end
