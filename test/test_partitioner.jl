@testset "partitioner.jl" begin
    # check partition checked is correct
    @test ispartition([1, [2, 4], [3], [5, 6, 7]], 7)
    @test !ispartition([1, [2, 4], [5, 6, 7]], 7)
    @test !ispartition([1, [2, 4], [4], [5, 6, 7]], 7)
    # check partition strategies
    loc = rand(10, 2) .* 100
    edm = build_edm(loc)
    new_loc, evals = euclid_embed(edm; centered=true)
    coords = size(new_loc)[2]
    for strat in STRATEGIES
        for num_par in 1:10
            par = partition(evals, num_par, strat)
            @test ispartition(par, coords)
            if strat != "all"
                @test length(par) <= num_par
            end
            partitioned_edms = build_edms(new_loc, par)
            @test first(size(partitioned_edms)) == length(par)
            aggregated_edms = first(sum(partitioned_edms; dims=1))
            @test maximum(abs.(edm - aggregated_edms)) <= FLOAT_TOL
        end
    end
end
