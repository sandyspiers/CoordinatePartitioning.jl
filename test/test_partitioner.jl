@testset "partitioner.jl" begin
    # check partition checked is correct
    @test ispartition([1, [2, 4], [3], [5, 6, 7]], 7)
    @test !ispartition([1, [2, 4], [5, 6, 7]], 7)
    @test !ispartition([1, [2, 4], [4], [5, 6, 7]], 7)
    # check evenish partitioner
    for n in 5:10
        for b in 1:10
            @test ispartition(evenish_partition(shuffle(1:n), b), n)
        end
    end
    # check partition strategies
    loc = rand(100, 2) .* 100
    edm = build_edm(loc)
    new_loc, evals = euclid_embed(edm; centered=true)
    coords = size(new_loc)[2]
    for strat in STRATEGIES_ALL
        for num_par in 1:10
            par = partition(evals, num_par, strat)
            @test ispartition(par, coords)
            if strat != "total" && strat != "none"
                @test length(par) <= num_par
            end
            partitioned_edms = build_edms(new_loc, par)
            @test all(isedm.(partitioned_edms))
            @test first(size(partitioned_edms)) == length(par)
            aggregated_edms = first(sum(partitioned_edms; dims=1))
            @test edm ≈ aggregated_edms
        end
    end
end
