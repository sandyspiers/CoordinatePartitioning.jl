using CoordinatePartitioning: FLOAT_TOL
using CoordinatePartitioning.Partitioner: STRATEGIES
using CoordinatePartitioning.Partitioner:
    build_edm,
    grammian,
    isedm,
    aggregated_matrix,
    euclid_embed,
    ispartition,
    partition,
    build_edms
using Test

@testset "Partitioner.jl" begin
    @testset "EDM tests" begin
        # test build edm
        locations = [0 1; 0 3]
        edm = build_edm(locations)
        @test edm == [0 2; 2 0]
        # test edmvalid
        for k in 1:5
            locations = rand(10, 2) .* 100
            edm = build_edm(locations)
            @test size(edm) == (10, 10)
            @test isedm(edm)
        end
        for k in 1:5
            bad_edm = rand(10, 9)
            @test !isedm(bad_edm)
            bad_edm = rand(10, 10)
            @test !isedm(bad_edm)
            bad_edm += bad_edm'
            @test !isedm(bad_edm)
        end
        # test aggregated matrix
        @test aggregated_matrix([1 2 2; 3 1 2]) == [9 8 9; 10 9 10]
        # test noncentered euclidean embeddings
        for k in 1:5
            edm = build_edm(rand(10, 2) .* 100)
            @test isedm(edm)
            # start with noncentered
            new_loc = euclid_embed(edm)
            new_edm = build_edm(new_loc) .^ 2
            diff = maximum(abs.(edm - new_edm))
            @test diff < FLOAT_TOL

            # now check centererd
            new_loc, _ = euclid_embed(edm; centered=true)
            new_edm = build_edm(new_loc) .^ 2
            diff = maximum(abs.(edm - new_edm))
            @test diff < FLOAT_TOL
        end
    end

    @testset "Partitioners" begin
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
                par = partition(evals, num_par; strategy=strat)
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
end
