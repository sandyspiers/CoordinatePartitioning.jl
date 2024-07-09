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
