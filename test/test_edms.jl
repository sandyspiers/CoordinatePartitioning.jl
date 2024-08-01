@testset "edms.jl" begin
    # test build edm
    locations = [0 1; 0 3]
    edm = build_edm(locations)
    @test edm == [0 2; 2 0]
    # test random edm
    @test isedm(rand_edm(10, 2))
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
        @test isedm(new_edm)
        @test edm ≈ new_edm

        # now check centererd
        new_loc, _ = euclid_embed(edm; centered=true)
        new_edm = build_edm(new_loc) .^ 2
        @test isedm(new_edm)
        @test edm ≈ new_edm

        # check we can re embed
        recovered_loc = re_embed(new_edm)
        @test build_edm(recovered_loc) ≈ edm
    end
end
