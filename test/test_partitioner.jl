using CoordinatePartitioning.Partitioner: build_edm, grammian, isedm
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
    end
end
