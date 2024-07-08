module Partitioner

using Base: findlastnot, nonnothingtype_checked
using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: issymmetric, eigen, eigvals
using LinearAlgebra: Diagonal

using ..CoordinatePartitioning: FLOAT_TOL

export !

# Build an edm base on set of locations
# By default, each row is a location
function build_edm(locations, locations_by_row=true)
    if locations_by_row
        return pairwise(Euclidean(), locations; dims=1)
    end
    return pairwise(Euclidean(), locations; dims=2)
end

# takes a matrix, and then return a new ones where every 
# element is given by the row and column sum of original
function aggregated_matrix(mat::Matrix{T}) where {T<:Real}
    nr, nc = size(mat)
    col_sum = repeat(sum(mat; dims=1), nr, 1)
    row_sum = repeat(sum(mat; dims=2), 1, nc)
    return col_sum + row_sum
end

# Returns the grammian of a given matrix
function grammian(edm::Matrix{T}; centered=false) where {T<:Real}
    n, _ = size(edm)
    magnitudes = edm[1, :] * ones(T, n)'
    gramm = (magnitudes + magnitudes' - edm) ./ 2
    if !centered
        return gramm
    end
    return gramm - 1 / n .* aggregated_matrix(gramm) .+ 1 / n^2 * sum(gramm)
end

# Check if given matrix is a valid EDM
function isedm(edm::Matrix{T})::Bool where {T<:Real}
    # check square
    n, n2 = size(edm)
    if n != n2
        # @warn "EDM is not square!"
        return false
    end

    # check symetry up to a tolerance
    if maximum(abs.(edm - edm')) > eps()
        # @warn "EDM is not symmetric!"
        return false
    end

    # check PSD of grammian
    min_eval = minimum(eigvals(grammian(edm)))
    if min_eval < -eps()
        # @warn "Grammian is not PSD! Min eva: $min_eval"
        return false
    end

    # Must be an edm
    return true
end

# Takes a distance matrix and returns a set of locations whose
# sqaured distances matches the given matrix
function euclid_embed(edm::Matrix{T}; centered=false) where {T<:Real}
    gramm = grammian(edm; centered=centered)
    vals, vecs = eigen(gramm)
    # get nonzero evals (coordinates)
    non_zero_evals = abs.(vals) .>= FLOAT_TOL
    @info non_zero_evals
    vals = vals[non_zero_evals]
    vecs = vecs[:, non_zero_evals]
    # recreate and return
    loc = vecs * Diagonal(sqrt.(max.(vals, 0)))
    return loc
end

end
