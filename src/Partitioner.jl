module Partitioner

using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: issymmetric, eigvals

export !

# Build an edm base on set of locations
# By default, each row is a location
function build_edm(locations, locations_by_row=true)
    if locations_by_row
        return pairwise(Euclidean(), locations; dims=1)
    end
    return pairwise(Euclidean(), locations; dims=2)
end

# Returns the grammian of a given matrix
function grammian(edm::Matrix{T}) where {T<:Real}
    n, _ = size(edm)
    magnitudes = edm[1, :] * ones(T, n)'
    return magnitudes + magnitudes' - edm
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

end
