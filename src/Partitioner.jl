module Partitioner

using Distances: pairwise
using Distances: Euclidean

using LinearAlgebra: issymmetric, eigen, eigvals
using LinearAlgebra: Diagonal

using ..CoordinatePartitioning: FLOAT_TOL

export !

const STRATEGIES = ["random"]

# Build an edm base on set of locations
# By default, each row is a location
function build_edm(locations, locations_by_row=true)
    if ndims(locations) <= 1
        return pairwise(Euclidean(), locations)
    end
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

# Conduct eigenval decomposition, remove evals below FLOAT_TOL
function nonzero_eigen(mat::Matrix{T}, tol::Real=FLOAT_TOL) where {T<:Real}
    vals, vecs = eigen(mat)
    non_zero_evals = abs.(vals) .>= tol
    return vals[non_zero_evals], vecs[:, non_zero_evals]
end

# Takes a distance matrix and returns a set of locations whose
# sqaured distances matches the given matrix
function euclid_embed(edm::Matrix{T}; centered=false) where {T<:Real}
    gramm = grammian(edm; centered=centered)
    # conduct nonzero eigenvalue decomposition
    vals, vecs = nonzero_eigen(gramm)
    # check that its PSD (otherwise raise error)
    if minimum(vals) < -FLOAT_TOL
        throw(DomainError("The EDM provided is not valid (Grammian not PSD)!"))
    end
    # recreate and return
    loc = vecs * Diagonal(sqrt.(vals))
    # if cenetered, assume they wants the evals as well
    if centered
        return loc, vals
    end
    return loc
end

# Remove a given item from a collention
function remove!(collection, item)
    return deleteat!(collection, findfirst(isequal(item), collection))
end

# Checks that a partition set is valid
function ispartition(partitions, elements::Int)
    remaining_elements = collect(1:elements)
    for element in vcat(partitions...)
        try
            remove!(remaining_elements, element)
        catch
            return false
        end
    end
    if length(remaining_elements) != 0
        return false
    end
    return true
end

# Return a partition set determined by evals and with a given strategy
function partition(evals::Vector{T}, num::Int; strategy::String="all") where {T<:Real}
    # TODO: Provide other strategies
    coords = length(evals)
    if strategy == "random"
        remaining = collect(1:coords)
        poprand!() = popat!(remaining, rand(1:length(remaining)))
        partitions = [[poprand!()] for _ in 1:min(num, coords)]
        p = 1
        while length(remaining) > 0
            push!(partitions[p], poprand!())
            p += 1
            if p > length(partitions)
                p = 1
            end
        end
        return partitions
    end
    throw(ArgumentError("$strategy is not a valid partition strategy"))
end

# takes the provided partition set, and returns a set a squared EDMs
function build_edms(locations::Matrix{T}, partition) where {T<:Real}
    return [build_edm(locations[:, par]) .^ 2 for par in partition]
end

end
