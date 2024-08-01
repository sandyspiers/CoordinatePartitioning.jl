"""
    build_edm(locations; locations_by_row=true)

Builds a (standard, not squared!) Euclidean distance matrix
between the given locations.
Locations are given by rows by default.
"""
function build_edm(locations; locations_by_row=true)
    if ndims(locations) <= 1
        return pairwise(Euclidean(), locations)
    end
    if locations_by_row
        return pairwise(Euclidean(), locations; dims=1)
    end
    return pairwise(Euclidean(), locations; dims=2)
end

"""
    rand_edm(num::Integer, coords::Integer; axis_limit::Integer=100)

Returns a random EDM of a given number of coordinates.
"""
function rand_edm(num::Integer, coords::Integer; axis_limit::Integer=100)
    locations = rand(num, coords) .* axis_limit
    return build_edm(locations)
end

"""
    grammian(edm::Matrix{T}; centered=false) where {T<:Real}

Returns the Gram matrix of the provided Euclidean distance matrix (EDM).
"""
function grammian(edm::Matrix{T}; centered=false) where {T<:Real}
    n, _ = size(edm)
    magnitudes = edm[1, :] * ones(T, n)'
    gramm = (magnitudes + magnitudes' - edm) ./ 2
    if !centered
        return gramm
    end
    col_sum = repeat(sum(gramm; dims=1), n, 1)
    row_sum = repeat(sum(gramm; dims=2), 1, n)
    return gramm - col_sum ./ n - row_sum ./ n .+ 1 / n^2 * sum(gramm)
end

"""
    isedm(edm::Matrix{T})::Bool where {T<:Real}

Checks if the matrix provided is in fact a EDM.
 1. Is it squared?
 2. Is it symmetric?
 3. Is it's Gramm matrix positive semidefinite?
Returns `true` if so, otherwise `false`.
"""
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

"""
    euclid_embed(edm::Matrix{T}; centered=false) where {T<:Real}

Takes a EDM and returns a set of locations whose **squared Euclidean distance**
matches the distances of the original matrix.
If `centered=true` then also returns the eigenvalues (for use in PCA).
"""
function euclid_embed(edm::Matrix{T}; centered=false) where {T<:Real}
    gramm = grammian(edm; centered=centered)
    # conduct nonzero eigenvalue decomposition
    vals, vecs = nonzero_eigen(gramm, FLOAT_TOL)
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

"""
    re_embed(edm::Matrix{T}; centered=false) where {T<:Real}

Takes a **squared** EDM and returns a set of locations whose **standard Euclidean distance**
matches the distances of the original matrix.
It therefore does the opposite of the above.
If `centered=true` then also returns the eigenvalues (for use in PCA).
"""
function re_embed(edm::Matrix{T}; kwargs...) where {T<:Real}
    return euclid_embed(edm .^ 2; kwargs...)
end
