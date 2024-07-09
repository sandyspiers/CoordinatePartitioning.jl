"""
    nonzero_eigen(mat::Matrix{T} where {T<:Real}, tol::Real)

conducts eigenvalue decompostion on the given matrix, returning
only the evals (and associated evecs) that are greater than a 
given tolerance
"""
function nonzero_eigen(mat::Matrix{T} where {T<:Real}, tol::Real)
    vals, vecs = eigen(mat)
    non_zero_evals = abs.(vals) .>= tol
    return vals[non_zero_evals], vecs[:, non_zero_evals]
end

"""
    remove!(collection, item)

Removes the first instance of `item` from the `collection`.
"""
function remove!(collection, item)
    return deleteat!(collection, findfirst(isequal(item), collection))
end
