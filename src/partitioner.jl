"""
    ispartition(partitions, elements::Integer)

Is the provided partition set a valid set (of depths at most 1)
of the number of elements?
"""
function ispartition(partitions, elements::Integer)::Bool
    # keep track of remaining elements
    remaining_elements = collect(1:elements)
    # flatten the partitions (should only be depths at most 1)
    # and iterate over
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
function partition(evals::Vector{T} where {T<:Real}, num::Integer, strategy::String="all")
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

"""
    build_edms(locations::Matrix{T}, partition) where {T<:Real}

Uses a partition set to build a set of **sqaured EDMs** of the given locations,
where the partition sets describes which coordinates to use.
"""
function build_edms(locations::Matrix{T}, partition) where {T<:Real}
    return [build_edm(locations[:, par]) .^ 2 for par in partition]
end
