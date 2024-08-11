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
function partition(evals::Vector{T} where {T<:Real}, bins::Integer, strategy=:all)
    coords = length(evals)
    bins = min(coords, bins)
    strategy = Symbol(strategy)
    if strategy == :random
        return evenish_partition(shuffle(1:coords), bins)
    elseif strategy == :greedy
        return evenish_partition(1:coords, bins)
    elseif strategy == :stratified
        partitions = Vector[Int[] for _ in 1:bins]
        p = 1
        for c in 1:coords
            push!(partitions[p], c)
            p = p == bins ? 1 : p + 1
        end
        return partitions
    elseif strategy == :stepped
        partitions = Vector[Int[] for _ in 1:bins]
        p = 1
        explained = 0
        for c in 1:coords
            push!(partitions[p], c)
            explained += evals[c] / coords
            if explained > 1 / bins && p < bins
                p += 1
                explained = 0
            end
        end
        return partitions
    elseif strategy == :total
        return [[c] for c in 1:coords]
    elseif strategy == :none
        return [collect(1:coords)]
    end
    throw(ArgumentError("$strategy is not a valid partition strategy"))
end

# creates evenish partitions preserving order
function evenish_partition(collection, bins)
    num = length(collection)
    s = Int(floor(num / bins))
    r = num % bins
    p1 = s + 1 == 0 ? [] : Iterators.partition(collection[1:(r * (s + 1))], s + 1)
    p2 = s == 0 ? [] : Iterators.partition(collection[(r * (s + 1) + 1):end], s)
    return vcat(collect(p1), collect(p2))
end

"""
    build_edms(locations::Matrix{T}, partition) where {T<:Real}

Uses a partition set to build a set of **sqaured EDMs** of the given locations,
where the partition sets describes which coordinates to use.
"""
function build_edms(locations::Matrix{T}, partition) where {T<:Real}
    return [build_edm(locations[:, par]) .^ 2 for par in partition]
end
