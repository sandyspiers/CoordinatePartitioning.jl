"""
    construct(
    edms::Vector{Matrix{T}}, cardinality::Integer, optimizer::Module
)::JumpModel where {T<:Real}

Takes a Euclidean distance matrices, a given cardinality,
a user provided optimizer and returns a diversity problem JuMP model
with outer approximate-type tangent callbacks
"""
function construct(edms::Matrix{T} where {T<:Real}, cardinality::Integer, optimizer::Module)
    return construct([edms], cardinality, optimizer)
end

"""
    construct(
    edms::Vector{Matrix{T}}, cardinality::Integer, optimizer::Module
)::JumpModel where {T<:Real}

Takes a set of partitioned Euclidean distance matrices, a given cardinality,
a user provided optimizer and returns a diversity problem JuMP model
with outer approximate-type tangent callbacks
"""
function construct(
    edms::Vector{Matrix{T}} where {T<:Real}, cardinality::Integer, optimizer::Module
)
    # instantiate the model
    n = size(edms[1])[1]
    model = JumpModel(optimizer.Optimizer)
    # add variables and cardinality constraint
    @variable(model, 0 <= location_vars[1:n] <= 1, Bin)
    @constraint(model, sum(location_vars) == cardinality)
    if cardinality <= 1
        return model
    end
    # add epigraph variables
    s = length(edms)
    @variable(model, 0 <= epigraphs[c=1:s] <= sum(edms[c]))
    @objective(model, Max, sum(epigraphs))
    # add tangent callback
    num_cuts = Ref(0)
    function _cb(cb_data)
        return tangent_callback(model, edms, location_vars, epigraphs, num_cuts, cb_data)
    end
    set_attribute(model, MOI.LazyConstraintCallback(), _cb)
    # return the model
    return model, num_cuts
end

function tangent_callback(
    model::JumpModel,
    edms::Vector{Matrix{T}} where {T},
    x::Vector{VariableRef},
    epigraphs::Vector{VariableRef},
    num_cuts::Ref{Int},
    callback_data,
)
    if callback_node_status(callback_data, model) != MOI.CALLBACK_NODE_STATUS_INTEGER
        return nothing
    end
    y = callback_value.(callback_data, x)
    for s in eachindex(edms)
        dfy = edms[s] * y
        fy = y' * dfy / 2
        epi_val = callback_value(callback_data, epigraphs[s])
        if epi_val > fy + FLOAT_TOL
            cut = @build_constraint(epigraphs[s] <= -fy + dfy' * x)
            MOI.submit(model, MOI.LazyConstraint(callback_data), cut)
            num_cuts[] += 1
        end
    end
    return nothing
end
