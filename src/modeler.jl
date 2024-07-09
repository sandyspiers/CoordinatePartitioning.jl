"""
    construct(
    edms::Vector{Matrix{T}}, cardinality::Integer, optimizer::Module
)::JumpModel where {T<:Real}

Takes a Euclidean distance matrices, a given cardinality,
a user provided optimizer and returns a diversity problem JuMP model
with outer approximate-type tangent callbacks
"""
function construct(
    edms::Matrix{T} where {T<:Real}, cardinality::Integer, optimizer::Module
)::JumpModel
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
    edms::Vector{Matrix{T}}, cardinality::Integer, optimizer::Module
)::JumpModel where {T<:Real}
    # instantiate the model
    n = size(edms[1])[1]
    model = JumpModel(optimizer.Optimizer)
    # add variables and cardinality constraint
    @variable(model, 0 <= location_vars[1:n] <= 1, Bin)
    @constraint(model, sum(location_vars) == cardinality)
    # add epigraph variables
    s = length(edms)
    @variable(model, epigraphs[1:s] >= 0)
    # add tangent callback
    _cb(cb_data) = tangent_callback(model, edms, location_vars, epigraphs, cb_data)
    set_attribute(model, MOI.LazyConstraintCallback(), _cb)
    # return the model
    return model
end

function tangent_callback(
    model::JumpModel,
    edms::Vector{Matrix{T}} where {T},
    x::Vector{VariableRef},
    epigraphs::Vector{VariableRef},
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
        end
    end
    return nothing
end
