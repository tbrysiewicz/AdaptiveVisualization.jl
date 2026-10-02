# Shared helpers for adaptive visualization.

#####################################
# Type Aliases
#####################################

const FunctionValues = Vector{Any}
const TriangleList = Vector{Vector{Int64}}
const PointKey = NTuple{2,Float64}
const TriangleKey = NTuple{3,Int64}

#####################################
#  Defaults
#####################################

const TRIANGULATION_CACHE_STRATEGIES = (:random, :sierpinski, :barycenter)
const TRIANGULATION_CACHE_DEFAULT_TOTAL_RESOLUTION = 1000
const TRIANGULATION_CACHE_DEFAULT_MIN_REFINEMENT_AREA = 1e-5
const TRIANGULATION_CACHE_DEFAULT_LEGEND_LIMIT = 20

###################################
# Point and Cache Helpers
###################################

# Convert a point to a cache key.
function point_key(point)::PointKey
    length(point) == 2 || error("TriangulationCache points must have dimension 2.")
    return (Float64(point[1]), Float64(point[2]))
end

# Build point-to-index lookup.
function point_index_dict(points::AbstractVector)
    indices = Dict{PointKey,Int64}()
    for (i, point) in pairs(points)
        key = point_key(point)
        haskey(indices, key) || (indices[key] = i)
    end
    return indices
end

#############################
# Oracle Adapters
#############################

# Keep the affine map alongside slice-based evaluators, including when they are
# passed directly to TriangulationCache instead of through visualize(System).
struct AffineParameterSlice
    center::Vector{Float64}
    basis::Matrix{Float64}
end

(slice::AffineParameterSlice)(point) = slice.center + slice.basis * collect(point_key(point))

function checked_slice_point(point, n, name)
    (point isa AbstractVector || point isa Tuple) && length(point) == n ||
        throw(ArgumentError("$name must have length $n."))
    all(x -> x isa Real && isfinite(x), point) ||
        throw(ArgumentError("$name must contain finite real values."))
    result = Float64.(collect(point))
    all(isfinite, result) || throw(ArgumentError("$name must contain finite Float64 values."))
    return result
end

function affine_parameter_slice(n::Integer; near=nothing, plane_points=nothing,
        zoomer=1.0, rng::AbstractRNG=default_rng())
    n >= 2 || throw(ArgumentError("A parameter slice needs at least two dimensions."))
    zoomer isa Real && isfinite(zoomer) && zoomer > 0 ||
        throw(ArgumentError("zoomer must be finite and positive."))
    if plane_points === nothing
        center = near === nothing ? zeros(n) : checked_slice_point(near, n, "near")
        basis = n == 2 ? Matrix{Float64}(I, 2, 2) : Matrix(qr(randn(rng, n, 2)).Q)[:, 1:2]
    else
        (plane_points isa AbstractVector || plane_points isa Tuple) && length(plane_points) == 3 ||
            throw(ArgumentError("plane_points must contain exactly three points [p,q,r]."))
        p, q, r = [checked_slice_point(point, n, "plane_points[$i]")
                   for (i, point) in enumerate(plane_points)]
        center, basis = p, hcat(q - p, r - p)
    end
    basis *= Float64(zoomer)
    all(isfinite, basis) && rank(basis) == 2 ||
        throw(ArgumentError("plane_points and zoomer must define two finite, independent directions."))
    return AffineParameterSlice(center, basis)
end

struct SlicedOracle{F<:Function} <: Function
    evaluate::F
    parameter_slice::AffineParameterSlice
    discrete::Union{Nothing,Bool}
end

SlicedOracle(evaluate::Function, slice::AffineParameterSlice) = SlicedOracle(evaluate, slice, nothing)
oracle_discrete(::Function) = nothing
oracle_discrete(oracle::SlicedOracle) = oracle.discrete

(oracle::SlicedOracle)(points) = oracle.evaluate(points)
oracle_parameter_slice(::Function) = nothing
oracle_parameter_slice(oracle::SlicedOracle) = oracle.parameter_slice

# Validate a batched oracle.
function checked_batch_oracle(function_oracle::Function)
    function checked_oracle(points)
        values = function_oracle(points)
        if !(values isa AbstractVector) || length(values) != length(points)
            error("Batched function oracle must return one output value for each input point.")
        end
        return collect(values)
    end
    return function_oracle isa SlicedOracle ?
        SlicedOracle(checked_oracle, function_oracle.parameter_slice, oracle_discrete(function_oracle)) :
        checked_oracle
end

# Convert a single-point oracle.
single_to_batched_oracle(function_oracle::Function) = points -> map(p -> function_oracle(p...), points)

#############################################
# Value Classification and Completeness
#############################################

# Detect real values.
is_real_value(value) = value isa Real

# Detect wildcard values.
is_wildcard_value(value) = value === :wildcard

# Drop wildcard values.
non_wildcard_values(values) = filter(!is_wildcard_value, values)

"""
    is_discrete(function_values)

Heuristically decide whether function values are discrete. Non-real values are
always treated as discrete; numeric values use a small-cardinality heuristic.
"""
function is_discrete(function_values::AbstractVector)
    values = non_wildcard_values(function_values)
    #if any value is not a real value or wildcard, declare 'discrete'
    any(value -> !is_real_value(value), values) && return true
    return length(unique(values)) < 50
end

# Check value agreement.
function values_are_complete(values::AbstractVector; tol = 0.0)
    vals = non_wildcard_values(values)
    isempty(vals) && return true
    all(is_real_value, vals) || return all(==(first(vals)), vals)
    vertex_function_values = sort(vals)
    return (vertex_function_values[end] - vertex_function_values[1]) <= tol
end

function validate_discrete(discrete)
    discrete === nothing || discrete isa Bool ||
        throw(ArgumentError("discrete must be true, false, or nothing."))
    return discrete
end

function resolve_discrete(oracle::Function, discrete)
    validate_discrete(discrete)
    return discrete === nothing ? oracle_discrete(oracle) : discrete
end

# Build the default triangle predicate.
function default_is_complete(function_values::AbstractVector; discrete=nothing)
    validate_discrete(discrete)
    tol = 0.0
    if !(discrete === nothing ? is_discrete(function_values) : discrete)
        values = non_wildcard_values(function_values)
        all(is_real_value, values) || error("Non-real function values must be handled as discrete values.")
        isempty(values) || (tol = (maximum(values) - minimum(values)) / 16)
    end
    return (vertices, values::AbstractVector; kwargs...) -> values_are_complete(values; tol=tol, kwargs...)
end

# Average numeric non-wildcard values.
function numeric_mean_or_nothing(values)
    numeric_values = Float64.(filter(is_real_value, non_wildcard_values(values)))
    isempty(numeric_values) && return nothing
    return sum(numeric_values) / length(numeric_values)
end

############################
# Geometry Helpers
############################

# Compute triangle area.
area_of_triangle(P::AbstractVector)::Float64 = 0.5 * abs((P[2][1] - P[1][1]) * (P[3][2] - P[1][2]) - (P[3][1] - P[1][1]) * (P[2][2] - P[1][2]))

#############################
# Console Reporting
#############################

# Print refinement totals.
function print_refinement_summary(resolution_used::Int, skipped_small_triangles::Int)
    println("Resolution used:", resolution_used, ". Small triangles skipped:", skipped_small_triangles)
end
