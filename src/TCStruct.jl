"""
    TriangulationCache(function_oracle; kwargs...)

Create an adaptive sampling object over a rectangular window. The Delaunay
triangulation stores the mesh, `function_values` stores
sampled values by vertex index, and `incomplete_triangles` stores only the triangles
still needing refinement.

`function_oracle` may be either batched, `f(points) -> values`, or single-point,
`f(point)` / `f(x, y)`. Batched oracles are preferred and are preserved when
detected.

Keyword arguments:
- `xlims`, `ylims`: domain limits, default `[-1, 1]`.
- `resolution`: initial mesh oracle-call target, default `1000`.
- `strategy`: one of `:random`, `:sierpinski`, `:barycenter`; default `:sierpinski`.
- `min_refinement_area`: finite, nonnegative minimum incomplete-triangle area,
  normalized by the current window area; defaults to `1e-5`, and refinement
  skips triangles at or below this threshold.
- `is_complete`: default `nothing` selects the built-in rule; otherwise a custom
  completeness predicate `(vertices, values; kwargs...)`,
  where `vertices` is an `NTuple{3,NTuple{2,Float64}}` and `values` contains the
  three corresponding oracle values.
- `verbose`: whether to print progress, default `false`.
- `batched`: `nothing` detects the call form, `false` calls the oracle once per
  point, and `true` requires a batched oracle.
- `input_dimension`, `near`, `plane_points`, `zoomer`, `rng`: optional affine
  two-plane slice for functions with more than two inputs. Defaults are
  `input_dimension=nothing`, `near=nothing`, `plane_points=nothing`,
  `zoomer=1.0`, and `rng=Random.default_rng()`.

The value `:wildcard` is special in the default completeness rule: it is treated
as equal to every other value. Non-real values are handled as discrete
categories; real values use a tolerance-based rule.

Important invariants:
- `function_values[i]` is the oracle value at Delaunay vertex `i`.
- `point_indices` maps a point coordinate key `(x, y)` to its Delaunay vertex
  index.
- `incomplete_triangles` stores sorted `TriangleKey`s, so triangle identity is
  independent of vertex order.

`covered_windows` records rectangular viewports that have already been seeded,
so zooming back out to a previously viewed region can redraw without spending
new oracle calls.

`parameter_slice` stores the affine coordinate map for caches built from HC
or ordinary-function slices, or `nothing` for unsliced functions. Call
`TC.parameter_slice([u, v])` to map plot coordinates to the original
input coordinates. [`retrieve_witnesses`](@ref) applies this map automatically.
"""
mutable struct TriangulationCache
    # Batched oracle for f: [(x1,y1)...(xn,yn)]->[f(x1,y1),...,f(xn,yn)]
    function_oracle::Function
    # Map plot coordinates back to the original system parameters, when present.
    parameter_slice::Union{Nothing,AffineParameterSlice}
    # Triangulation type inside DelaunayTriangulation.
    triangulation::Triangulation
    # Function values in lock-step with the point order in the triangulation.
    function_values::FunctionValues

    # Quick lookup from point coordinates to triangulation vertex indices.
    point_indices::Dict{PointKey,Int64}
    incomplete_triangles::Set{TriangleKey}

    # Completeness predicate.
    is_complete::Function
    # Refinement strategy.
    strategy::Symbol

    total_oracle_calls::Int64 #record keeping

    oracle_budget::Union{Nothing,Int64} #user bound
    min_refinement_area::Float64 #user bound

    xlims::Vector{Float64}#user visualization bound
    ylims::Vector{Float64}
    covered_windows::Vector{NTuple{4,Float64}} #keeps track of seen windows
    plot_value_order::Dict{Bool,Vector{Any}}

    # Verbose variable used mostly for debugging.
    verbose::Bool
end

#######################
# Display
#######################

function Base.show(io::IO, TC::TriangulationCache)
    n_incomplete = length(incomplete_triangle_keys(TC))
    n_triangles = count_solid_triangles(TC)
    counts = (length(function_values(TC)), n_triangles)
    count_width = maximum(length, string.(counts))
    print(io, "TriangulationCache with:\n",
        "  ", lpad(string(counts[1]), count_width), " function values\n",
        "  ", lpad(string(counts[2]), count_width), " (", n_incomplete,
        ") triangles (incomplete)")
end

#######################
# Getters
#######################

function_oracle(TC::TriangulationCache) = TC.function_oracle
function_values(TC::TriangulationCache) = TC.function_values
triangulation(TC::TriangulationCache) = TC.triangulation
point_indices(TC::TriangulationCache) = TC.point_indices
incomplete_triangle_keys(TC::TriangulationCache) = TC.incomplete_triangles
incomplete_triangles(TC::TriangulationCache) = [collect(key) for key in incomplete_triangle_keys(TC)]

input_points(TC::TriangulationCache) = [collect(get_point(triangulation(TC), i)) for i in 1:num_points(triangulation(TC))]
output_values(TC::TriangulationCache) = function_values(TC)
is_discrete(TC::TriangulationCache) = is_discrete(function_values(TC))
dimension(::TriangulationCache) = 2
strategy(TC::TriangulationCache) = TC.strategy
is_verbose(TC::TriangulationCache) = TC.verbose
remaining_oracle_budget(TC::TriangulationCache) = TC.oracle_budget

"""
    retrieve_witnesses(TC::TriangulationCache)

Return a vector of parameter vectors, one sampled witness for each distinct
non-`:wildcard` value in `TC`, in the order those values first appear in
`function_values(TC)`. Uses all cached samples, including outside the current
window, without evaluating the oracle again.

For a cache built from a sliced evaluator, witnesses are in the original
parameter coordinates. HC system parameters follow the order of
`HomotopyContinuation.parameters(F)`. Unsliced caches return two-dimensional
plot coordinates. An all-wildcard cache returns an empty vector.

Throws `ArgumentError` when the cached values are classified as continuous.
Classification ignores `:wildcard`: any non-real value makes the cache discrete;
otherwise, fewer than 50 distinct values are required for a discrete cache.
"""
function retrieve_witnesses(TC::TriangulationCache)
    is_discrete(TC) || throw(ArgumentError("retrieve_witnesses requires a discrete TriangulationCache."))
    witnesses = Vector{Float64}[]
    seen = Set{Any}()
    for (i, value) in pairs(function_values(TC))
        (is_wildcard_value(value) || value in seen) && continue
        point = collect(get_point(triangulation(TC), i))
        push!(witnesses, TC.parameter_slice === nothing ? point : TC.parameter_slice(point))
        push!(seen, value)
    end
    return witnesses
end

#######################
# Triangle Helpers
#######################

triangle_indices(T) = Int64[triangle_vertices(T)...]
canonical_triangle(triangle::Vector{Int64}) = Tuple(sort(triangle))::TriangleKey
count_solid_triangles(TC::TriangulationCache) = count(_ -> true, each_solid_triangle(triangulation(TC)))

#######################
# Triangle Queries
#######################

function complete_triangles(TC::TriangulationCache)
    incomplete = incomplete_triangle_keys(TC)
    triangles = TriangleList()
    for T in each_solid_triangle(triangulation(TC))
        triangle = triangle_indices(T)
        canonical_triangle(triangle) in incomplete || push!(triangles, triangle)
    end
    return triangles
end

#######################
# Completeness
#######################

"""
    is_complete(triangle, TC::TriangulationCache; kwargs...)

Evaluate `TC`'s stored completeness predicate on a triangle represented by
vertex indices. The stored predicate receives `(vertices, values)`, where
`vertices` are the triangle's coordinate tuples and `values` are the
corresponding function values.
"""
function is_complete(triangle::Vector{Int64}, TC::TriangulationCache; kwargs...)
    vertices = ntuple(i -> point_key(get_point(triangulation(TC), triangle[i])), Val(3))
    values = function_values(TC)[triangle]
    return TC.is_complete(vertices, values; kwargs...)
end

function recompute_incomplete_triangles!(TC::TriangulationCache)
    empty!(TC.incomplete_triangles)
    for T in each_solid_triangle(triangulation(TC))
        triangle = triangle_indices(T)
        is_complete(triangle, TC) || push!(TC.incomplete_triangles, canonical_triangle(triangle))
    end
    return TC
end
