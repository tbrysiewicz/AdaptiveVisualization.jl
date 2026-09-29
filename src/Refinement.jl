# Refinement strategies and refinement loops.

############################
# Triangle Geometry
############################

# Compute the area of a cached triangle.
function area_of_triangle(TC::TriangulationCache, triangle::Vector{Int64})
    return area_of_triangle(triangle_points(TC, triangle))
end

# Scale the minimum area threshold.
scaled_min_refinement_area(TC::TriangulationCache) = TC.min_refinement_area * bounding_box_area(TC)

# Get the points of a cached triangle.
function triangle_points(TC::TriangulationCache, triangle::Vector{Int64})
    return [get_point(triangulation(TC), i) for i in triangle]
end

############################
# Refinement Strategies
############################

# Compute the triangle barycenter.
function triangle_barycenter(points)
    n = length(points)
    return [sum(p[1] for p in points) / n, sum(p[2] for p in points) / n]
end

# Sample a random triangle point.
function random_triangle_point(points)
    a, b, c = rand(), rand(), rand()
    s = a + b + c
    weights = (a / s, b / s, c / s)
    return [
        sum(weights[i] * points[i][1] for i in 1:3),
        sum(weights[i] * points[i][2] for i in 1:3),
    ]
end

# Compute edge-midpoint refinement points.
function sierpinski_triangle_points(points)
    return [
        [(points[1][1] + points[2][1]) / 2, (points[1][2] + points[2][2]) / 2],
        [(points[2][1] + points[3][1]) / 2, (points[2][2] + points[3][2]) / 2],
        [(points[3][1] + points[1][1]) / 2, (points[3][2] + points[1][2]) / 2],
    ]
end

# Select refinement points for a triangle.
function refinement_points(TC::TriangulationCache, triangle::Vector{Int64})
    points = triangle_points(TC, triangle)
    if strategy(TC) == :random
        return [random_triangle_point(points)]
    elseif strategy(TC) == :sierpinski
        return sierpinski_triangle_points(points)
    elseif strategy(TC) == :barycenter
        return [triangle_barycenter(points)]
    end
    error("Invalid strategy $(strategy(TC)). Use one of $(TRIANGULATION_CACHE_STRATEGIES).")
end

############################
# Candidate Selection
############################

# Keep candidate points that are not cached or already queued.
# Mutates `queued_keys` so later triangles in the pass skip duplicates.
function uncached_refinement_points(TC::TriangulationCache, points, queued_keys::Set{PointKey})
    new_points = Vector{Vector{Float64}}()
    for p in points
        key = point_key(p)
        if !haskey(point_indices(TC), key) && !(key in queued_keys)
            push!(new_points, Float64.(p))
            push!(queued_keys, key)
        end
    end
    return new_points
end

# Check whether a triangle touches the window.
function triangle_intersects_window(TC::TriangulationCache, triangle::Vector{Int64})
    points = triangle_points(TC, triangle)
    xs = map(p -> p[1], points)
    ys = map(p -> p[2], points)
    return maximum(xs) >= TC.xlims[1] &&
        minimum(xs) <= TC.xlims[2] &&
        maximum(ys) >= TC.ylims[1] &&
        minimum(ys) <= TC.ylims[2]
end

# List incomplete triangles worth refining.
function candidate_triangles(TC::TriangulationCache)
    min_area = scaled_min_refinement_area(TC)
    triangle_areas = [(T, area_of_triangle(TC, T)) for T in incomplete_triangles(TC) if triangle_intersects_window(TC, T)]
    selected = [item for item in triangle_areas if item[2] > min_area]
    return first.(sort(selected; by=last, rev=true))
end

# Count incomplete triangles below the area threshold.
function count_small_refinement_triangles(TC::TriangulationCache)
    min_area = scaled_min_refinement_area(TC)
    return count(T -> area_of_triangle(TC, T) <= min_area, incomplete_triangles(TC))
end

# Collect new points within the budget, window, and minimum area cutoff.
function collect_refinement_points(TC::TriangulationCache, budget::Union{Nothing,Int})
    budget === 0 && return Vector{Vector{Float64}}()
    new_points = Vector{Vector{Float64}}()
    queued_keys = Set{PointKey}()
    for T in candidate_triangles(TC)
        candidates = uncached_refinement_points(TC, refinement_points(TC, T), queued_keys)
        isempty(candidates) && continue
        allowed = budget === nothing ? length(candidates) : min(length(candidates), budget - length(new_points))
        append!(new_points, candidates[1:allowed])
        budget !== nothing && length(new_points) >= budget && break
    end
    return new_points
end

############################
# Point Insertion
############################

# Evaluate and insert new points.
function evaluate_and_insert_points!(TC::TriangulationCache, points::Vector{Vector{Float64}}; verbose=is_verbose(TC))
    isempty(points) && return 0
    verbose && println("Evaluating oracle on ", length(points), " new point(s).")
    #evaluate at all the points
    values = function_oracle(TC)(points)
    length(values) == length(points) || error("Did not evaluate at each parameter.")
    inserted = 0
    #put that information into TC
    for i in eachindex(points)
        point = points[i]
        key = point_key(point)
        haskey(point_indices(TC), key) && continue
        #This adds the point to the DelaunayTriangulation object in a way
        #  which locally changes the triangulation.
        add_point!(triangulation(TC), point[1], point[2])
        index = num_points(triangulation(TC))
        push!(TC.function_values, values[i])
        index == length(function_values(TC)) || error("Function values and Delaunay point indices are out of sync.")
        TC.point_indices[key] = index
        inserted += 1
    end
    verbose && println("Delaunay update: inserted ", inserted, " point(s) incrementally with local add_point! updates; no global retriangulation.")
    #recompute incomplete_triangles
    recompute_incomplete_triangles!(TC)
    return inserted
end

# Spend oracle calls from the budget.
function spend_oracle_budget!(TC::TriangulationCache, amount::Integer)
    TC.oracle_budget === nothing || (TC.oracle_budget = max(0, TC.oracle_budget - amount))
    return TC.oracle_budget
end

# Resolve the budget for this refinement call.
function effective_refinement_budget(TC::TriangulationCache, budget)
    budget === nothing || budget isa Integer || error("Refinement budget must be an integer.")
    budget === nothing || budget >= 0 || error("Refinement budget must be nonnegative.")
    return budget === nothing ? remaining_oracle_budget(TC) : budget
end

############################
# Refinement Loops
############################

# Run one refinement pass.
function refine_one_pass!(TC::TriangulationCache; budget=nothing, verbose=is_verbose(TC))
    pass_budget = effective_refinement_budget(TC, budget)
    points = collect_refinement_points(TC, pass_budget)
    isempty(points) && (print_refinement_summary(0, count_small_refinement_triangles(TC)); return 0)
    verbose && println("Refinement pass: selected ", length(points), " new point(s) from incomplete triangles.")
    inserted = evaluate_and_insert_points!(TC, points; verbose=verbose)
    TC.total_oracle_calls += inserted
    budget === nothing && spend_oracle_budget!(TC, inserted)
    print_refinement_summary(inserted, count_small_refinement_triangles(TC))
    return inserted
end

# Refine until the minimum area cutoff stops progress.
function refine_to_min_area!(TC::TriangulationCache, min_refinement_area::Real; budget=nothing, verbose=is_verbose(TC))
    budget === nothing || budget isa Integer || error("Refinement budget must be an integer.")
    budget === nothing || budget >= 0 || error("Refinement budget must be nonnegative.")
    resolved_min_refinement_area = Float64(min_refinement_area)
    isfinite(resolved_min_refinement_area) || error("min_refinement_area must be finite.")
    resolved_min_refinement_area > 0 || error("Iterative min_refinement_area must be positive.")
    TC.min_refinement_area = resolved_min_refinement_area
    resolution_used = 0
    pass = 0
    while true
        pass_budget = budget === nothing ? nothing : budget - resolution_used
        pass_budget === 0 && break
        pass += 1
        verbose && println("Iterative refinement pass ", pass, ": refining until all incomplete triangles have area <= ", scaled_min_refinement_area(TC), ".")
        inserted = refine_one_pass!(TC; budget=pass_budget, verbose=verbose)
        resolution_used += inserted
        inserted == 0 && break
    end
    budget === nothing || resolution_used == budget || error("Unable to spend refinement budget; no eligible refinement points remain.")
    return resolution_used
end

# Refine until a local call budget is spent.
function refine_with_budget!(TC::TriangulationCache, budget::Integer; verbose=is_verbose(TC))
    budget >= 0 || error("Refinement budget must be nonnegative.")
    remaining = Int(budget)
    resolution_used = 0
    pass = 0
    while remaining > 0
        pass += 1
        verbose && println("Iterative refinement pass ", pass, ": remaining call budget before pass is ", remaining, ".")
        inserted = refine_one_pass!(TC; budget=remaining, verbose=verbose)
        resolution_used += inserted
        remaining -= inserted
        if inserted == 0
            println("Unable to spend refinement budget; no eligible refinement points remain.")
            remaining = 0
        end
    end
    return resolution_used
end

# Refine a TriangulationCache.
"""
    refine!(TC::TriangulationCache; min_refinement_area=nothing, budget=nothing, verbose=TC.verbose)

Refine a `TriangulationCache`.

By default, this performs one pass over the current incomplete triangles, adding
new points according to `TC.strategy` and skipping triangles whose area is at or
below `TC.min_refinement_area * window_area`.

If `budget` is supplied, refinement repeats until exactly that many new oracle
calls are spent. If no eligible refinement points remain before the budget is
spent, an error is thrown. Otherwise, one refinement pass uses `TC.oracle_budget`
as its cap.

If `min_refinement_area` is supplied, `TC.min_refinement_area` is updated and
refinement repeats until every incomplete triangle in the current window is at
or below the resulting area cutoff. Iterative minimum-area refinement requires
a positive, finite value.
"""
function refine!(TC::TriangulationCache; min_refinement_area=nothing, budget=nothing, verbose=is_verbose(TC))
    min_refinement_area === nothing || return refine_to_min_area!(TC, min_refinement_area; budget=budget, verbose=verbose)
    budget === nothing || return refine_with_budget!(TC, budget; verbose=verbose)
    return refine_one_pass!(TC; verbose=verbose)
end

# Refine until the oracle budget is spent.
function refine_until_budget_exhausted!(TC::TriangulationCache; verbose=is_verbose(TC))
    resolution_used = 0
    pass = 0
    while remaining_oracle_budget(TC) !== 0
        budget = remaining_oracle_budget(TC)
        budget === nothing && error("refine_until_budget_exhausted! requires a finite oracle_budget.")
        pass += 1
        verbose && println("Iterative refinement pass ", pass, ": remaining oracle budget before pass is ", budget, ".")
        inserted = refine_one_pass!(TC; verbose=verbose)
        resolution_used += inserted
        inserted == 0 && break
    end
    return resolution_used
end
