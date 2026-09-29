# Mathematical window operations for TriangulationCache.

window_width(TC::TriangulationCache) = TC.xlims[2] - TC.xlims[1]
window_height(TC::TriangulationCache) = TC.ylims[2] - TC.ylims[1]
bounding_box_area(TC::TriangulationCache) = window_width(TC) * window_height(TC)

function translate_window!(TC::TriangulationCache, dx_fraction::Real, dy_fraction::Real)
    dx = Float64(dx_fraction) * window_width(TC)
    dy = Float64(dy_fraction) * window_height(TC)
    TC.xlims .+= dx
    TC.ylims .+= dy
    return TC
end

function zoom_window!(TC::TriangulationCache, factor::Real)
    factor = Float64(factor)
    factor > 0 || error("Zoom factor must be positive.")
    xcenter = sum(TC.xlims) / 2
    ycenter = sum(TC.ylims) / 2
    half_width = window_width(TC) * factor / 2
    half_height = window_height(TC) * factor / 2
    TC.xlims .= [xcenter - half_width, xcenter + half_width]
    TC.ylims .= [ycenter - half_height, ycenter + half_height]
    return TC
end

current_window(TC::TriangulationCache) = (TC.xlims[1], TC.xlims[2], TC.ylims[1], TC.ylims[2])

window_contains(outer::NTuple{4,Float64}, inner::NTuple{4,Float64}) = outer[1] <= inner[1] && outer[2] >= inner[2] && outer[3] <= inner[3] && outer[4] >= inner[4]
window_is_covered(TC::TriangulationCache, window=current_window(TC)) = any(covered -> window_contains(covered, window), TC.covered_windows)

function seed_window_mesh!(TC::TriangulationCache; resolution::Integer, verbose=is_verbose(TC))
    points = triangulation_initial_points(TC.xlims, TC.ylims, resolution)
    new_points = Vector{Vector{Float64}}()
    queued_keys = Set{PointKey}()
    for p in points
        key = point_key(p)
        if !haskey(point_indices(TC), key) && !(key in queued_keys)
            push!(new_points, p)
            push!(queued_keys, key)
        end
    end
    inserted = evaluate_and_insert_points!(TC, new_points; verbose=verbose)
    TC.total_oracle_calls += inserted
    return inserted
end

function change_window_and_refine!(
        TC::TriangulationCache;
        translate=(0.0, 0.0),
        zoom_factor=nothing,
        seed_and_refine=true,
        navigation_initial_resolution=250,
        navigation_refinement_budget=1000,
        verbose=is_verbose(TC))

    zoom_factor === nothing || zoom_window!(TC, zoom_factor)
    translate_window!(TC, translate[1], translate[2])
    needs_new_sampling = seed_and_refine && !window_is_covered(TC)
    if needs_new_sampling
        seed_window_mesh!(TC; resolution=navigation_initial_resolution, verbose=verbose)
        refine_with_budget!(TC, navigation_refinement_budget; verbose=verbose)
        #record the window
        push!(TC.covered_windows, current_window(TC))
    end
    return TC
end
