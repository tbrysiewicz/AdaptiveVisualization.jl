# GLMakie rendering and controls for TriangulationCache.

############################
# Color and Legend Helpers
############################

# Return the default colormap.
function default_continuous_colormap()
    return [
        GLMakie.RGBf(0.015, 0.015, 0.035),
        GLMakie.RGBf(0.200, 0.020, 0.250),
        GLMakie.RGBf(0.600, 0.040, 0.180),
        GLMakie.RGBf(0.900, 0.250, 0.070),
        GLMakie.RGBf(0.990, 0.720, 0.200),
        GLMakie.RGBf(1.000, 0.980, 0.760),
    ]
end

# Pad flat color ranges.
function finite_color_range(values::AbstractVector{<:Real})
    lo, hi = extrema(values)
    if lo == hi
        pad = max(abs(lo), 1.0) / 2
        return (lo - pad, hi + pad)
    end
    return (lo, hi)
end

# Build categorical colors.
function categorical_palette(n::Integer)
    base_colors = GLMakie.Makie.wong_colors()
    n <= length(base_colors) && return base_colors[1:n]

    extra_colors = map(GLMakie.Makie.to_color, [:purple, :brown, :pink, :gray, :olive, :cyan])
    colors = vcat(base_colors, extra_colors)
    while length(colors) < n
        append!(colors, colors)
    end
    return colors[1:n]
end

function value_label(value::Real)
    return value == round(value) ? string(Int(round(value))) : string(value)
end

function value_label(value::AbstractVector)
    return "[" * join(value_label.(value), ", ") * "]"
end

function value_label(value)
    return value isa Real ? value_label(value) : string(value)
end

# Add a value legend.
function add_value_legend!(fig, elements, labels, title)
    isempty(elements) && return nothing
    legend = GLMakie.Legend(fig[1, 2], elements, labels, title;
        tellwidth=true, tellheight=false,
        margin=(10, 10, 10, 10))
    GLMakie.colsize!(fig.layout, 2, GLMakie.Auto())
    return legend
end

############################
# Plot Value Helpers
############################

function triangle_plot_value(TC::TriangulationCache, triangle::Vector{Int64}; plot_log_transform=false)
    vertex_values = non_wildcard_values(function_values(TC)[triangle])
    isempty(vertex_values) && return nothing
    if all(is_real_value, vertex_values)
        value = numeric_mean_or_nothing(vertex_values)
        if plot_log_transform
            value <= -1 && return nothing
            value = log(value + 1)
        end
        return isfinite(value) ? value : nothing
    end

    all(==(first(vertex_values)), vertex_values) || error("Cannot plot a complete triangle with non-equal non-real vertex values: $(vertex_values).")
    return first(vertex_values)
end

function cached_plot_value(value; plot_log_transform=false)
    is_wildcard_value(value) && return nothing
    if value isa Real
        if plot_log_transform
            value <= -1 && return nothing
            value = log(value + 1)
        end
        return isfinite(value) ? value : nothing
    end
    return value
end

function categorical_unique_values_any(values)
    categories = unique(filter(value -> !isnothing(value) && !is_wildcard_value(value), values))
    return sort(categories; by=string)
end

function append_missing_values!(values::Vector{Any}, candidates)
    for value in candidates
        if !isnothing(value) && !is_wildcard_value(value) && !(value in values)
            push!(values, value)
        end
    end
    return values
end

function sort_numeric_values!(values::Vector{Any})
    isempty(values) && return values

    # Values of unrelated types (for example vector-valued invariants and
    # status strings) need not be mutually comparable. Keep type-groups in
    # discovery order and sort only within each homogeneous group.
    groups = Vector{Vector{Any}}()
    group_indices = Dict{DataType,Int}()
    for value in values
        group_index = get!(group_indices, typeof(value)) do
            push!(groups, Any[])
            length(groups)
        end
        push!(groups[group_index], value)
    end

    for group in groups
        try
            sorted_group = all(value -> value isa AbstractVector, group) ?
                sort(group; by=Tuple) : sort(group)
            copyto!(group, sorted_group)
        catch
            # Preserve discovery order when a homogeneous group is unsortable.
        end
    end

    copyto!(values, reduce(vcat, groups))
    return values
end

function stable_plot_value_order!(TC::TriangulationCache, visible_values; plot_log_transform=false)
    values = get!(TC.plot_value_order, Bool(plot_log_transform), Any[])
    cached_values = (cached_plot_value(value; plot_log_transform=plot_log_transform) for value in output_values(TC))
    append_missing_values!(values, cached_values)
    append_missing_values!(values, visible_values)
    sort_numeric_values!(values)
    return values
end

function vertex_plot_values(triangle_values)
    colors = Vector{Any}(undef, 3 * length(triangle_values))
    for (i, value) in pairs(triangle_values)
        # A category may itself be iterable (for example, a vector-valued
        # invariant). Treat it as one atomic value rather than broadcasting
        # its entries across the triangle's three vertices.
        colors[(3i - 2):(3i)] .= Ref(value)
    end
    return colors
end

function continuous_vertex_values(values, categories)
    if all(value -> value === nothing || value isa Real, values)
        return [value === nothing ? NaN : Float64(value) for value in values]
    end
    category_index = Dict(value => Float64(i) for (i, value) in pairs(categories))
    return [value === nothing ? NaN : category_index[value] for value in values]
end

function numeric_continuous_color_range(TC::TriangulationCache, visible_values; plot_log_transform=false)
    cached_values = [
        cached_plot_value(value; plot_log_transform=plot_log_transform)
        for value in output_values(TC)
    ]
    values = [
        Float64(value)
        for value in vcat(cached_values, visible_values)
        if value isa Real && isfinite(value)
    ]
    return isempty(values) ? (0.0, 1.0) : finite_color_range(values)
end

############################
# Mesh Assembly
############################

function selected_triangles(TC::TriangulationCache, plot_all_triangles::Bool)
    triangles = plot_all_triangles ?
        vcat(complete_triangles(TC), incomplete_triangles(TC)) : complete_triangles(TC)
    return [T for T in triangles if triangle_intersects_window(TC, T)]
end

function triangulation_mesh_data(TC::TriangulationCache, triangles::TriangleList, triangle_values)
    total_vertices = 3 * length(triangles)
    vertices = Matrix{Float64}(undef, total_vertices, 2)
    faces = Matrix{Int64}(undef, length(triangles), 3)

    for (triangle_index, triangle) in enumerate(triangles)
        vertex_offset = 3 * (triangle_index - 1)
        for local_index in 1:3
            point = get_point(triangulation(TC), triangle[local_index])
            row = vertex_offset + local_index
            vertices[row, 1] = point[1]
            vertices[row, 2] = point[2]
            faces[triangle_index, local_index] = row
        end
    end

    return vertices, faces
end

function triangle_edge_coordinates(vertices, faces)
    xs = Float64[]
    ys = Float64[]
    for face_index in axes(faces, 1)
        for vertex_index in (faces[face_index, 1], faces[face_index, 2], faces[face_index, 3], faces[face_index, 1])
            push!(xs, vertices[vertex_index, 1])
            push!(ys, vertices[vertex_index, 2])
        end
        push!(xs, NaN)
        push!(ys, NaN)
    end
    return xs, ys
end

function add_triangle_edges!(ax, decorations, vertices, faces; color, linewidth)
    edge_xs, edge_ys = triangle_edge_coordinates(vertices, faces)
    push!(decorations, GLMakie.lines!(ax, edge_xs, edge_ys; color=color, linewidth=linewidth))
    return decorations
end

############################
# Drawing Lifecycle
############################

function set_axis_window!(ax, TC::TriangulationCache)
    GLMakie.xlims!(ax, TC.xlims[1], TC.xlims[2])
    GLMakie.ylims!(ax, TC.ylims[1], TC.ylims[2])
    return ax
end

function draw_triangulation!(fig, ax, TC::TriangulationCache; kwargs...)
    plot_log_transform = get(kwargs, :plot_log_transform, false)
    plot_all_triangles = get(kwargs, :plot_all_triangles, !is_discrete(TC))
    colormap = get(kwargs, :colormap, default_continuous_colormap())
    show_legend = get(kwargs, :show_legend, true)
    legend_max_values = get(kwargs, :legend_max_values, TRIANGULATION_CACHE_DEFAULT_LEGEND_LIMIT)
    discrete_legend = get(kwargs, :discrete_legend, nothing)
    legend_title = get(kwargs, :legend_title, "value")
    plot_triangle_edges = get(kwargs, :edges, get(kwargs, :plot_triangle_edges, false))
    triangle_edge_color = get(kwargs, :triangle_edge_color, GLMakie.RGBAf(0, 0, 0, 0.35))
    triangle_edge_linewidth = get(kwargs, :triangle_edge_linewidth, 0.5)
    decorations = Any[]

    triangles = selected_triangles(TC, plot_all_triangles)
    triangle_values = [triangle_plot_value(TC, T; plot_log_transform=plot_log_transform) for T in triangles]
    vertices, faces = triangulation_mesh_data(TC, triangles, triangle_values)
    vertex_values = vertex_plot_values(triangle_values)

    if size(faces, 1) == 0
        return (plot=nothing, decorations=decorations)
    end

    visible_categories = categorical_unique_values_any(vertex_values)
    if isempty(visible_categories)
        plt = GLMakie.mesh!(ax, vertices, faces; color=:black, shading=false)
        plot_triangle_edges && add_triangle_edges!(ax, decorations, vertices, faces; color=triangle_edge_color, linewidth=triangle_edge_linewidth)
        return (plot=plt, decorations=decorations)
    end

    categories = stable_plot_value_order!(TC, vertex_values; plot_log_transform=plot_log_transform)
    use_discrete_legend = discrete_legend === nothing ? 0 < length(categories) <= legend_max_values : discrete_legend
    if use_discrete_legend
        category_colors = categorical_palette(length(categories))
        color_map = Dict(value => color for (value, color) in zip(categories, category_colors))
        mesh_colors = [value === nothing ? GLMakie.RGBAf(0, 0, 0, 1) : color_map[value] for value in vertex_values]
        plt = GLMakie.mesh!(ax, vertices, faces; color=mesh_colors, shading=false)
        plot_triangle_edges && add_triangle_edges!(ax, decorations, vertices, faces; color=triangle_edge_color, linewidth=triangle_edge_linewidth)
        elements = [GLMakie.PolyElement(color=color, strokecolor=color) for color in category_colors]
        labels = value_label.(categories)
        if show_legend
            legend = add_value_legend!(fig, elements, labels, legend_title)
            legend === nothing || push!(decorations, legend)
        end
        return (plot=plt, decorations=decorations)
    end

    continuous_values = continuous_vertex_values(vertex_values, categories)
    colorrange = all(value -> value === nothing || value isa Real, vertex_values) ?
        numeric_continuous_color_range(TC, vertex_values; plot_log_transform=plot_log_transform) :
        finite_color_range(collect(1.0:length(categories)))
    plt = GLMakie.mesh!(ax, vertices, faces; color=continuous_values, colormap=colormap, colorrange=colorrange, nan_color=:black, shading=false)
    plot_triangle_edges && add_triangle_edges!(ax, decorations, vertices, faces; color=triangle_edge_color, linewidth=triangle_edge_linewidth)
    show_legend && push!(decorations, GLMakie.Colorbar(fig[1, 2], plt))
    return (plot=plt, decorations=decorations)
end

function delete_drawn_triangulation!(ax, drawn)
    drawn.plot === nothing || GLMakie.delete!(ax, drawn.plot)
    for decoration in drawn.decorations
        try
            GLMakie.delete!(ax, decoration)
        catch
            try
                GLMakie.delete!(decoration)
            catch
            end
        end
    end
    return nothing
end

function redraw_triangulation!(fig, ax, TC::TriangulationCache, drawn_ref; kwargs...)
    set_axis_window!(ax, TC)
    delete_drawn_triangulation!(ax, drawn_ref[])
    drawn_ref[] = draw_triangulation!(fig, ax, TC; kwargs...)
    return drawn_ref[]
end

############################
# Interactive Controls
############################

function min_refinement_area_from_slider_exponent(exponent::Real)
    return 10.0 ^ (-Int(round(exponent)))
end

function default_min_refinement_area_slider_range()
    return 2:6
end

function default_min_refinement_area_slider_start(TC::TriangulationCache, slider_range)
    if !isfinite(TC.min_refinement_area) || TC.min_refinement_area <= 0
        return last(slider_range)
    end
    default_start = Int(round(-log10(TC.min_refinement_area)))
    return slider_range[argmin(abs.(slider_range .- default_start))]
end

function min_refinement_area_slider_label(exponent::Real)
    return "Min area: window * 1e-" * string(Int(round(exponent)))
end

function navigate_and_refine!(
        fig,
        ax,
        TC::TriangulationCache,
        drawn_ref;
        translate=(0.0, 0.0),
        zoom_factor=nothing,
        seed_and_refine=true,
        navigation_initial_resolution=250,
        navigation_refinement_budget=1000,
        verbose=is_verbose(TC),
        kwargs...)

    change_window_and_refine!(
        TC;
        translate=translate,
        zoom_factor=zoom_factor,
        seed_and_refine=seed_and_refine,
        navigation_initial_resolution=navigation_initial_resolution,
        navigation_refinement_budget=navigation_refinement_budget,
        verbose=verbose,
    )
    redraw_triangulation!(fig, ax, TC, drawn_ref; kwargs...)
    return TC
end

function add_refine_button!(
        fig,
        ax,
        TC::TriangulationCache,
        drawn_ref;
        button_refinement_passes=1,
        navigation_step=0.20,
        zoom_step=0.20,
        min_refinement_area_controls=true,
        min_refinement_area_slider_range=default_min_refinement_area_slider_range(),
        min_refinement_area_slider_start=nothing,
        navigation_initial_resolution=250,
        navigation_refinement_budget=1000,
        plot_triangle_edges=false,
        triangle_edge_color=GLMakie.RGBAf(0, 0, 0, 0.35),
        triangle_edge_linewidth=0.5,
        verbose=is_verbose(TC),
        kwargs...)

    controls = fig[2, :] = GLMakie.GridLayout(tellwidth=false, tellheight=true)
    refine_button = GLMakie.Button(controls[1, 1]; label="Refine", tellwidth=false, width=120, height=30)
    left_button = GLMakie.Button(controls[1, 2]; label="←", tellwidth=false, width=42, height=30)
    right_button = GLMakie.Button(controls[1, 3]; label="→", tellwidth=false, width=42, height=30)
    up_button = GLMakie.Button(controls[1, 4]; label="↑", tellwidth=false, width=42, height=30)
    down_button = GLMakie.Button(controls[1, 5]; label="↓", tellwidth=false, width=42, height=30)
    zoom_in_button = GLMakie.Button(controls[1, 6]; label="Zoom +", tellwidth=false, width=84, height=30)
    zoom_out_button = GLMakie.Button(controls[1, 7]; label="Zoom -", tellwidth=false, width=84, height=30)
    edge_button = GLMakie.Button(controls[1, 8]; label="Edges", tellwidth=false, width=84, height=30)
    edge_visible = Ref(Bool(get(kwargs, :edges, plot_triangle_edges)))
    current_plot_kwargs() = merge((; kwargs...),
        (; edges=edge_visible[], triangle_edge_color, triangle_edge_linewidth))
    GLMakie.rowgap!(fig.layout, 8)
    GLMakie.rowsize!(fig.layout, 2, GLMakie.Fixed(min_refinement_area_controls ? 92 : 52))
    GLMakie.rowsize!(controls, 1, GLMakie.Fixed(40))
    GLMakie.rowgap!(controls, 4)

    redraw() = redraw_triangulation!(
        fig,
        ax,
        TC,
        drawn_ref;
        current_plot_kwargs()...)

    if min_refinement_area_controls
        slider_range = collect(min_refinement_area_slider_range)
        isempty(slider_range) && error("min_refinement_area_slider_range must contain at least one value.")
        slider_start = min_refinement_area_slider_start === nothing ? default_min_refinement_area_slider_start(TC, slider_range) : min_refinement_area_slider_start
        min_refinement_area_slider = GLMakie.Slider(controls[2, 2:4]; range=slider_range, startvalue=slider_start, tellwidth=true, width=260)
        min_refinement_area_label = GLMakie.Label(
            controls[2, 1],
            GLMakie.lift(min_refinement_area_slider_label, min_refinement_area_slider.value);
            tellwidth=false,
            width=160,
        )
        fully_refine_button = GLMakie.Button(controls[2, 5]; label="Fully Refine", tellwidth=false, width=120, height=30)
        GLMakie.rowsize!(controls, 2, GLMakie.Fixed(40))

        GLMakie.on(fully_refine_button.clicks) do _
            refine!(TC; min_refinement_area=min_refinement_area_from_slider_exponent(min_refinement_area_slider.value[]), verbose=verbose)
            redraw()
        end
    end

    GLMakie.on(refine_button.clicks) do _
        for _ in 1:button_refinement_passes
            refine!(TC)
        end
        redraw()
    end

    navigate(; translate=(0.0, 0.0), zoom_factor=nothing, seed_and_refine=true) = navigate_and_refine!(
        fig,
        ax,
        TC,
        drawn_ref;
        translate=translate,
        zoom_factor=zoom_factor,
        seed_and_refine=seed_and_refine,
        navigation_initial_resolution=navigation_initial_resolution,
        navigation_refinement_budget=navigation_refinement_budget,
        verbose=verbose,
        current_plot_kwargs()...)

    GLMakie.on(edge_button.clicks) do _
        edge_visible[] = !edge_visible[]
        redraw()
    end

    GLMakie.on(left_button.clicks) do _
        navigate(translate=(-navigation_step, 0.0))
    end
    GLMakie.on(right_button.clicks) do _
        navigate(translate=(navigation_step, 0.0))
    end
    GLMakie.on(up_button.clicks) do _
        navigate(translate=(0.0, navigation_step))
    end
    GLMakie.on(down_button.clicks) do _
        navigate(translate=(0.0, -navigation_step))
    end
    GLMakie.on(zoom_in_button.clicks) do _
        navigate(zoom_factor=1 - zoom_step, seed_and_refine=false)
    end
    GLMakie.on(zoom_out_button.clicks) do _
        navigate(zoom_factor=1 + zoom_step)
    end

    return refine_button
end

############################
# Figure API
############################

"""
    visualize(TC::TriangulationCache; kwargs...) -> GLMakie.Figure

Render a `TriangulationCache` using GLMakie.

Useful keyword arguments:
- `buttons`: add interactive refinement controls, default `true`.
- `button_refinement_passes`: number of refinement passes per button click.
- `navigation_step`: arrow-button pan amount as a fraction of window size.
- `zoom_step`: zoom amount as a fraction of window size.
- `navigation_refinement_budget`: maximum new refinement samples after pan/zoom,
  in addition to samples used to seed the new window.
- `navigation_initial_resolution`: coarse mesh size seeded after pan/zoom.
- `min_refinement_area_controls`: add a slider that selects the normalized
  minimum used by the Fully Refine button.
- `min_refinement_area_slider_range`: integer exponents for a normalized
  minimum area of `1e-X`, default `2:6`.
- `figure_size`: Makie figure size, default `(1300, 900)` with buttons enabled
  or `(900, 900)` when `buttons=false`.
- `figure_padding`: outer margins `(left, right, bottom, top)`, default
  `(20, 20, 20, 40)`. A single number applies to all sides.
- `xlabel`: axis x label, default `"x"`.
- `ylabel`: axis y label, default `"y"`.
- `xlabelsize`: axis x label font size, default Makie axis label size.
- `ylabelsize`: axis y label font size, default Makie axis label size.
- `xticklabelsize`: axis x tick label font size, default Makie tick label size.
- `yticklabelsize`: axis y tick label font size, default Makie tick label size.
- `title`: axis title, default `""`.
- `titlesize`: axis title font size, default Makie axis title size.
- `plot_all_triangles`: include incomplete triangles in the colored mesh.
  Defaults to `false` for discrete caches and `true` for continuous caches.
- `edges`: overlay thin triangle edges, default `false`.
- `plot_triangle_edges`: deprecated alias for `edges`.
- `triangle_edge_color`: edge overlay color.
- `triangle_edge_linewidth`: edge overlay line width.
- `show_legend`: show legends and colorbars, default `true`.
- `legend_max_values`: categorical legend threshold, default `20`.
- `discrete_legend`: `true` selects a categorical legend, `false` selects a
  continuous colorbar, and the default `nothing` chooses using `legend_max_values`.
  `show_legend=false` hides either; no legend is drawn if no non-wildcard values
  are visible.
"""
function visualize(TC::TriangulationCache; kwargs...)::GLMakie.Figure
    buttons = get(kwargs, :buttons, true)
    button_refinement_passes = get(kwargs, :button_refinement_passes, 1)

    figure_size = get(kwargs, :figure_size, buttons ? (1300, 900) : (900, 900))
    xlabel = get(kwargs, :xlabel, "x")
    ylabel = get(kwargs, :ylabel, "y")
    title = get(kwargs, :title, "")
    axis_kwargs = Dict{Symbol,Any}(
        :xlabel => xlabel,
        :ylabel => ylabel,
        :title => title,
    )
    for key in (:xlabelsize, :ylabelsize, :xticklabelsize, :yticklabelsize, :titlesize)
        haskey(kwargs, key) && (axis_kwargs[key] = kwargs[key])
    end
    figure_padding = get(kwargs, :figure_padding, (20, 20, 20, 40))
    fig = GLMakie.Figure(size=figure_size, figure_padding=figure_padding)
    ax = GLMakie.Axis(fig[1, 1]; axis_kwargs..., aspect=GLMakie.DataAspect(), backgroundcolor=:black)
    GLMakie.colsize!(fig.layout, 1, GLMakie.Relative(0.82))
    # Use the space remaining after the fixed-height controls and margins.
    GLMakie.rowsize!(fig.layout, 1, GLMakie.Auto())
    set_axis_window!(ax, TC)

    drawn_ref = Ref{Any}(draw_triangulation!(fig, ax, TC; kwargs...))

    buttons && add_refine_button!(fig, ax, TC, drawn_ref; button_refinement_passes=button_refinement_passes, kwargs...)

    return fig
end

const TRIANGULATION_CACHE_VISUALIZE_KEYWORDS = Set([
    :xlims,
    :ylims,
    :strategy,
    :is_complete,
])

"""
    visualize(function_oracle::Function; kwargs...) -> (TriangulationCache, GLMakie.Figure)

Construct, refine, display, and return a `TriangulationCache` and its Makie
figure. `total_resolution` is the maximum number of sampled input points
(default `1000`), counting initialization and refinement, unless
`min_refinement_area` is supplied. In that case, refinement continues until no
incomplete triangle in the current window exceeds
`min_refinement_area * window_area`. If `initial_resolution` is supplied, it
controls the initialization mesh size; otherwise initialization uses one
quarter of `total_resolution`. An explicit `initial_resolution` larger than the
default total raises the effective total unless `total_resolution` was explicitly
supplied. The remaining budget uses the number of points actually initialized.

The `strategy` keyword selects refinement points for incomplete triangles:
`:sierpinski` (default) samples the three edge midpoints, `:barycenter` samples
the centroid, and `:random` samples a random interior point. The selected strategy
is stored in the returned cache and used for subsequent refinement.

For ordinary functions, `input_dimension=n` specifies the number of coordinates
in each input point and selects a random two-dimensional plane through the
origin. Its default is `nothing`; supply `near` or `plane_points=[p,q,r]` to
infer the input dimension from those points instead. For two-dimensional inputs,
the original axis directions are retained. The plot always has two axes.
`zoomer` scales the plane and `rng` selects its random directions. The original coordinates are retained in
`TC.parameter_slice` for `retrieve_witnesses(TC)`. Set `batched=false` for a
single-point function whose batch call could look valid by accident; use
`batched=true` to require a batch oracle. Omission keeps automatic detection.
Sampling budgets count individual points, even when the function evaluates
several points in a single batched call.

Figure options such as `discrete_legend`, `figure_padding`, and `buttons` are
forwarded to `visualize(TC; kwargs...)`; see that method’s keyword descriptions.
"""
function visualize(function_oracle::Function; total_resolution=nothing, initial_resolution=nothing, resolution=nothing, min_refinement_area=nothing, verbose=false, batched=nothing, input_dimension=nothing, near=nothing, plane_points=nothing, zoomer=1.0, rng::AbstractRNG=default_rng(), kwargs...)
    resolution === nothing || total_resolution === nothing ||
        error("Use `total_resolution`, not both `resolution` and `total_resolution`.")
    explicit_total = total_resolution !== nothing || resolution !== nothing
    resolved_total = Int(total_resolution === nothing ?
        (resolution === nothing ? TRIANGULATION_CACHE_DEFAULT_TOTAL_RESOLUTION : resolution) : total_resolution)
    resolved_total >= 4 || error("total_resolution must be at least 4.")
    resolved_initial = initial_resolution === nothing ? max(4, round(Int, 0.25 * resolved_total)) : Int(initial_resolution)
    resolved_initial >= 4 || error("initial_resolution must be at least 4.")
    if resolved_initial > resolved_total
        explicit_total && error("initial_resolution must be less than or equal to total_resolution.")
        resolved_total = resolved_initial
    end

    cache_kwargs = Dict{Symbol,Any}()
    plot_kwargs = Dict{Symbol,Any}()
    for (key, value) in kwargs
        if key in TRIANGULATION_CACHE_VISUALIZE_KEYWORDS
            cache_kwargs[key] = value
        else
            plot_kwargs[key] = value
        end
    end

    min_refinement_area === nothing || (cache_kwargs[:min_refinement_area] = Float64(min_refinement_area))
    TC = TriangulationCache(function_oracle; resolution=resolved_initial, verbose=verbose,
        batched, input_dimension, near, plane_points, zoomer, rng, cache_kwargs...)
    refinement_budget = resolved_total - TC.total_oracle_calls
    if min_refinement_area !== nothing
        refine!(TC; min_refinement_area=Float64(min_refinement_area), verbose=verbose)
        TC.oracle_budget = nothing
    else
        TC.oracle_budget = refinement_budget
        refine_until_budget_exhausted!(TC; verbose=verbose)
        TC.oracle_budget = nothing
    end
    fig = visualize(TC; plot_kwargs...)
    display(fig)
    return TC, fig
end

"""
    save(fig::GLMakie.Figure, filename::String; file_extension="png", dpi=300)

Save a GLMakie figure under `OutputFiles/` unless `filename` already starts
with that directory. The path is relative to the current working directory;
missing parent directories are created. If the filename does not end with
`file_extension`, append `"." * file_extension` (an existing different extension
is not replaced). Supply the extension without a leading dot.

`dpi` controls raster scaling through Makie's `px_per_unit=dpi/150`; the default
`dpi=300` uses two pixels per figure unit. Returns the final filename.
"""
function save(fig::GLMakie.Figure, filename::String; file_extension = "png", dpi = 300)
    if !startswith(filename, "OutputFiles/")
        filename = "OutputFiles/" * filename
    end
    if !endswith(filename, file_extension)
        filename *= "." * file_extension
    end
    mkpath(dirname(filename))
    GLMakie.save(filename, fig; px_per_unit=dpi/150)
    println("Plot saved to $filename")
    return filename
end
