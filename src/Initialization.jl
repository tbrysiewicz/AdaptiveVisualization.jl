# Construction and oracle initialization.

function try_evaluate_triangulation_oracle(function_oracle::Function, points; batched=nothing)
    errors = Tuple{Symbol,Any}[]
    if batched !== false
        try
            oracle = checked_batch_oracle(function_oracle)
            return oracle, oracle(points)
        catch err
            push!(errors, (:batched_points, err))
            batched === true && return nothing, errors
        end
    end
    try
        oracle = pts -> map(p -> function_oracle(p), pts)
        return oracle, oracle(points)
    catch err
        push!(errors, (:single_point_vector, err))
        try
            oracle = single_to_batched_oracle(function_oracle)
            return oracle, oracle(points)
        catch err
            push!(errors, (:single_point_coordinates, err))
            return nothing, errors
        end
    end
end
function oracle_error_message(errors)
    message = IOBuffer()
    println(message, "Function oracle must either accept a vector of points and return a vector of values, accept one point as f(point), or accept one point as f(x, y).")
    println(message, "Tried these call forms:")
    for (form, err) in errors
        print(message, "  ", form, ": ")
        showerror(message, err)
        println(message)
    end
    return String(take!(message))
end

# Adapt the original oracle once, then apply the same affine map to every batch.
function sliced_function_oracle(function_oracle::Function; input_dimension=nothing,
        near=nothing, plane_points=nothing, zoomer=1.0, rng::AbstractRNG=default_rng(),
        batched=nothing)
    if function_oracle isa SlicedOracle &&
            (input_dimension !== nothing || near !== nothing || plane_points !== nothing || zoomer != 1.0)
        throw(ArgumentError("Cannot apply another slice to an existing SlicedOracle."))
    end
    if input_dimension === nothing && near === nothing && plane_points === nothing && zoomer == 1.0
        return function_oracle
    end
    input_dimension === nothing || input_dimension isa Integer ||
        throw(ArgumentError("input_dimension must be an integer."))
    if plane_points !== nothing
        (plane_points isa AbstractVector || plane_points isa Tuple) && length(plane_points) == 3 ||
            throw(ArgumentError("plane_points must contain exactly three points [p,q,r]."))
    end
    n = input_dimension === nothing ?
        (plane_points === nothing ? (near === nothing ? 2 : length(near)) : length(first(plane_points))) :
        input_dimension
    slice = affine_parameter_slice(n; near, plane_points, zoomer, rng)
    adapted = Ref{Any}(nothing)
    function evaluate(points)
        mapped = slice.(points)
        if adapted[] === nothing
            oracle, values = try_evaluate_triangulation_oracle(function_oracle, mapped; batched)
            oracle === nothing && error(oracle_error_message(values))
            adapted[] = oracle
            return values
        end
        return adapted[](mapped)
    end
    return SlicedOracle(evaluate, slice)
end

function TriangulationCache(
        function_oracle::Function;
        xlims = [-1, 1],
        ylims = [-1, 1],
        resolution = TRIANGULATION_CACHE_DEFAULT_TOTAL_RESOLUTION,
        strategy::Symbol = :sierpinski,
        min_refinement_area = TRIANGULATION_CACHE_DEFAULT_MIN_REFINEMENT_AREA,
        is_complete = nothing,
        verbose::Bool = false,
        batched=nothing,
        input_dimension=nothing,
        near=nothing,
        plane_points=nothing,
        zoomer=1.0,
        rng::AbstractRNG=default_rng(),
        kwargs...)

    batched === nothing || batched isa Bool ||
        throw(ArgumentError("batched must be true, false, or nothing."))
    already_sliced = function_oracle isa SlicedOracle
    function_oracle = sliced_function_oracle(function_oracle;
        input_dimension, near, plane_points, zoomer, rng, batched)
    already_sliced && batched === false &&
        throw(ArgumentError("An existing SlicedOracle requires batched evaluation."))
    isempty(kwargs) || error("Unsupported TriangulationCache keyword(s): $(join(keys(kwargs), ", ")).")
    strategy in TRIANGULATION_CACHE_STRATEGIES || error("Invalid strategy $strategy. Use one of $(TRIANGULATION_CACHE_STRATEGIES).")

    resolution = Int(resolution)
    resolution >= 4 || error("resolution must be at least 4 so the initial triangulation can be two-dimensional.")

    xlimits = Float64.(collect(xlims))
    ylimits = Float64.(collect(ylims))
    length(xlimits) == 2 || error("xlims must have two entries.")
    length(ylimits) == 2 || error("ylims must have two entries.")
    resolved_min_refinement_area = Float64(min_refinement_area)
    isfinite(resolved_min_refinement_area) || error("min_refinement_area must be finite.")
    resolved_min_refinement_area >= 0 || error("min_refinement_area must be nonnegative.")

    parameters = triangulation_initial_points(xlimits, ylimits, resolution)
    evaluation = try_evaluate_triangulation_oracle(function_oracle, parameters;
        batched=function_oracle isa SlicedOracle ? true : batched)
    first(evaluation) === nothing && error(oracle_error_message(last(evaluation)))
    batched_oracle, values = evaluation

    function_values = FunctionValues(values)
    verbose && println("Building initial Delaunay triangulation globally from ", length(parameters), " sampled points.")
    tri = triangulate(Tuple.(parameters))
    tc_is_complete = is_complete === nothing ? default_is_complete(function_values) : is_complete

    TC = TriangulationCache(
        batched_oracle,
        oracle_parameter_slice(function_oracle),
        tri,
        function_values,
        point_index_dict(parameters),
        Set{TriangleKey}(),
        tc_is_complete,
        strategy,
        length(parameters),
        nothing,
        resolved_min_refinement_area,
        xlimits,
        ylimits,
        [(xlimits[1], xlimits[2], ylimits[1], ylimits[2])],
        Dict{Bool,Vector{Any}}(),
        verbose,
    )
    recompute_incomplete_triangles!(TC)
    return TC
end

function initial_parameter_distribution(; kwargs...)
    xlims = get(kwargs, :xlims, [-1, 1])
    ylims = get(kwargs, :ylims, [-1, 1])
    resolution = get(kwargs, :resolution, 1000)
    xlength = abs(xlims[2] - xlims[1])
    ylength = abs(ylims[2] - ylims[1])

    nx = clamp(floor(Int, sqrt(resolution * xlength / ylength)), 2, div(resolution, 2))
    ny = clamp(floor(Int, sqrt(resolution * ylength / xlength)), 2, div(resolution, nx))

    return range(xlims[1], xlims[2], length=nx), range(ylims[1], ylims[2], length=ny)
end

function triangulation_initial_points(xlims::Vector{Float64}, ylims::Vector{Float64}, resolution::Integer)
    x_values, y_values = initial_parameter_distribution(; xlims=xlims, ylims=ylims, resolution=resolution)
    return [[Float64(x), Float64(y)] for x in x_values for y in y_values]
end
