module HomotopyContinuationExt

import AdaptiveVisualization as AV
import AdaptiveVisualization: kuramoto_model
import HomotopyContinuation as HC
using HomotopyContinuation: @var
using Random: AbstractRNG, default_rng, randn, MersenneTwister

# One validator for real slice coordinates and complex continuation data.
function checked_vector(::Type{T}, values, n, name) where {T}
    (values isa AbstractVector || values isa Tuple) && length(values) == n ||
        throw(ArgumentError("$name must be a vector of length $n."))
    element_type = T <: Real ? Real : Number
    all(x -> x isa element_type && isfinite(x), values) ||
        throw(ArgumentError("$name must contain finite $element_type values."))
    result = T.(collect(values))
    all(isfinite, result) || throw(ArgumentError("$name must be representable as finite $T values."))
    return result
end

function parameter_slice(F::HC.System; near=nothing, plane_points=nothing,
        zoomer=1.0, rng::AbstractRNG=default_rng())
    return AV.affine_parameter_slice(HC.nparameters(F); near, plane_points, zoomer, rng)
end

parameter_values(slice, points) =
    [slice(checked_vector(Float64, point, 2, "plot point")) for point in points]

function start_fibre(F, start_parameters, start_solutions, rng, solver_options,
        monodromy_options, verbose)
    P = start_parameters === nothing ? nothing :
        checked_vector(ComplexF64, start_parameters, HC.nparameters(F), "start_parameters")
    if start_solutions !== nothing
        P === nothing && throw(ArgumentError("start_solutions requires start_parameters."))
        isempty(start_solutions) && throw(ArgumentError("start_solutions must not be empty."))
        S = [checked_vector(ComplexF64, s, length(HC.variables(F)), "start solution") for s in start_solutions]
        return (; parameters=P, solutions=S)
    end

    options = merge((; show_progress=verbose, seed=rand(rng, UInt32)),
                    monodromy_options, (; catch_interrupt=false))
    M = HC.monodromy_solve(F; options...)
    S = Vector{ComplexF64}.(HC.solutions(M))
    isempty(S) && error("HC found no start solutions. Supply start_parameters and start_solutions explicitly.")
    source = ComplexF64.(HC.parameters(M))
    if P !== nothing
        result = HC.solve(F, S; solver_options..., start_parameters=source, target_parameters=P)
        moved = HC.solutions(result)
        HC.nfailed(result) == 0 && length(moved) == length(S) ||
            error("Could not transport all start solutions to start_parameters. Choose a generic start parameter.")
        return (; parameters=P, solutions=Vector{ComplexF64}.(moved))
    end
    return (; parameters=source, solutions=S)
end

function solution_counter(F::HC.System, mode::Symbol;
        near=nothing, plane_points=nothing, zoomer=1.0, rng::AbstractRNG=default_rng(),
        start_parameters=nothing, start_solutions=nothing, max_retries=nothing,
        real_tol=1e-8, positivity_tol=1e-7, imaginary_zero_atol=1e-10,
        wildcard_parity_mismatch=true, solver_options=(;), monodromy_options=(;),
        certification_options=(;), verbose=false)
    slice = parameter_slice(F; near, plane_points, zoomer, rng)
    retries = max_retries === nothing ? (mode === :certify_real ? 5 : 2) : max_retries
    retries isa Integer && retries >= 0 || throw(ArgumentError("max_retries must be a nonnegative integer."))
    for (name, tolerance) in ((:real_tol, real_tol), (:positivity_tol, positivity_tol),
                              (:imaginary_zero_atol, imaginary_zero_atol))
        tolerance isa Real && isfinite(tolerance) && tolerance >= 0 && isfinite(Float64(tolerance)) ||
            throw(ArgumentError("$name must be finite and nonnegative."))
    end
    real_tol, positivity_tol, imaginary_zero_atol = Float64.((real_tol, positivity_tol, imaginary_zero_atol))
    # Keep HC tracking options configurable without allowing the batch shape to change.
    reserved = (:start_parameters, :target_parameters, :transform_result,
                :transform_parameters, :flatten, :iterator_only, :target_subspaces,
                :bitmask, :stop_early_cb, :catch_interrupt)
    any(key -> haskey(solver_options, key), reserved) &&
        throw(ArgumentError("solver_options cannot override parameter batching or result handling."))
    options = merge((; show_progress=verbose), solver_options, (; catch_interrupt=false))
    certification = merge((; show_progress=verbose, max_precision=16384), certification_options)
    base = start_fibre(F, start_parameters, start_solutions, rng, options, monodromy_options, verbose)
    degree = length(base.solutions)
    fibres = Union{Nothing,typeof(base)}[base]
    verbose && @info "HC start solutions prepared" count=degree

    function sample_value(result, parameters)
        HC.nfailed(result) == 0 || return :wildcard
        solutions = HC.solutions(result)
        length(solutions) == degree || return :wildcard

        if mode === :dietmaier
            norms = filter(>(imaginary_zero_atol), [sum(abs, imag.(s)) for s in solutions])
            return isempty(norms) ? 0.0 : minimum(norms)
        elseif mode === :certify_real
            certificates = HC.distinct_certificates(HC.certify(F, solutions, parameters; certification...))
            total = count(HC.is_certified, certificates)
            nreal = count(c -> HC.is_certified(c) && HC.is_real(c), certificates)
            # Nonreal requires an imaginary interval excluding zero, not just a
            # failed real certificate. HC.is_complex performs this interval check.
            nonreal = count(c -> HC.is_certified(c) && HC.is_complex(c), certificates)
            total == degree && nreal + nonreal == total || return :wildcard
        else
            real_solutions = filter(s -> all(z -> abs(imag(z)) <= real_tol, s), solutions)
            nreal = length(real_solutions)
        end
        wildcard_parity_mismatch && isodd(nreal) != isodd(degree) && return :wildcard
        # Parity applies to the total real count, never to the positive subset.
        return mode === :positive ?
            count(s -> all(z -> real(z) > positivity_tol, s), real_solutions) : nreal
    end

    # The closure keeps the slice and lazily prepared retry fibres across batches.
    function evaluate(points)
        isempty(points) && return Union{Int,Float64,Symbol}[]
        targets = parameter_values(slice, points)
        values = Union{Int,Float64,Symbol}[:wildcard for _ in targets]
        pending = collect(eachindex(targets))
        for attempt in 0:retries
            isempty(pending) && break
            if attempt + 1 > length(fibres)
                P = randn(rng, ComplexF64, length(base.parameters))
                result = HC.solve(F, base.solutions; options...,
                                  start_parameters=base.parameters, target_parameters=P)
                S = HC.solutions(result)
                fibre = HC.nfailed(result) == 0 && length(S) == degree ?
                    (; parameters=P, solutions=Vector{ComplexF64}.(S)) : nothing
                push!(fibres, fibre)
            end
            fibre = fibres[attempt + 1]
            fibre === nothing && continue
            verbose && attempt > 0 && @info "Retrying unresolved HC samples" attempt count=length(pending)
            results = HC.solve(F, fibre.solutions; options...,
                start_parameters=fibre.parameters, target_parameters=targets[pending],
                transform_result=(result, _) -> result, flatten=false)
            for (i, result) in zip(pending, results)
                values[i] = sample_value(result, targets[i])
            end
            filter!(i -> values[i] === :wildcard, pending)
        end
        return values
    end
    return AV.SlicedOracle(evaluate, slice)
end

AV.real_solution_function(F::HC.System; near=nothing, plane_points=nothing, kwargs...) =
    solution_counter(F, :real; near, plane_points, kwargs...)

AV.positive_solution_function(F::HC.System; near=nothing, plane_points=nothing, kwargs...) =
    solution_counter(F, :positive; near, plane_points, kwargs...)

AV.certify_real(F::HC.System; near=nothing, plane_points=nothing, kwargs...) =
    solution_counter(F, :certify_real; near, plane_points, kwargs...)

AV.dietmaier_function(F::HC.System; near=nothing, plane_points=nothing, kwargs...) =
    solution_counter(F, :dietmaier; near, plane_points, kwargs...)

"""
    visualize(F::HomotopyContinuation.System; func=:real, near=nothing,
              plane_points=nothing, kwargs...) -> (TriangulationCache, GLMakie.Figure)

Visualize solution counts or an imaginary-part diagnostic on an affine parameter plane.
Load both AdaptiveVisualization and HomotopyContinuation to enable this method.

`func` defaults to `:real` (numerical real counts) and also accepts `:certify_real`
(soft-certified real counts), `:positive` (numerical strictly positive real counts),
and `:dietmaier` (minimum nonzero imaginary L1 norm). The latter uses
`imaginary_zero_atol=1e-10` and returns zero when no norm exceeds that tolerance.
`plane_points=[p,q,r]` overrides `near` and maps plot coordinates `(u,v)` to
`p + zoomer*(u*(q-p) + v*(r-p))`. Without `plane_points`, two-parameter systems
keep their original axis directions, centered at `near` or zero. Larger parameter
spaces use random orthonormal directions through `near`, or the origin
if omitted. At least two parameters are required. `xlims` and `ylims` refer to
plot coordinates. The returned cache retains this map in `TC.parameter_slice`;
`retrieve_witnesses(TC)` returns sampled witnesses in the original parameter
coordinates for discrete caches.

The legend title defaults to `n_real` for `:real` and `:certify_real`, and
`n_pos` for `:positive`. Pass `legend_title` to override it.

The `:certify_real` mode provides soft certificates for floating-point input,
not guarantees about an intended exact model or input rounding. See
[`AdaptiveVisualization.certify_real`](@ref) and [`AdaptiveVisualization.real_solution_function`](@ref) for numerical options.
Remaining keywords are forwarded to the ordinary function-based `visualize`.
"""
function AV.visualize(F::HC.System; func=:real, near=nothing, plane_points=nothing,
        zoomer=1.0, rng::AbstractRNG=default_rng(), start_parameters=nothing,
        start_solutions=nothing, max_retries=nothing, real_tol=1e-8, positivity_tol=1e-7,
        imaginary_zero_atol=1e-10, wildcard_parity_mismatch=true,
        solver_options=(;), monodromy_options=(;),
        certification_options=(;), verbose=false, kwargs...)
    func in (:real, :certify_real, :positive, :dietmaier) ||
        throw(ArgumentError("func must be :real, :certify_real, :positive, or :dietmaier."))
    oracle = solution_counter(F, func; near, plane_points, zoomer, rng,
        start_parameters, start_solutions, max_retries, real_tol, positivity_tol,
        imaginary_zero_atol, wildcard_parity_mismatch, solver_options, monodromy_options,
        certification_options, verbose)
    default_title = func === :positive ? "n_pos" :
        func in (:real, :certify_real) ? "n_real" : "value"
    plot_kwargs = merge((; legend_title=default_title), (; kwargs...))
    return AV.visualize(oracle; verbose, plot_kwargs...)
end


#==============================================================================#
#                                 EXAMPLES                                     #
#==============================================================================#
# Example systems and an opt-in visualization helper follow. Loading the
# extension defines them; it never solves a system or opens an example figure.
# After loading both packages, kuramoto_model is available directly:
#
#   visualize(kuramoto_model(3))
#
# Other examples are accessed through:
#   HCExamples = Base.get_extension(AdaptiveVisualization, :HomotopyContinuationExt)
#==============================================================================#

function kuramoto_model(n)
    @var w[1:(n-1)], s[1:n], c[1:n]
    equations = []
    for i in 1:(n-1)
        coupling_sum = 0
        for j in 1:n
            coupling_sum += s[i] * c[j] - s[j] * c[i]
        end
        f_1 = w[i] - (1 / n) * coupling_sum
        f_2 = c[i]^2 + s[i]^2 - 1
        f_1 = HC.subs(f_1, [s[n], c[n]] => [0, 1])
        f_2 = HC.subs(f_2, [s[n], c[n]] => [0, 1])
        f_1 == 0 || push!(equations, f_1)
        f_2 == 0 || push!(equations, f_2)
    end
    return HC.System(equations; variables=[s[1:n-1]..., c[1:n-1]...], parameters=[w[1:n-1]...])
end

"""
    TwentySevenLines()

Construct the system of lines on an affine cubic surface, parameterized by its
20 coefficients. A generic cubic surface has 27 complex lines.
"""
function TwentySevenLines()
    @var x y z
    @var a[1:4, 1:4, 1:4]

    terms = []
    for i in 0:3
        for j in 0:3
            for k in 0:3
                i + j + k <= 3 && push!(terms, [i, j, k])
            end
        end
    end

    f = sum(a[c[1]+1, c[2]+1, c[3]+1] * x^c[1] * y^c[2] * z^c[3] for c in terms)
    Params = [a[c[1]+1, c[2]+1, c[3]+1] for c in terms]

    @var t b[1:2] c[1:2]

    lx = t
    ly = b[1] * t + b[2]
    lz = c[1] * t + c[2]

    g = HC.subs(f, [x, y, z] => [lx, ly, lz])
    Eqs = HC.coefficients(g, [t])

    return HC.System(Eqs; variables=[b[1], b[2], c[1], c[2]], parameters=Params)
end

# Plane conics in P^3 meeting points and lines; see Hauenstein-Sottile,
# "Algorithm XXX: alphaCertified", Section 5.3.
const SPACE_CONIC_DEGREES = Dict(
    (npoints=0, nlines=8) => 92,
    (npoints=1, nlines=6) => 18,
    (npoints=2, nlines=4) => 4,
    (npoints=3, nlines=2) => 1,
    (npoints=4, nlines=0) => 0,
)

"""
    SpaceConics(; npoints, nlines)
    SpaceConics(nlines, npoints)

Construct the incidence system for plane conics in projective three-space
meeting `npoints` points and `nlines` lines, with `2npoints + nlines == 8`.
The parameters describe the unfixed incidence conditions.
"""
function SpaceConics(; npoints::Int, nlines::Int)
    @assert 2npoints + nlines == 8 "Plane conics in P^3 have dimension 8, so need 2npoints + nlines == 8."
    @assert 0 <= npoints <= 4 "This Schubert problem has npoints in 0:4."
    @assert 0 <= nlines <= 8 "This Schubert problem has nlines in 0:8."

    @var A B C
    @var c11 c22 c12 c13 c23

    fixed_points = npoints >= 4 ? 4 : 0
    fixed_lines = fixed_points == 0 && nlines >= 3 ? 3 : 0
    free_points = npoints - fixed_points
    free_lines = nlines - fixed_lines

    # Kill the ambient projective symmetry in the parameter space by fixing
    # standard general-position incidence conditions whenever this Schubert
    # problem contains enough of them.
    e = [[i == j ? 1 : 0 for i in 1:4] for j in 1:4]

    # Free point parameters: p[:, i] is the i-th non-normalized point in P^3.
    @var p[1:4, 1:free_points]

    # Free line parameters: the j-th non-normalized line is spanned by u[:, j] and v[:, j].
    @var u[1:4, 1:free_lines] v[1:4, 1:free_lines]

    points = vcat(e[1:fixed_points], [p[:, i] for i in 1:free_points])
    # Use skew lines that are compatible with the conic chart q = ... + x3^2.
    # In particular, avoid fixing span(e3,e4), which would force the incidence
    # point [0:0:1] where this chart cannot see a conic.
    fixed_line_us = fixed_lines == 0 ? [] : [e[1], e[2], e[1] + e[2] + e[3]]
    fixed_line_vs = fixed_lines == 0 ? [] : [e[3], e[4], e[1] - e[2] + e[4]]
    line_us = vcat(fixed_line_us, [u[:, j] for j in 1:free_lines])
    line_vs = vcat(fixed_line_vs, [v[:, j] for j in 1:free_lines])

    variables = [A, B, C, c11, c22, c12, c13, c23]
    parameters = vcat(vec(p), vec(u), vec(v))

    q(x1, x2, x3) =
        c11*x1^2 +
        c22*x2^2 +
        x3^2 +
        c12*x1*x2 +
        c13*x1*x3 +
        c23*x2*x3

    h(x1, x2, x3, x4) = x4 - A*x1 - B*x2 - C*x3

    eqs = []

    for i in 1:npoints
        point = points[i]
        push!(eqs, h(point[1], point[2], point[3], point[4]))
        push!(eqs, q(point[1], point[2], point[3]))
    end

    for j in 1:nlines
        line_u = line_us[j]
        line_v = line_vs[j]

        hu = h(line_u[1], line_u[2], line_u[3], line_u[4])
        hv = h(line_v[1], line_v[2], line_v[3], line_v[4])

        z1 = hv*line_u[1] - hu*line_v[1]
        z2 = hv*line_u[2] - hu*line_v[2]
        z3 = hv*line_u[3] - hu*line_v[3]

        push!(eqs, q(z1, z2, z3))
    end

    return HC.System(eqs; variables=variables, parameters=parameters)
end

# Compatibility with the original positional convention: lines first, points second.
SpaceConics(nlines::Int, npoints::Int) = SpaceConics(; npoints=npoints, nlines=nlines)

"""
    expected_space_conic_degree(; npoints, nlines)
    expected_space_conic_degree(nlines, npoints)

Return the generic complex solution count for the corresponding conic problem.
"""
expected_space_conic_degree(; npoints::Int, nlines::Int) =
    SPACE_CONIC_DEGREES[(npoints=npoints, nlines=nlines)]

expected_space_conic_degree(nlines::Int, npoints::Int) =
    expected_space_conic_degree(; npoints=npoints, nlines=nlines)

"""Explore a random affine slice through the family of cubic surfaces."""
function run_twenty_seven_lines_example(; total_resolution=1000, kwargs...)
    return AV.visualize(TwentySevenLines();
        xlims=[-5, 5],
        ylims=[-5, 5],
        total_resolution=total_resolution,
        strategy=:sierpinski,
        rng=MersenneTwister(27),
        kwargs...,
    )
end

#=
Small polynomial example: run these commands explicitly in a Julia session.

using AdaptiveVisualization
using HomotopyContinuation

@var x y a b
F = System([x^2 + y^2 - a, x - y + b^3]; variables=[x, y], parameters=[a, b])
TC, fig = visualize(F;
    plane_points=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    xlims=[-1, 3], ylims=[-1, 1],
    total_resolution=144,
    xlabel="a", ylabel="b",
)
=#

end # module HomotopyContinuationExt
