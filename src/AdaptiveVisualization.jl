#==============================================================================#
#  ADAPTIVE VISUALIZATION
#==============================================================================#

module AdaptiveVisualization

using DelaunayTriangulation: Triangulation, triangulate, each_solid_triangle, triangle_vertices, add_point!, get_point, num_points
import GLMakie
using LinearAlgebra: I, qr, rank
using Random: AbstractRNG, default_rng, randn

export
    visualize,
    TriangulationCache,
    refine!,
    retrieve_witnesses,
    real_solution_function,
    positive_solution_function,
    certify_real,
    dietmaier_function,
    kuramoto_model,
    save

"""
    real_solution_function(F; near=nothing, plane_points=nothing, kwargs...)

Construct a batched numerical real-solution counter for a `HomotopyContinuation.System`.
Load `HomotopyContinuation` to activate this method. No certification is performed.
Caches and visualizations built from this evaluator default to `discrete=true`.
Pass `discrete=false` to use continuous behavior instead.

All solution-counter constructors use the same rules for choosing a real
parameter slice. Three `plane_points=[p,q,r]` override `near` and map `(u,v)` to `p + zoomer*(u*(q-p) + v*(r-p))`.
Without `plane_points`, two-parameter systems keep their original axis directions,
centered at `near` or zero. Larger parameter spaces use two random orthonormal
directions through `near`, or the origin if omitted. Systems need at
least two parameters. Use `rng` to reproduce random slices; `zoomer=1.0` controls
scale in all cases. The returned evaluator carries its slice into caches built
by `TriangulationCache` or `visualize`, so [`retrieve_witnesses`](@ref) returns
parameters in the original coordinates.

The evaluator reuses generic start solutions and retries unresolved samples up to
`max_retries=2` times, returning `:wildcard` if they remain unresolved.
`real_tol=1e-8` is an absolute imaginary-part tolerance. Solver keywords belong in
`solver_options=(; ...)`, monodromy keywords in `monodromy_options=(; ...)`.
Pass `start_parameters` and `start_solutions` together to supply a known start fibre.

The returned evaluator accepts a vector of two-dimensional slice coordinates
and returns one value per point, in the same order. For example:

```julia
using AdaptiveVisualization, HomotopyContinuation
F = kuramoto_model(3)
evaluate = real_solution_function(F)
counts = evaluate([[0.1, 0.2], [0.2, 0.1]])
```
"""
function real_solution_function end

"""
    positive_solution_function(F; near=nothing, plane_points=nothing, kwargs...)

Construct a batched counter of real solutions whose every variable is strictly
greater than `positivity_tol=1e-7`. Uses the same slicing rules, retries, and numerical
realness tolerance as [`real_solution_function`](@ref), without certification.
Load `HomotopyContinuation` to activate this method.
Caches and visualizations built from this evaluator default to `discrete=true`.
Pass `discrete=false` to use continuous behavior instead.

The returned evaluator accepts a vector of two-dimensional slice coordinates
and returns one value per point, in the same order. For example:

```julia
using AdaptiveVisualization, HomotopyContinuation
F = kuramoto_model(3)
evaluate = positive_solution_function(F)
counts = evaluate([[0.1, 0.2], [0.2, 0.1]])
```
"""
function positive_solution_function end

"""
    certify_real(F; near=nothing, plane_points=nothing, kwargs...)

Construct a batched real-solution counter using HomotopyContinuation interval certification.
Caches and visualizations built from this evaluator default to `discrete=true`.
Pass `discrete=false` to use continuous behavior instead.
Uses the same slicing rules as [`real_solution_function`](@ref), with `max_retries=5`.
Each accepted sample accounts for the start fibre's solutions as distinct
certified real or certified nonreal solutions; unresolved samples are `:wildcard`.
Certification options belong in `certification_options=(; ...)`; the default
maximum precision is 16384 bits. Progress is disabled unless `verbose=true`.

These are **soft certificates** for the floating-point system and parameter
values supplied to HomotopyContinuation. They do not certify an intended exact model or robustness
to input rounding, nor independently prove completeness of the starting set.
Load `HomotopyContinuation` to activate this method.

The returned evaluator accepts a vector of two-dimensional slice coordinates
and returns one value per point, in the same order. For example:

```julia
using AdaptiveVisualization, HomotopyContinuation
F = kuramoto_model(3)
evaluate = certify_real(F)
counts = evaluate([[0.1, 0.2], [0.2, 0.1]])
```
"""
function certify_real end

"""
    dietmaier_function(F; near=nothing, plane_points=nothing, imaginary_zero_atol=1e-10, kwargs...)

Construct a batched imaginary-part diagnostic using the same slicing rules as
[`real_solution_function`](@ref). At each point, the evaluator returns the
smallest imaginary L1 norm strictly greater than `imaginary_zero_atol` among
the tracked solutions, or `0.0` if no norm exceeds that threshold. Unresolved
samples are `:wildcard`. Load `HomotopyContinuation` to activate this method.

The returned evaluator accepts a vector of two-dimensional slice coordinates
and returns one value per point, in the same order. For example:

```julia
using AdaptiveVisualization, HomotopyContinuation
F = kuramoto_model(3)
evaluate = dietmaier_function(F)
values = evaluate([[0.1, 0.2], [0.2, 0.1]])
```
"""
function dietmaier_function end

"""
    kuramoto_model(n)

Construct the polynomial Kuramoto equilibrium system for `n` oscillators, fixing
the last oscillator's phase (`s[n]=0`, `c[n]=1`). Return a
`HomotopyContinuation.System` with variables ordered as
`[s[1], …, s[n-1], c[1], …, c[n-1]]` and frequency parameters ordered as
`[w[1], …, w[n-1]]`. Parameter vectors supplied to the solution-counter helpers
must follow this order. The `n - 1` frequency parameters are real when
exploring real equilibria. Use `n >= 3` for a two-dimensional visualization.
Load `HomotopyContinuation` to activate this constructor.
"""
function kuramoto_model end

include("Utilities.jl")
include("TCStruct.jl")
include("Initialization.jl")
include("WindowChange.jl")
include("Refinement.jl")
include("MakieInterface.jl")

end # module AdaptiveVisualization
