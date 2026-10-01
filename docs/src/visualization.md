# Visualization

```@meta
CurrentModule = AdaptiveVisualization
```

[`visualize`](@ref) accepts an ordinary function, a polynomial system, or an
existing [`TriangulationCache`](@ref). Functions and systems return `(TC, fig)`;
plotting an existing cache returns only the figure.

## Refinement strategy

Choose where new samples are placed with `strategy`: `:sierpinski` (the default)
uses edge midpoints, `:barycenter` uses centroids, and `:random` uses interior points.

```julia
TC, fig = visualize((x, y) -> x^2 + y^2; strategy=:barycenter)
```

## Slices of higher-dimensional functions

For an ordinary function, `input_dimension=n` tells `visualize` how many
coordinates each input point has. It chooses a random two-dimensional plane
through the origin in that input space. For example, each input to `f` below
has four coordinates, while the plot has two axes:

```julia
using AdaptiveVisualization

f(p) = sum(abs2, p)
TC, fig = visualize(f; input_dimension=4, batched=false)
```

Alternatively, supply `near=p` for a plane through `p`, or
`plane_points=[p, q, r]` to specify the plane. Their point lengths determine the
input dimension, so `input_dimension` can be omitted. HomotopyContinuation
systems infer it from their parameters automatically.

`TC.parameter_slice([u, v])` maps plot coordinates back to the original input
space. `zoomer` scales the plane, and `rng` controls random plane selection.

## Figure options

Use `discrete=true` for categorical values, including integer solution counts:

```julia
TC, fig = visualize((x, y) -> floor(Int, 4x); discrete=true)
```

This requires exact agreement of non-wildcard vertex values, hides incomplete
triangles by default, and uses a categorical legend regardless of how many
values occur. Mixed-value triangles are never colored with an averaged value
in this mode. A custom `is_complete` predicate still controls refinement.

Use `discrete=false` for continuous numeric values and a colorbar. Omitting
`discrete` keeps automatic detection. The setting is retained by the cache during
refinement and navigation. For an existing cache, `visualize(TC; discrete=true)`
updates its classification without resampling; `discrete=nothing` restores
the evaluator's default. Ordinary functions use automatic detection, while
HomotopyContinuation solution counters default to discrete behavior. An explicit
`discrete=false` overrides either default.

Categorical legends show only values on colored triangles in the current window.
They update when you zoom or pan, including with the mouse. Colors remain assigned
to the same values across window changes and refinement; newly discovered values
do not change existing assignments. Legend entries are sorted independently of
those assignments. Values found only on uncolored triangles are omitted.

`discrete_legend` remains available as a legend-only override: `true` forces a
categorical legend and `false` a colorbar. `show_legend=false` hides either.
Use `edges=true` to show triangle boundaries, and
`figure_padding=(20, 20, 20, 40)` to set the left, right, bottom, and top margins. Set `buttons=false` to omit interactive controls.

## Methods

```@docs
visualize
```
