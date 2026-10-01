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

Use `discrete_legend=true` for a categorical legend or `false` for a continuous
colorbar. Omitting it selects automatically; `show_legend=false` hides either.
Use `edges=true` to show triangle boundaries, and
`figure_padding=(20, 20, 20, 40)` to set the left, right, bottom, and top margins. Set `buttons=false` to omit interactive controls.

## Methods

```@docs
visualize
```
