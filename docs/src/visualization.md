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

Pass `input_dimension` for a random two-dimensional plane through the origin,
`near` for a plane through a specified point, or `plane_points=[p, q, r]` to
choose the plane explicitly. These options also work with ordinary functions:

```julia
f(p) = sum(abs2, p)
TC, fig = visualize(f; input_dimension=4, batched=false)
```

`TC.parameter_slice([u, v])` maps plot coordinates back to the original input
space. `zoomer` scales the plane, and `rng` controls random plane selection.

## Figure options

Use `discrete_legend=true` for a categorical legend, `edges=true` to show
triangle boundaries, and `figure_padding=(20, 20, 20, 40)` to set the left,
right, bottom, and top margins. Set `buttons=false` to omit interactive controls.

## Methods

```@docs
visualize
```
