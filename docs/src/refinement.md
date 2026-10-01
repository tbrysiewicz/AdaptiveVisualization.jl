# Sampling and refinement

```@meta
CurrentModule = AdaptiveVisualization
```

Construct a cache to sample a function without opening a figure. Refinement
adds samples where triangles are incomplete according to the cache's predicate.

```julia
TC = TriangulationCache((x, y) -> x^2 + y^2; resolution=100)
inserted = refine!(TC; budget=200, min_refinement_area=1e-4)
```

`budget` limits new samples; `min_refinement_area` sets the normalized triangle
area cutoff. When both are given, refinement stops when the budget is spent or
no eligible points remain. The return value is the number of inserted samples.

## Sampling cache

```@docs
TriangulationCache
```

## Refining a cache

```@docs
refine!
```

## Witnesses for discrete values

Retrieve one sampled input for each distinct non-wildcard value. For sliced
functions, witnesses use the original input coordinates.

```@docs
retrieve_witnesses
```
