# Getting started

## Installation

Requires Julia 1.11 or later. Once the package is registered:

```julia
using Pkg
Pkg.add("AdaptiveVisualization")
```

Before registration, or for development, use `Pkg.develop(path="/path/to/AdaptiveVisualization")`
with a local checkout. For polynomial-system examples, also install
HomotopyContinuation with `Pkg.add("HomotopyContinuation")`.
Interactive figures use GLMakie and require a working OpenGL display.

## Ordinary functions

```julia
using AdaptiveVisualization

TC, fig = visualize((x, y) -> x^2 + y^2; total_resolution=1000)
```

To sample without opening a figure, construct a cache and refine it directly:

```julia
TC = TriangulationCache((x, y) -> x^2 + y^2; resolution=100)
inserted = refine!(TC; budget=200, min_refinement_area=1e-4)
```

The budget limits new samples. Refinement can stop earlier when no eligible
points remain. See [`refine!`](@ref) for the stopping rules.

## Polynomial systems

Loading both packages activates the integration:

```julia
using AdaptiveVisualization, HomotopyContinuation

@var x b c
F = System([x^2 + b*x + c]; variables=[x], parameters=[b, c])
TC, fig = visualize(F; func=:real, xlims=[-3, 3], ylims=[-2, 3])
```

The `func` keyword selects numerical real counts (`:real`), strictly positive
real counts (`:positive`), soft-certified real counts (`:certify_real`), or an
imaginary-part diagnostic (`:dietmaier`). See [`visualize`](@ref) and the
[Polynomial systems](polynomial_systems.md) reference for the available evaluators.
