# Polynomial systems

```@meta
CurrentModule = AdaptiveVisualization
```

Load HomotopyContinuation alongside AdaptiveVisualization to activate
polynomial-system methods. [`visualize`](@ref) selects its evaluator with `func`:

| `func` | Quantity |
| --- | --- |
| `:real` (default) | Numerical real-solution count |
| `:positive` | Numerical strictly positive real-solution count |
| `:certify_real` | Soft-certified real-solution count |
| `:dietmaier` | Smallest imaginary L1 norm above `imaginary_zero_atol`, or zero |

Solution counts default to `discrete=true`, including when you pass a counter
from `real_solution_function`, `positive_solution_function`, or `certify_real`
to `TriangulationCache` or `visualize`. This prevents averaging different counts
and keeps a categorical legend even with many distinct counts. Pass
`discrete=false` to use continuous behavior instead. Dietmaier is a continuous diagnostic and retains
automatic classification unless you supply `discrete`.

The constructors below return evaluators that can be used without plotting.
Each returned evaluator accepts a batch of two-dimensional slice coordinates
and returns one value per point.

## Numerical solution counts

```@docs
real_solution_function
positive_solution_function
```

## Certification

```@docs
certify_real
```

## Imaginary-part diagnostic

```@docs
dietmaier_function
```

## Example model

```@docs
kuramoto_model
```
