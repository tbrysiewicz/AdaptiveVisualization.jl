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
| `:dietmaier` | Minimum nonzero imaginary L1 norm |

The constructors below create evaluators independently of plotting. Each accepts
a batch of two-dimensional slice coordinates and returns one value per point.

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
