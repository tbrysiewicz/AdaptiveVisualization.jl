# Saving figures

```@meta
CurrentModule = AdaptiveVisualization
```

Save the figure returned by [`visualize`](@ref):

```julia
AdaptiveVisualization.save(fig, "example"; dpi=300)
```

This writes `OutputFiles/example.png` relative to the working directory.

```@docs
AdaptiveVisualization.save
```
