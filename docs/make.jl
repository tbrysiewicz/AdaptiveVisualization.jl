using AdaptiveVisualization
using HomotopyContinuation
using Documenter

# Load the extension so visualize(::System) is included in the API reference.
hc_extension = Base.get_extension(AdaptiveVisualization, :HomotopyContinuationExt)

# Documenter checks inclusion of existing docstrings; also catch missing docstrings.
for name in names(AdaptiveVisualization)
    name === :AdaptiveVisualization && continue
    Docs.hasdoc(AdaptiveVisualization, name) || error("Missing public docstring: $name")
end

makedocs(
    root=@__DIR__,
    sitename="AdaptiveVisualization.jl",
    modules=[AdaptiveVisualization, hc_extension],
    checkdocs=:public,
    # Explicit filenames support both local file browsing and web hosting.
    format=Documenter.HTML(prettyurls=false),
    pages=[
        "Introduction" => "index.md",
        "Getting started" => "getting_started.md",
        "Visualizing functions" => [
            "Visualization" => "visualization.md",
            "Sampling and refinement" => "refinement.md",
        ],
        "Polynomial systems" => "polynomial_systems.md",
        "Saving figures" => "saving.md",
    ],
)
