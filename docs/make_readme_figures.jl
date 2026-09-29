using AdaptiveVisualization
using GLMakie
using HomotopyContinuation
using LinearAlgebra
using Random

const ASSET_DIR = joinpath(@__DIR__, "assets")
const README_FIGURE_RESOLUTION = 144
const README_FIGURE_REFINEMENT_PASSES = 5
const KURAMOTO_GRID_RESOLUTIONS = [5^2, 10^2, 40^2, 50^2, 100^2]
const KURAMOTO_GRID_RS = 0:5
mkpath(ASSET_DIR)

function save_readme_figure(fig, filename)
    path = joinpath(ASSET_DIR, filename)
    GLMakie.save(path, fig; px_per_unit=2)
    println("Saved ", path)
    return path
end

function refine_readme_figure!(TC)
    for _ in 1:README_FIGURE_REFINEMENT_PASSES
        refine!(TC; verbose=false)
    end
    return TC
end

function disk_indicator_figure()
    f(x, y) = x^2 + y^2 < 1 ? 1 : 2
    TC = TriangulationCache(f;
        xlims=[-2, 2],
        ylims=[-2, 2],
        resolution=README_FIGURE_RESOLUTION,
        strategy=:sierpinski,
        verbose=false,
    )
    refine_readme_figure!(TC)
    return visualize(TC;
        buttons=false,
        edges=true,
        figure_size=(900, 700),
        legend_title="region",
        show_legend=false,
    )
end

function categorical_figure()
    f(x, y) = x^2 + y^2 < 1 ? "inside" : "outside"
    TC = TriangulationCache(f;
        xlims=[-2, 2],
        ylims=[-2, 2],
        resolution=README_FIGURE_RESOLUTION,
        strategy=:sierpinski,
        verbose=false,
    )
    refine_readme_figure!(TC)
    return visualize(TC;
        buttons=false,
        edges=true,
        figure_size=(900, 700),
        legend_title="class",
        show_legend=false,
    )
end

function continuous_quadratic_figure()
    f(x, y) = x^2 + y^2 - x
    TC = TriangulationCache(f;
        xlims=[-2, 2],
        ylims=[-2, 2],
        resolution=README_FIGURE_RESOLUTION,
        strategy=:sierpinski,
        verbose=false,
    )
    refine_readme_figure!(TC)
    return visualize(TC;
        buttons=false,
        edges=false,
        figure_size=(900, 700),
        legend_title="value",
        show_legend=true,
        xlabel="",
        ylabel="",
        title="Continuous scalar field",
    )
end

function kuramoto_figure(; resolution=README_FIGURE_RESOLUTION, refinement_passes=README_FIGURE_REFINEMENT_PASSES, title=nothing)
    Random.seed!(3)
    F = kuramoto_model(3)
    f = real_solution_function(F;
        plane_points=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    )
    TC = TriangulationCache(f;
        xlims=[-1, 1],
        ylims=[-1, 1],
        resolution=resolution,
        strategy=:sierpinski,
        min_refinement_area=0,
        verbose=false,
    )
    for _ in 1:refinement_passes
        refine!(TC; verbose=false)
    end
    figure_title = if title === :oracle_summary
        "($(resolution),$(refinement_passes)) - total oracle evaluations = $(length(AdaptiveVisualization.function_values(TC)))"
    elseif title === nothing
        ""
    else
        title
    end
    return visualize(TC;
        buttons=false,
        edges=false,
        figure_size=(900, 700),
        legend_title="nreal",
        show_legend=false,
        xlabel="",
        ylabel="",
        title=figure_title,
        titlesize=27,
    )
end

function dietmaier_kuramoto_figure(;
        resolution=72,
        refinement_passes=4)
    Random.seed!(3)
    F = kuramoto_model(3)
    f = dietmaier_function(F;
        plane_points=[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    )
    TC = TriangulationCache(f;
        xlims=[-1, 1],
        ylims=[-1, 1],
        resolution=resolution,
        strategy=:sierpinski,
        min_refinement_area=0,
        verbose=false,
    )
    for _ in 1:refinement_passes
        refine!(TC; verbose=false)
    end
    return visualize(TC;
        buttons=false,
        edges=false,
        figure_size=(900, 700),
        legend_title="minimum nonzero imaginary L1 norm",
        show_legend=true,
        xlabel="",
        ylabel="",
        title="Kuramoto: Dietmaier Function",
        plot_log_transform=true,
    )
end

function kuramoto_grid()
    paths = String[]
    for resolution in KURAMOTO_GRID_RESOLUTIONS
        for r in KURAMOTO_GRID_RS
            filename = "kuramoto-n$(resolution)-r$(r).png"
            push!(paths, save_readme_figure(kuramoto_figure(;
                resolution=resolution,
                refinement_passes=r,
                title=:oracle_summary,
            ), filename))
        end
    end
    return paths
end

const FIGURE_BUILDERS = Dict(
    "disk" => () -> save_readme_figure(disk_indicator_figure(), "disk-indicator.png"),
    "categorical" => () -> save_readme_figure(categorical_figure(), "categorical-inside-outside.png"),
    "continuous" => () -> save_readme_figure(continuous_quadratic_figure(), "continuous-quadratic.png"),
    "kuramoto" => () -> save_readme_figure(kuramoto_figure(), "kuramoto-real-solutions.png"),
    "dietmaier-kuramoto" => () -> save_readme_figure(dietmaier_kuramoto_figure(), "dietmaier-kuramoto.png"),
    "kuramoto-grid" => kuramoto_grid,
)

const DEFAULT_README_FIGURES = ["disk", "continuous", "dietmaier-kuramoto", "kuramoto-grid"]

requested_figures = isempty(ARGS) ? DEFAULT_README_FIGURES : ARGS
for name in requested_figures
    haskey(FIGURE_BUILDERS, name) || error("Unknown README figure '$name'. Use one of $(sort(collect(keys(FIGURE_BUILDERS)))).")
    FIGURE_BUILDERS[name]()
end
