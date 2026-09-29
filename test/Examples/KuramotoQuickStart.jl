using AdaptiveVisualization
using HomotopyContinuation

K = kuramoto_model(6)
plane_points = [
    zeros(5),
    [-1/3, -1/3, 0.0, 1/3, 1/3],
    [1/3, 1/3, 0.0, 1/3, 1/3],
]
TC, fig = visualize(K;
    func=:real,
    plane_points=plane_points,
    xlims=[-1.5, 1.5],
    ylims=[-1.5, 1.5],
    strategy=:sierpinski,
    verbose=false,
    total_resolution=576,
    initial_resolution=144,
)
