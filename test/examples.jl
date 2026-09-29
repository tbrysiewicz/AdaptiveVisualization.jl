
using AdaptiveVisualization

#Disk: make sure construction/initialization works
f(x, y) = x^2 + y^2 < 1 ? 1 : 0
TC = TriangulationCache(f;
    xlims = [-2, 2],
    ylims = [-2, 2],
    resolution = 100,
    verbose = false,
)

#Visualizatoin works with some kwargs given
fig = visualize(TC; buttons = false)
fig = visualize(TC; buttons = true, edges = true)
#Display runs
display(fig)

#User defined complete function works. Triangles are given mean values of vertices
TC = TriangulationCache(f;
           xlims = [-2, 2],
           ylims = [-2, 2],
           resolution = 100,
           verbose = false,
           is_complete=(V,vals)->true,
           )


#Non-numerical values work
f(x,y) = x^2+y^2< 1 ? "inside" : "outside"
TC = TriangulationCache(f;
           xlims = [-3, 3],
           ylims = [-2, 2],
           verbose = true,
           )
visualize(TC, edges = true)


#continuous values work
f(x,y) = x^2+(y^2)*sin(x)
TC = TriangulationCache(f;
           xlims = [-2, 2],
           ylims = [-2, 2],
           verbose = true,
           )
visualize(TC, edges = true)
