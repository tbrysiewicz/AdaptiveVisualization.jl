using Test
using Random
using AdaptiveVisualization

@testset "1.0.0 regressions" begin
    AV = AdaptiveVisualization
    M = AV.GLMakie

    @testset "Incomplete categorical regions stay blank" begin
        TC = TriangulationCache((x, y) -> x < 0 ? "left" : "right"; resolution=4)
        @test isempty(AdaptiveVisualization.complete_triangles(TC))
        @test !isempty(AdaptiveVisualization.incomplete_triangles(TC))
        @test isempty(AV.selected_triangles(TC, false))
        fig = M.Figure()
        ax = M.Axis(fig[1, 1])
        @test AV.draw_triangulation!(fig, ax, TC).plot === nothing
    end

    @testset "Explicit discrete mode" begin
        calls = Ref(0)
        count_at(x, y) = floor(Int, 4 * (x + 1)) + 9 * floor(Int, 4 * (y + 1))
        oracle(points) = begin
            calls[] += length(points)
            [count_at(p...) for p in points]
        end
        TC = TriangulationCache(oracle; resolution=289, discrete=false, strategy=:barycenter)
        before = calls[]
        continuous_fig = visualize(TC; buttons=false)
        @test any(x -> x isa M.Colorbar, continuous_fig.content)
        @test any(v -> !isinteger(v), TC.plot_value_order[false])

        fig = visualize(TC; discrete=true, buttons=false)
        @test calls[] == before
        @test AV.is_discrete(TC)
        @test length(retrieve_witnesses(TC)) == 81
        @test any(x -> x isa M.Legend, fig.content)
        @test Set(TC.plot_value_order[false]) == Set(AV.function_values(TC))
        @test !TC.is_complete(nothing, [42, 44, 44])
        @test TC.is_complete(nothing, Any[42, :wildcard, 42])
        @test TC.is_complete(nothing, fill(:wildcard, 3))
        @test !isempty(AV.incomplete_triangles(TC))
        @test all(T -> AV.triangle_plot_value(TC, T) === nothing,
            AV.incomplete_triangles(TC))
        all_fig = visualize(TC; plot_all_triangles=true, buttons=false)
        @test any(x -> x isa M.Legend, all_fig.content)
        @test all(isinteger, TC.plot_value_order[false])
        override_fig = visualize(TC; discrete_legend=false, buttons=false)
        @test any(x -> x isa M.Colorbar, override_fig.content)
        @test TC.discrete === true
        @test refine!(TC; budget=1) == 1
        @test AV.is_discrete(TC)
        @test !TC.is_complete(nothing, [42, 44, 44])

        visualize(TC; discrete=nothing, buttons=false)
        @test TC.discrete === nothing
        @test !AV.is_discrete(TC)
        @test TC.is_complete(nothing, [42, 44, 44])

        forced = TriangulationCache(oracle; resolution=289, discrete=true)
        @test AV.is_discrete(forced)
        @test !forced.is_complete(nothing, [42, 44, 44])
        @test length(retrieve_witnesses(forced)) == 81
        tagged = AV.SlicedOracle(oracle, AV.affine_parameter_slice(2), true)
        inherited = TriangulationCache(tagged; resolution=289)
        @test inherited.discrete === true
        @test length(retrieve_witnesses(inherited)) == 81
        @test !inherited.is_complete(nothing, [42, 44, 44])
        @test any(x -> x isa M.Legend, visualize(inherited; buttons=false).content)
        shown, shown_fig = visualize((x, y) -> x < 0 ? 42 : 44;
            discrete=true, initial_resolution=9, total_resolution=9, buttons=false)
        @test shown.discrete === true
        @test any(x -> x isa M.Legend, shown_fig.content)
        @test all(isinteger, shown.plot_value_order[false])

        smooth = TriangulationCache((x, y) -> x + y; resolution=9, discrete=false)
        @test !AV.is_discrete(smooth)
        @test smooth.is_complete(nothing, [0.0, 0.1, 0.2])
        @test_throws ArgumentError retrieve_witnesses(smooth)
        @test any(x -> x isa M.Colorbar, visualize(smooth; buttons=false).content)

        custom(vertices, values; kwargs...) = true
        custom_TC = TriangulationCache((x, y) -> x; resolution=4, is_complete=custom)
        visualize(custom_TC; discrete=true, buttons=false)
        @test custom_TC.is_complete === custom
        @test isempty(AV.incomplete_triangles(custom_TC))
        @test all(T -> AV.triangle_plot_value(custom_TC, T) === nothing,
            AV.complete_triangles(custom_TC))
        @test_throws ArgumentError TriangulationCache(oracle; discrete=:yes)
        @test_throws ArgumentError visualize(forced; discrete=1)
        @test forced.discrete === true
    end

    @testset "Visible legend and persistent colors" begin
        region(x, y) = x < -1 ? 3 : (x < 0 ? 5 : 19)
        TC = TriangulationCache(region; resolution=81, discrete=true)
        fig = M.Figure()
        ax = M.Axis(fig[1, 1])
        AV.set_axis_window!(ax, TC)
        drawn = Ref{Any}(AV.draw_triangulation!(fig, ax, TC))
        legend() = only(filter(x -> x isa M.Legend, fig.content))
        entries() = isempty(legend().entrygroups[]) ? [] : only(legend().entrygroups[])[2]
        labels() = [entry.label[] for entry in entries()]
        colors() = Dict(entry.label[] => M.RGBAf(first(entry.elements).attributes.polycolor[])
            for entry in entries())
        initial_colors = colors()
        @test labels() == ["5", "19"]
        initial_calls = TC.total_oracle_calls
        listener_count = length(M.Makie.listeners(ax.finallimits))

        M.xlims!(ax, -0.9, -0.3)
        @test labels() == ["5"]
        @test colors()["5"] == initial_colors["5"]
        M.xlims!(ax, 0.3, 0.9)
        @test labels() == ["19"]
        @test colors()["19"] == initial_colors["19"]
        M.xlims!(ax, 3, 4)
        @test isempty(labels())
        @test !legend().blockscene.visible[]
        M.xlims!(ax, -1, 1)
        @test labels() == ["5", "19"]
        @test legend().blockscene.visible[]
        @test colors() == initial_colors
        @test TC.total_oracle_calls == initial_calls

        AV.navigate_and_refine!(fig, ax, TC, drawn;
            translate=(-0.75, 0.0), navigation_initial_resolution=81,
            navigation_refinement_budget=0)
        @test labels() == ["3", "5"]
        @test colors()["5"] == initial_colors["5"]
        @test TC.plot_value_order[false] == [5, 19, 3]
        new_color = colors()["3"]
        @test Set(M.RGBAf.(drawn[].plot.color[])) == Set(values(colors()))
        @test length(M.Makie.listeners(ax.finallimits)) == listener_count

        AV.navigate_and_refine!(fig, ax, TC, drawn;
            translate=(0.75, 0.0), seed_and_refine=false)
        @test labels() == ["5", "19"]
        @test colors() == initial_colors
        AV.navigate_and_refine!(fig, ax, TC, drawn;
            translate=(-1.0, 0.0), zoom_factor=0.5, seed_and_refine=false)
        @test all(label -> label in ["3", "5"], labels())
        @test colors()["3"] == new_color
        @test length(M.Makie.listeners(ax.finallimits)) == listener_count

        triangle = ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        @test !AV.triangle_overlaps_window(triangle, (0.8, 1.0, 0.8, 1.0))
        @test AV.triangle_overlaps_window(triangle, (0.1, 0.2, 0.1, 0.2))
        @test AV.triangle_overlaps_window(triangle, (-1.0, 2.0, -1.0, 2.0))
        @test !AV.triangle_overlaps_window(triangle, (1.0, 2.0, 0.0, 1.0))
    end

    @testset "Edges button and navigation" begin
        for initial_edges in (true, false)
            TC = TriangulationCache((x, y) -> x < 0 ? 0 : 1; resolution=9)
            fig = visualize(TC; buttons=true, edges=initial_edges,
                min_refinement_area_controls=false,
                navigation_initial_resolution=4, navigation_refinement_budget=0)
            ax = only(filter(x -> x isa M.Axis, fig.content))
            button(label) = only(filter(x -> x isa M.Button && x.label[] == label, fig.content))
            count_edges() = count(x -> x isa M.Lines, ax.scene.plots)
            @test count_edges() == Int(initial_edges)
            edge_button = button("Edges")
            edge_button.clicks[] += 1
            @test count_edges() == Int(!initial_edges)
            right_button = button("→")
            right_button.clicks[] += 1
            @test count_edges() == Int(!initial_edges)
            edge_button.clicks[] += 1
            @test count_edges() == Int(initial_edges)
        end
    end

    @testset "Plot stays inside figure margins" begin
        TC = TriangulationCache((x, y) -> x + y; resolution=9)
        for (buttons, slider) in ((true, true), (true, false), (false, false))
            fig = visualize(TC; buttons, min_refinement_area_controls=slider,
                title="Refinement", show_legend=false)
            ax = only(filter(x -> x isa M.Axis, fig.content))
            bounds = ax.layoutobservables.computedbbox[]
            height = fig.scene.viewport[].widths[2]
            @test bounds.origin[2] >= 20
            @test bounds.origin[2] + bounds.widths[2] <= height - 40
        end
    end

    @testset "Saving creates the destination directory" begin
        TC = TriangulationCache((x, y) -> x + y; resolution=4)
        fig = visualize(TC; buttons=false)
        mktempdir() do directory
            cd(directory) do
                filename = AV.save(fig, "nested/figure"; file_extension="png", dpi=150)
                @test filename == "OutputFiles/nested/figure.png"
                @test isfile(filename)
            end
        end
    end

    @testset "Initialization and total evaluation budgets" begin
        points = AV.triangulation_initial_points([0.0, 100.0], [0.0, 1.0], 4)
        @test length(points) == 4
        @test length(unique(first.(points))) == 2
        @test length(unique(last.(points))) == 2
        TC = TriangulationCache(points -> fill(0, length(points));
            xlims=[0, 100], ylims=[0, 1], resolution=4)
        @test TC.total_oracle_calls == 4

        calls = Ref(0)
        oracle(points) = begin
            calls[] += length(points)
            fill(0, length(points))
        end
        never_complete(vertices, values; kwargs...) = false
        TC, _ = visualize(oracle; xlims=[0, 3], ylims=[0, 1],
            initial_resolution=9, total_resolution=10,
            is_complete=never_complete, strategy=:barycenter, buttons=false)
        @test TC.total_oracle_calls == 10
        @test calls[] == 10

        # An omitted total grows to accommodate an explicitly requested mesh.
        TC, _ = visualize(points -> fill(0, length(points));
            initial_resolution=1001, buttons=false)
        @test 4 <= TC.total_oracle_calls <= 1001

        rejected_calls = Ref(0)
        rejected(points) = begin
            rejected_calls[] += length(points)
            fill(0, length(points))
        end
        @test_throws ErrorException visualize(rejected;
            initial_resolution=1001, total_resolution=1000, buttons=false)
        @test rejected_calls[] == 0
        @test_throws ErrorException visualize(rejected;
            initial_resolution=1001, resolution=1000, buttons=false)
        @test rejected_calls[] == 0
    end

    @testset "Nonnumeric categories above the legend limit" begin
        label(x, y) = x < -0.25 ? "left" : (x > 0.25 ? "right" : "middle")
        TC = TriangulationCache(label; resolution=25)
        @test !isempty(AdaptiveVisualization.complete_triangles(TC))
        fig = M.Figure()
        ax = M.Axis(fig[1, 1])
        drawn = AV.draw_triangulation!(fig, ax, TC;
            legend_max_values=1, show_legend=false)
        @test drawn.plot !== nothing
        @test Tuple(drawn.plot.colorrange[]) == (1.0, 3.0)
    end

    @testset "Explicit oracle call form" begin
        plane = [[1.0, 2.0, 3.0, 4.0], [2.0, 2.0, 3.0, 4.0],
                 [1.0, 3.0, 3.0, 4.0]]
        TC = TriangulationCache(sum; plane_points=plane, resolution=4, batched=false)
        mapped = TC.parameter_slice.(AV.input_points(TC))
        expected = [sum(point) for point in mapped]
        @test AV.output_values(TC) == expected
        @test all(point -> length(point) == 4, mapped)

        shown, _ = visualize(sum; plane_points=plane, batched=false,
            initial_resolution=4, total_resolution=4, buttons=false)
        @test AV.output_values(shown) == expected

        batch_sizes = Int[]
        batch(points) = begin
            push!(batch_sizes, length(points))
            [sum(point) for point in points]
        end
        batch_TC = TriangulationCache(batch; plane_points=plane,
            resolution=4, batched=true)
        @test batch_sizes == [4]
        @test AV.output_values(batch_TC) == expected
    end

    @testset "Ordinary-function parameter slices" begin
        category(p) = p[1] > 0 ? :positive : :negative
        TC = TriangulationCache(category; input_dimension=3,
            rng=MersenneTwister(29), resolution=4)
        @test TC.parameter_slice([0.0, 0.0]) == zeros(3)
        mapped = TC.parameter_slice.(AV.input_points(TC))
        @test category.(mapped) == AV.output_values(TC)
        @test category.(retrieve_witnesses(TC)) == unique(AV.output_values(TC))
        @test all(p -> length(p) == 3, retrieve_witnesses(TC))
        repeated = TriangulationCache(category; input_dimension=3,
            rng=MersenneTwister(29), resolution=4)
        @test repeated.parameter_slice.basis == TC.parameter_slice.basis

        plane = [[1.0, 2.0, 3.0], [3.0, 2.0, 3.0], [1.0, 6.0, 3.0]]
        batch_calls = Ref(0)
        batch(points) = begin
            batch_calls[] += 1
            [sum(p) > 6 ? :high : :low for p in points]
        end
        explicit = TriangulationCache(batch; plane_points=plane,
            near=[99.0, 99.0, 99.0], zoomer=0.5, resolution=4)
        @test batch_calls[] == 1
        @test explicit.parameter_slice([1.0, 1.0]) == [2.0, 4.0, 3.0]
        @test all(p -> length(p) == 3, retrieve_witnesses(explicit))

        shifted = TriangulationCache(category; near=[1.0, 2.0, 3.0], resolution=4)
        @test shifted.parameter_slice([0.0, 0.0]) == [1.0, 2.0, 3.0]
        @test_throws ArgumentError TriangulationCache(category; input_dimension=1)

        shown, _ = visualize(category; input_dimension=3,
            rng=MersenneTwister(17), initial_resolution=4,
            total_resolution=4, buttons=false)
        @test shown.parameter_slice !== nothing
        @test all(p -> length(p) == 3, retrieve_witnesses(shown))

        plain_calls = Ref(0)
        plain_batch(points) = begin
            plain_calls[] += 1
            [p[1] + p[2] for p in points]
        end
        plain = TriangulationCache(plain_batch; resolution=4)
        @test plain.parameter_slice === nothing
        @test plain_calls[] == 1
    end
end
