using Test
using AdaptiveVisualization

@testset "AdaptiveVisualization" begin
    @testset "Discrete/continuous detection" begin
        discrete_values = Any[1, 1, 2, :wildcard]
        continuous_values = [float(i) for i in 1:60]
        many_string_values = ["value-$i" for i in 1:60]

        @test is_discrete(discrete_values)
        @test !is_discrete(continuous_values)
        @test is_discrete(many_string_values)
        @test AdaptiveVisualization.values_are_complete(Any["a", :wildcard, "a"])
        @test !AdaptiveVisualization.values_are_complete(Any["a", :wildcard, "b"])
    end

    @testset "Batched TriangulationCache construction" begin
        batch_calls = Ref(0)
        batch_sizes = Int[]
        f(points) = begin
            batch_calls[] += 1
            push!(batch_sizes, length(points))
            [p[1] + p[2] for p in points]
        end

        TC = TriangulationCache(f;
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=16,
            verbose=false,
        )

        @test AdaptiveVisualization.dimension(TC) == 2
        @test length(AdaptiveVisualization.function_values(TC)) == 16
        @test batch_calls[] == 1
        @test first(batch_sizes) == 16
        @test length(complete_triangles(TC)) + length(incomplete_triangles(TC)) > 0
        @test AdaptiveVisualization.remaining_oracle_budget(TC) === nothing
    end

    @testset "Discrete witnesses" begin
        batch_calls = Ref(0)
        value_at(p) = p[1] < 0 ? :left : (p[2] < 0 ? :wildcard : :right)
        TC = TriangulationCache(points -> begin
            batch_calls[] += 1
            value_at.(points)
        end; resolution=9)

        @test TC.parameter_slice === nothing
        witnesses = retrieve_witnesses(TC)
        @test value_at.(witnesses) == [:left, :right]
        @test witnesses == [[-1.0, -1.0], [0.0, 0.0]]
        @test batch_calls[] == 1

        # Retrieval includes cached values even after the visible window moves.
        AdaptiveVisualization.translate_window!(TC, 2.0, 0.0)
        @test retrieve_witnesses(TC) == witnesses
        @test batch_calls[] == 1

        continuous = TriangulationCache((x, y) -> x + 10y; resolution=64)
        @test !is_discrete(continuous)
        @test_throws ArgumentError retrieve_witnesses(continuous)
    end

    @testset "Single-point oracle wrapping" begin
        TC = TriangulationCache((x, y) -> x^2 + y^2;
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            verbose=false,
        )

        @test length(AdaptiveVisualization.function_values(TC)) == 9
        @test all(value isa Real for value in AdaptiveVisualization.output_values(TC))
    end

    @testset "Wildcard completeness" begin
        TC = TriangulationCache(points -> fill(:wildcard, length(points));
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            verbose=false,
        )

        @test isempty(incomplete_triangles(TC))
        @test !isempty(complete_triangles(TC))
        @test isempty(retrieve_witnesses(TC))
    end

    @testset "Refinement" begin
        never_complete(vertices, values; kwargs...) = false
        TC = TriangulationCache(points -> [p[1] - p[2] for p in points];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            strategy=:barycenter,
            min_refinement_area=0.0,
            is_complete=never_complete,
            verbose=false,
        )

        before = length(AdaptiveVisualization.function_values(TC))
        inserted = refine!(TC; verbose=false)

        @test inserted > 0
        @test length(AdaptiveVisualization.function_values(TC)) > before
        @test TC.total_oracle_calls == length(AdaptiveVisualization.function_values(TC))
        @test !isempty(AdaptiveVisualization.candidate_triangles(TC))
    end

    @testset "Refinement call budget" begin
        never_complete(vertices, values; kwargs...) = false
        TC = TriangulationCache(points -> [p[1] - p[2] for p in points];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            strategy=:barycenter,
            min_refinement_area=0.0,
            is_complete=never_complete,
            verbose=false,
        )

        inserted = refine!(TC; budget=1, verbose=false)

        @test inserted == 1
        @test AdaptiveVisualization.remaining_oracle_budget(TC) === nothing

        TC.oracle_budget = 1
        inserted = refine!(TC; budget=2, verbose=false)

        @test inserted == 2
        @test AdaptiveVisualization.remaining_oracle_budget(TC) == 1
    end

    @testset "Iterative refinement by normalized area" begin
        never_complete(vertices, values; kwargs...) = false
        TC = TriangulationCache(points -> [0 for _ in points];
            xlims=[0, 2],
            ylims=[0, 3],
            resolution=9,
            strategy=:barycenter,
            min_refinement_area=1.0,
            is_complete=never_complete,
            verbose=false,
        )

        inserted = refine!(TC; min_refinement_area=0.1, verbose=false)

        @test inserted > 0
        @test TC.min_refinement_area == 0.1
        @test isapprox(AdaptiveVisualization.scaled_min_refinement_area(TC), 0.6)
        @test isempty(AdaptiveVisualization.candidate_triangles(TC))
        @test_throws ErrorException refine!(TC; min_refinement_area=0.0, verbose=false)
    end

    @testset "Minimum-area slider exponent mapping" begin
        TC = TriangulationCache(points -> [p[1] + p[2] for p in points];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            min_refinement_area=4e-4,
            verbose=false,
        )

        slider_range = AdaptiveVisualization.default_min_refinement_area_slider_range()

        @test collect(slider_range) == collect(2:6)
        @test AdaptiveVisualization.default_min_refinement_area_slider_start(TC, collect(slider_range)) == 3
        @test AdaptiveVisualization.min_refinement_area_from_slider_exponent(5) == 1e-5
        @test AdaptiveVisualization.min_refinement_area_slider_label(6) == "Min area: window * 1e-6"
    end

    @testset "Visualization refinement to minimum area" begin
        never_complete(vertices, values; kwargs...) = false
        TC, fig = visualize(points -> zeros(Int, length(points));
            xlims=[0, 2],
            ylims=[0, 3],
            total_resolution=9,
            initial_resolution=9,
            strategy=:barycenter,
            min_refinement_area=0.1,
            is_complete=never_complete,
            buttons=false,
            verbose=false,
        )

        @test fig isa AdaptiveVisualization.GLMakie.Figure
        @test TC.min_refinement_area == 0.1
        @test TC.total_oracle_calls > 9
        @test isempty(AdaptiveVisualization.candidate_triangles(TC))
    end

    @testset "Stable categorical color value order" begin
        TC = TriangulationCache(points -> [p[1] < 0 ? 5 : 9 for p in points];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            verbose=false,
        )

        order = AdaptiveVisualization.stable_plot_value_order!(TC, [5, 9])
        @test order == [5, 9]

        order = AdaptiveVisualization.stable_plot_value_order!(TC, [3, 5])
        @test order == [3, 5, 9]

        label_TC = TriangulationCache(points -> [p[1] < 0 ? "inside" : "outside" for p in points];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=9,
            verbose=false,
        )

        order = AdaptiveVisualization.stable_plot_value_order!(label_TC, ["border", "inside"])
        @test order == ["border", "inside", "outside"]

        rational_triples = [
            Rational{Int}[0, 2, 0],
            Rational{Int}[0, 1, 1],
        ]
        vertex_values = AdaptiveVisualization.vertex_plot_values(rational_triples)
        @test vertex_values == [
            rational_triples[1], rational_triples[1], rational_triples[1],
            rational_triples[2], rational_triples[2], rational_triples[2],
        ]
        @test all(value -> value isa Vector{Rational{Int}}, vertex_values)
        @test AdaptiveVisualization.value_label(rational_triples[1]) == "[0, 2, 0]"

        triple_TC = TriangulationCache(points -> [
                p[1] < 0 ? Rational{Int}[1, 0, 0] : Rational{Int}[0, 2, 0]
                for p in points
            ];
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=4,
            verbose=false,
        )
        triple_order = AdaptiveVisualization.stable_plot_value_order!(
            triple_TC,
            [Rational{Int}[0, 1, 1]],
        )
        @test triple_order == [
            Rational{Int}[0, 1, 1],
            Rational{Int}[0, 2, 0],
            Rational{Int}[1, 0, 0],
        ]

        mixed_order = Any[
            Rational{Int}[1, 0, 0],
            "unknown",
            Rational{Int}[0, 2, 0],
            "failed",
            Rational{Int}[0, 1, 1],
        ]
        AdaptiveVisualization.sort_numeric_values!(mixed_order)
        @test mixed_order == Any[
            Rational{Int}[0, 1, 1],
            Rational{Int}[0, 2, 0],
            Rational{Int}[1, 0, 0],
            "failed",
            "unknown",
        ]
    end

    @testset "Categorical complete triangle plotting" begin
        always_complete(vertices, values; kwargs...) = true
        TC = TriangulationCache(points -> fill("inside", length(points));
            xlims=[-1, 1],
            ylims=[-1, 1],
            resolution=4,
            is_complete=always_complete,
            verbose=false,
        )
        triangle = first(complete_triangles(TC))

        @test AdaptiveVisualization.triangle_plot_value(TC, triangle) == "inside"

        AdaptiveVisualization.function_values(TC)[triangle[1]] = "outside"
        @test_throws ErrorException AdaptiveVisualization.triangle_plot_value(TC, triangle)
    end

    @testset "Example visualization workflows" begin
        disk_indicator(x, y) = x^2 + y^2 < 1 ? 1 : 0
        TC = TriangulationCache(disk_indicator;
            xlims=[-2, 2],
            ylims=[-2, 2],
            resolution=25,
            verbose=false,
        )

        fig = visualize(TC; buttons=false)
        @test fig isa AdaptiveVisualization.GLMakie.Figure

        fig = visualize(TC; buttons=true, edges=true)
        @test fig isa AdaptiveVisualization.GLMakie.Figure

        complete_TC = TriangulationCache(disk_indicator;
            xlims=[-2, 2],
            ylims=[-2, 2],
            resolution=25,
            verbose=false,
            is_complete=(vertices, values) -> true,
        )
        @test isempty(incomplete_triangles(complete_TC))
        @test !isempty(complete_triangles(complete_TC))

        region_label(x, y) = x^2 + y^2 < 1 ? "inside" : "outside"
        categorical_TC = TriangulationCache(region_label;
            xlims=[-3, 3],
            ylims=[-2, 2],
            resolution=25,
            verbose=false,
        )
        @test categorical_TC.min_refinement_area == 1e-5
        @test is_discrete(AdaptiveVisualization.output_values(categorical_TC))

        fig = visualize(categorical_TC; buttons=false, edges=true)
        @test fig isa AdaptiveVisualization.GLMakie.Figure
    end

    @testset "Invalid strategy" begin
        @test_throws Exception TriangulationCache((x, y) -> x + y; strategy=:quadtree, verbose=false)
        @test_throws Exception TriangulationCache((x, y) -> x + y; initial_resolution=9, verbose=false)
        @test_throws Exception TriangulationCache((x, y) -> x + y; total_resolution=9, verbose=false)
        @test_throws Exception TriangulationCache((x, y) -> x + y; min_refinement_area=-1, verbose=false)
        @test_throws Exception TriangulationCache((x, y) -> x + y; min_refinement_area=NaN, verbose=false)
    end
end

using Random
using LinearAlgebra
using HomotopyContinuation

include("regressions.jl")

@testset "HomotopyContinuation integration" begin
    extension = Base.get_extension(AdaptiveVisualization, :HomotopyContinuationExt)
    @test extension !== nothing
    @test nparameters(kuramoto_model(3)) == 2

    @var x y z a b c
    F = System([x^2 - a, y - b]; variables=[x, y], parameters=[a, b])
    F3 = System([x^2 - a, y - b, z - c];
        variables=[x, y, z], parameters=[a, b, c])
    identity_plane = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

    # Known start solutions avoid stochastic monodromy discovery in most tests.
    # A nonreal start parameter also keeps paths to negative a away from a = 0.
    start_parameters = ComplexF64[1 + im, 1 + im]
    start_solutions = [
        ComplexF64[sqrt(1 + im), 1 + im],
        ComplexF64[-sqrt(1 + im), 1 + im],
    ]
    starts = (; start_parameters, start_solutions, max_retries=0)

    @testset "Affine parameter plane" begin
        coordinates = [[0.0, 0.0], [1.0, -2.0], [3.5, 4.0]]
        for seed in (11, 27)
            default_slice = extension.parameter_slice(F; rng=MersenneTwister(seed))
            @test extension.parameter_values(default_slice, coordinates) == coordinates
            centered_slice = extension.parameter_slice(F;
                near=[2.0, 3.0], rng=MersenneTwister(seed))
            @test extension.parameter_values(centered_slice, coordinates) ==
                [p + [2.0, 3.0] for p in coordinates]
        end
        plane_points = [[2.0, 3.0], [4.0, 3.0], [3.0, 6.0]]
        slice = extension.parameter_slice(F;
            near=[91.0, 92.0], plane_points, zoomer=0.5, rng=MersenneTwister(1))
        coordinates = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [-1.0, 2.0]]
        expected = [[2.0, 3.0], [3.0, 3.0], [2.5, 4.5], [2.0, 6.0]]
        @test all(isapprox.(extension.parameter_values(slice, coordinates), expected))
        @test isempty(extension.parameter_values(slice, Vector{Float64}[]))

        center = [2.0, 3.0, 5.0]
        slice1 = extension.parameter_slice(F3;
            near=center, plane_points=nothing, zoomer=0.25, rng=MersenneTwister(19))
        slice2 = extension.parameter_slice(F3;
            near=center, plane_points=nothing, zoomer=0.25, rng=MersenneTwister(19))
        mapped1 = extension.parameter_values(slice1, coordinates[1:3])
        mapped2 = extension.parameter_values(slice2, coordinates[1:3])
        @test mapped1 == mapped2
        @test mapped1[1] == center
        directions = hcat(mapped1[2] - center, mapped1[3] - center)
        @test directions' * directions ≈ 0.25^2 * Matrix{Float64}(I, 2, 2)

        random_slice1 = extension.parameter_slice(F3;
            near=nothing, plane_points=nothing, zoomer=1.0, rng=MersenneTwister(29))
        random_slice2 = extension.parameter_slice(F3;
            near=nothing, plane_points=nothing, zoomer=1.0, rng=MersenneTwister(29))
        @test extension.parameter_values(random_slice1, coordinates) ==
              extension.parameter_values(random_slice2, coordinates)
    end

    @testset "Parameter validation" begin
        one_parameter = System([x^2 - a]; variables=[x], parameters=[a])
        @test_throws ArgumentError extension.parameter_slice(one_parameter;
            near=nothing, plane_points=nothing, zoomer=1.0, rng=MersenneTwister(1))
        @test_throws ArgumentError extension.parameter_slice(F;
            near=[1.0], plane_points=nothing, zoomer=1.0, rng=MersenneTwister(1))
        @test_throws ArgumentError extension.parameter_slice(F;
            near=[Inf, 1.0], plane_points=nothing, zoomer=1.0, rng=MersenneTwister(1))
        for bad_plane in (
            [[0.0, 0.0], [1.0, 0.0]],
            [[0.0, 0.0], [1.0], [0.0, 1.0]],
            [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]],
        )
            @test_throws ArgumentError extension.parameter_slice(F;
                near=nothing, plane_points=bad_plane, zoomer=1.0, rng=MersenneTwister(1))
        end
        @test_throws ArgumentError extension.parameter_slice(F;
            near=nothing, plane_points=identity_plane, zoomer=NaN, rng=MersenneTwister(1))
        @test_throws ArgumentError visualize(F; func=:unsupported)
        @test_throws ArgumentError real_solution_function(F;
            plane_points=identity_plane, starts..., real_tol=-1.0)
        @test_throws ArgumentError real_solution_function(F;
            plane_points=identity_plane, start_parameters, start_solutions, max_retries=-1)
    end

    @testset "Numerical and soft-certified counts" begin
        points = [[4.0, 2.0], [-1.0, 2.0], [4.0, -2.0], [4.0, 0.0], [9.0, 3.0]]
        expected_real = [2, 0, 2, 2, 2]
        expected_positive = [1, 0, 0, 0, 1]

        @test real_solution_function(F; starts...)(points) == expected_real
        for evaluator in (real_solution_function, certify_real, positive_solution_function, dietmaier_function)
            centered = evaluator(F; near=[1.0, 2.0], starts...)
            identity = evaluator(F; plane_points=identity_plane, starts...)
            @test centered([[3.0, 0.0], [-2.0, 0.0]]) == identity([[4.0, 2.0], [-1.0, 2.0]])
        end

        # This option would be rejected by HC.certify. Numerical modes must never call it.
        unused_certification_options = (this_is_not_a_certification_option=true,)
        numerical = real_solution_function(F;
            plane_points=identity_plane, starts...,
            certification_options=unused_certification_options)
        positive = positive_solution_function(F;
            plane_points=identity_plane, starts...,
            certification_options=unused_certification_options)
        certified = certify_real(F;
            plane_points=identity_plane, starts...,
            certification_options=(max_precision=256,))

        @test numerical isa Function
        @test positive isa Function
        @test certified isa Function
        @test numerical(points) == expected_real
        @test positive(points) == expected_positive
        @test certified(points) == expected_real
        @test numerical(reverse(points)) == reverse(expected_real)
        @test isempty(numerical(Vector{Float64}[]))
        @test isempty(positive(Vector{Float64}[]))
        @test isempty(certified(Vector{Float64}[]))

        # At a = 0, the double root cannot satisfy the nonsingular certification check.
        @test certified([[0.0, 2.0]]) == [:wildcard]
        retry_rng = MersenneTwister(41)
        retrying_certified = certify_real(F;
            plane_points=identity_plane, start_parameters, start_solutions,
            max_retries=1, rng=retry_rng,
            certification_options=(max_precision=256,))
        mixed_points = [[4.0, 2.0], [0.0, 2.0], [-1.0, 2.0]]
        draws_before_retry = rand(copy(retry_rng), UInt32, 4)
        @test retrying_certified(mixed_points) == Any[2, :wildcard, 0]
        draws_after_retry = rand(copy(retry_rng), UInt32, 4)
        @test draws_after_retry != draws_before_retry
        @test retrying_certified(mixed_points) == Any[2, :wildcard, 0]
        # A second evaluation reuses its retry start rather than drawing another one.
        @test rand(copy(retry_rng), UInt32, 4) == draws_after_retry

        # An odd positive count is valid even though the system has degree two.
        @test positive([[4.0, 2.0]]) == [1]
        @test positive([[4.0, 5e-8], [4.0, 2e-7]]) == [0, 1]
        thresholded = positive_solution_function(F;
            plane_points=identity_plane, starts..., positivity_tol=2.5)
        @test thresholded([[4.0, 3.0], [9.0, 3.0]]) == [0, 1]

        translated = real_solution_function(F;
            near=[100.0, 100.0],
            plane_points=[[4.0, 2.0], [5.0, 2.0], [4.0, 3.0]], starts...)
        @test translated([[0.0, 0.0], [-5.0, 0.0]]) == [2, 0]

        diagnostic = dietmaier_function(F; plane_points=identity_plane, starts...)
        @test diagnostic([[4.0, 2.0], [-1.0, 2.0]]) ≈ [0.0, 1.0]
    end

    @testset "Witnesses in original parameter coordinates" begin
        plane_points = [[0.25, 2.0, 5.0], [1.25, 2.0, 7.0], [0.25, 5.0, 6.0]]
        numerical = real_solution_function(F3; plane_points, zoomer=2.0,
            near=[99.0, 99.0, 99.0], max_retries=0,
            start_parameters=vcat(start_parameters, 2 + im),
            start_solutions=[vcat(s, 2 + im) for s in start_solutions])
        @test_throws ArgumentError TriangulationCache(numerical;
            input_dimension=2, resolution=4)
        @test_throws ArgumentError visualize(numerical;
            input_dimension=2, initial_resolution=4, total_resolution=4, buttons=false)
        TC = TriangulationCache(numerical; resolution=4, min_refinement_area=0.0)
        slice = TC.parameter_slice
        @test slice === numerical.parameter_slice
        @test slice([0.0, 0.0]) == plane_points[1]
        @test slice([1.0, -1.0]) == [2.25, -4.0, 7.0]
        @test retrieve_witnesses(TC) == [[-1.75, -4.0, -1.0], [2.25, -4.0, 7.0]]

        refine!(TC; budget=2)
        @test TC.parameter_slice === slice
        witnesses = retrieve_witnesses(TC)
        # These systems have two real solutions precisely when a > 0.
        @test [p[1] > 0 ? 2 : 0 for p in witnesses] ==
              unique(AdaptiveVisualization.function_values(TC))
        mapped_points = slice.(AdaptiveVisualization.input_points(TC))
        @test all(p -> p in mapped_points, witnesses)

        random_oracle = real_solution_function(F3; rng=MersenneTwister(29),
            max_retries=0, start_parameters=vcat(start_parameters, 2 + im),
            start_solutions=[vcat(s, 2 + im) for s in start_solutions])
        random_TC = TriangulationCache(random_oracle; resolution=4)
        @test random_TC.parameter_slice([0.0, 0.0]) == zeros(3)
        random_parameters = random_TC.parameter_slice.(AdaptiveVisualization.input_points(random_TC))
        @test [p[1] > 0 ? 2 : 0 for p in random_parameters] ==
              AdaptiveVisualization.function_values(random_TC)
        @test [p[1] > 0 ? 2 : 0 for p in retrieve_witnesses(random_TC)] ==
              unique(AdaptiveVisualization.function_values(random_TC))
    end

    @testset "Automatic start solution preparation" begin
        discovered = real_solution_function(F;
            plane_points=identity_plane, max_retries=0,
            rng=MersenneTwister(53), monodromy_options=(seed=UInt32(53),))
        @test discovered([[4.0, 2.0], [-1.0, 2.0]]) == [2, 0]
    end

    @testset "Dietmaier visualization smoke test" begin
        # For a < 0, x = ±im*sqrt(-a) and y = b, so the imaginary L1 norm is sqrt(-a).
        for imaginary_zero_atol in (1e-10, 3.0)
            TC, fig = visualize(F;
                func=:dietmaier, imaginary_zero_atol, starts...,
                xlims=[-4.0, -1.0], ylims=[1.0, 4.0],
                total_resolution=4, initial_resolution=4, buttons=false,
                certification_options=(this_is_not_a_certification_option=true,))
            expected = [sqrt(-p[1]) for p in AdaptiveVisualization.input_points(TC)]
            expected = [value > imaginary_zero_atol ? value : 0.0 for value in expected]
            @test fig isa AdaptiveVisualization.GLMakie.Figure
            @test length(expected) == 4
            @test all(isapprox.(AdaptiveVisualization.function_values(TC), expected))
        end
    end

    @testset "System visualization dispatch" begin
        # The whole displayed square has a > 0 and b > 0, so no refinement is needed.
        visual_options = (;
            plane_points=[[4.0, 2.0], [5.0, 2.0], [4.0, 3.0]],
            starts...,
            xlims=[0.0, 1.0], ylims=[0.0, 1.0],
            initial_resolution=4, total_resolution=4,
            buttons=false, verbose=false,
        )
        for (mode, expected) in ((:real, 2), (:certify_real, 2), (:positive, 1), (:dietmaier, 0.0))
            TC, fig = visualize(F; func=mode, visual_options...)
            @test fig isa AdaptiveVisualization.GLMakie.Figure
            @test length(AdaptiveVisualization.function_values(TC)) == 4
            @test all(==(expected), AdaptiveVisualization.function_values(TC))
            @test TC.parameter_slice([0.0, 0.0]) == [4.0, 2.0]
            @test retrieve_witnesses(TC) == [[4.0, 2.0]]
            legend = only(filter(x -> x isa AdaptiveVisualization.GLMakie.Legend, fig.content))
            @test legend.entrygroups[][1][1] == (mode === :positive ? "n_pos" :
                mode === :dietmaier ? "value" : "n_real")
        end
        _, custom_fig = visualize(F; func=:real, visual_options..., legend_title="custom")
        custom_legend = only(filter(x -> x isa AdaptiveVisualization.GLMakie.Legend,
            custom_fig.content))
        @test custom_legend.entrygroups[][1][1] == "custom"
        positive_TC, _ = visualize(F; func=:positive, plane_points=identity_plane,
            starts..., xlims=[4.0, 5.0], ylims=[5e-8, 2e-7],
            initial_resolution=4, total_resolution=4, buttons=false)
        @test sort(unique(AdaptiveVisualization.function_values(positive_TC))) == [0, 1]
        # A default visualization must not invoke certification.
        TC, fig = visualize(F; visual_options...,
            certification_options=(this_is_not_a_certification_option=true,))
        @test fig isa AdaptiveVisualization.GLMakie.Figure
        @test all(==(2), AdaptiveVisualization.function_values(TC))
    end
end
