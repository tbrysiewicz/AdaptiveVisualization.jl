# Compatibility entry point for existing notebooks and scripts.
# HC functionality and example systems live in the package extension.
using AdaptiveVisualization
using HomotopyContinuation

HCExamples = Base.get_extension(AdaptiveVisualization, :HomotopyContinuationExt)
using .HCExamples: kuramoto_model, TwentySevenLines, SpaceConics,
    expected_space_conic_degree, run_twenty_seven_lines_example, SPACE_CONIC_DEGREES

if get(ENV, "ADAPTIVE_SAMPLING_RUN_HC_EXAMPLES", "false") == "true"
    run_twenty_seven_lines_example()
end
