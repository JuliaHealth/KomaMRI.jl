using TestItems, TestItemRunner

@run_package_tests filter=ti->!(:skipci in ti.tags)&&(:plots in ti.tags) #verbose=true

#GUI tests
include("GUI_PlotlyBase_backend_test.jl")
include("GUI_PlutoPlotly_backend_test.jl")
