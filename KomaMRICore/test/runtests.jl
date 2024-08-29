using TestItems, TestItemRunner

### NOTE: by default, tests are run on the CPU with the number of threads set to
#   Threads.nthreads(). To run on a specific GPU backend, add the name of the
#   backend package ("AMDGPU", "CUDA", "Metal", or "oneAPI") to the test/Project.toml
#   file in KomaMRICore and pass the name as a test argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; test_args=["CUDA"])
#
#   To run on the cpu with a specific number of threads, pass the number of threads
#   as a julia argument.
#
#   Example:
#
#   import Pkg
#   Pkg.test("KomaMRICore"; julia_args=`--threads=4`)
#
#   For changing the default backend used for testing,
#   modify the [preferences.KomaMRICore] section in the test/Project.toml:
#
#     [preferences.KomaMRICore]
#     test_backend = "CPU"
#
#   For the backend preference to take effect, you need to:
#   - REPL testing: No action needed. `] test` should pick up the preference right away.
#   - VSCode testing: You need to restart VSCode.
#
#   Sadly, LocalPreferences.toml are not picked up by VScode (that could be .gitignore'd),
#   so we had put them into the test/Project.toml.
#
###

group     = get(ENV, "TEST_GROUP", :all) |> Symbol
test_file = get(ENV, "TEST_FILE", :none) |> Symbol

# if we are testing just a single file then group = :none
# to skip the full test suite
if test_file != :none
    group = :none
end

@testset "KomaMRICore" begin
    if test_file != :none
        @testset "Single file test" begin
            include(String(test_file))
        end
    end
    
    if group == :core || group == :all
        @testset "Core" begin
            include("test_core.jl")
        end
    end

    if group == :motion || group == :all
        @testset "Motion" begin
            include("test_motion.jl")
        end
    end
end

#Environment variable set by CI
const CI = get(ENV, "CI", nothing)

@run_package_tests filter=ti->(:core in ti.tags)&&(isnothing(CI) || :skipci ∉ ti.tags) #verbose=true

