# Magnus Paper Examples

## Speed/Accuracy Benchmarks

`0_acc_bench.jl`: HS pulse, 400,001 native waveform samples, 100k spins, Magnus1/2/4/6, CPU Float64 Magnus6 reference.
Benchmark files save NRMSE, median runtime, runtime standard deviation, and raw timing samples.

```sh
cd examples/6.magnus_paper/speed_acc_benchmarks
julia --project=.. --threads=auto 1_run_convergence.jl
julia --project=.. --threads=auto 1_run_metal.jl 100
julia --project=.. --threads=auto 1_run_cuda.jl 100
julia --project=.. 2_fig_acc.jl
```

Output: `output/convergence_acc_results.jld2`, `output/metal_acc_results.jld2`, `output/cuda_acc_results.jld2`, `output/fig_acc.svg`.
