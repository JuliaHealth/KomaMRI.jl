# Pulseq Sequence-Generation Benchmarks

Benchmarks Pulseq-style sequence generation for repeated TRs:

- RF block
- rotated radial or spiral readout with ADC
- rotated spoiler
- alternating RF/ADC phase cycling
- write to `.seq`
- read the written `.seq`

Each case warms up construction, writing, and reading on a tiny sequence before
timing, then reports the median over `SEQGEN_REPEATS`.

Run the full table from the repo root:

```bash
julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
```

Useful options:

```bash
SEQGEN_NTR=100000 SEQGEN_NS=128 julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
SEQGEN_REPEATS=5 julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
SEQGEN_RUN_MATLAB=false julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
SEQGEN_RUN_PYPULSEQ=false julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
PULSEQ_MATLAB_PATH=/tmp/pulseq-dev/matlab julia --project=benchmarks benchmarks/pulseq_seqgen/run_suite.jl
```

Requirements:

- KomaMRI local checkout
- `uv` for the PyPulseq benchmark
- MATLAB R2025b by default at `/Applications/MATLAB_R2025b.app/bin/matlab`
- MATLAB Pulseq dev path via `PULSEQ_MATLAB_PATH`, or `/tmp/pulseq-dev/matlab`

Generated `.seq` files are written to `SEQGEN_OUTDIR`, defaulting to a temp
directory.
