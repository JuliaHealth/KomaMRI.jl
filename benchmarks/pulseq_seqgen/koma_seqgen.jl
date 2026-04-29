using KomaMRIBase
using KomaMRIFiles
using LinearAlgebra
using Statistics

const NTR = parse(Int, get(ENV, "SEQGEN_NTR", "100000"))
const NS = parse(Int, get(ENV, "SEQGEN_NS", "128"))
const REPEATS = parse(Int, get(ENV, "SEQGEN_REPEATS", "3"))
const WARMUP_NTR = parse(Int, get(ENV, "SEQGEN_WARMUP_NTR", "1"))
const OUTDIR = get(ENV, "SEQGEN_OUTDIR", joinpath(tempdir(), "koma_seqgen_bench"))
mkpath(OUTDIR)

const sys = Scanner(Gmax=40e-3, Smax=170.0, B1=Inf, ADC_Δt=100e-9, DUR_Δt=10e-6, GR_Δt=10e-6, RF_Δt=1e-6)
const rf = RF(1e-6, 0.5e-3)
const radial_gx = Grad(1e-3, 1e-3, 10e-6)
const spoiler_gx = Grad(-0.5e-3, 0.5e-3, 10e-6)
const radial_adc = ADC(64, (64 - 1) * 10e-6, 5e-6)

const spiral_t = range(0, 1; length=NS)
const spiral_amp = 2e-5
const spiral_gx = Grad(collect(spiral_amp .* sin.(2π .* 4 .* spiral_t) .* spiral_t), (NS - 1) * 10e-6, 0.0, 0.0, 0.0, 0.0, 0.0)
const spiral_gy = Grad(collect(spiral_amp .* cos.(2π .* 4 .* spiral_t) .* spiral_t), (NS - 1) * 10e-6, 0.0, 0.0, 0.0, 0.0, 0.0)
const spiral_adc = ADC(NS, (NS - 1) * 10e-6, 5e-6)

function radial_sequence(n=NTR)
    tr = Sequence(sys)
    @addblock tr += (rf, z=spoiler_gx) + (x=radial_gx, radial_adc) + (x=spoiler_gx)
    seq = Sequence(sys)
    @addblocks for k in 0:n-1
        θ = π * k / n
        phase = isodd(k) ? -1 + 0im : 1 + 0im
        seq += phase * (rotz(θ) * tr)
    end
    return seq
end

function spiral_sequence(n=NTR)
    tr = Sequence(sys)
    @addblock tr += (rf, z=spoiler_gx) + (x=spiral_gx, y=spiral_gy, spiral_adc) + (x=spoiler_gx)
    seq = Sequence(sys)
    @addblocks for k in 0:n-1
        θ = 2π * k / n
        phase = isodd(k) ? -1 + 0im : 1 + 0im
        seq += phase * (rotz(θ) * tr)
    end
    return seq
end

function run_case(name, build)
    warm = build(WARMUP_NTR)
    warmfile = joinpath(OUTDIR, "koma_$(name)_warmup.seq")
    write_seq(warm, warmfile; sys, verbose=false)
    read_seq(warmfile; verbose=false)

    construct_times = Float64[]
    write_times = Float64[]
    read_times = Float64[]
    blocks = 0
    bytes = 0
    for rep in 1:REPEATS
        GC.gc()
        push!(construct_times, @elapsed seq = build(NTR))
        filename = joinpath(OUTDIR, "koma_$(name)_$(NTR)tr_$(NS)ns_rep$(rep).seq")
        rm(filename; force=true)
        GC.gc()
        push!(write_times, @elapsed write_seq(seq, filename; sys, verbose=false))
        blocks = length(seq)
        bytes = filesize(filename)
        seq = nothing
        GC.gc()
        push!(read_times, @elapsed read_seq(filename; verbose=false))
    end
    t_construct = median(construct_times)
    t_write = median(write_times)
    t_read = median(read_times)
    println("koma,$name,construct,$NTR,$blocks,$NS,$t_construct,$bytes")
    println("koma,$name,write_seq,$NTR,$blocks,$NS,$t_write,$bytes")
    println("koma,$name,read_seq,$NTR,$blocks,$NS,$t_read,$bytes")
end

println("framework,case,operation,ntr,blocks,ns,seconds,bytes")
run_case("radial", radial_sequence)
run_case("spiral", spiral_sequence)
