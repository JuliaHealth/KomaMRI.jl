using Printf

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, "..", ".."))
const NTR = get(ENV, "SEQGEN_NTR", "100000")
const NS = get(ENV, "SEQGEN_NS", "128")
const OUTDIR = get(ENV, "SEQGEN_OUTDIR", joinpath(tempdir(), "koma_seqgen_bench"))
const MATLAB = get(ENV, "MATLAB", "/Applications/MATLAB_R2025b.app/bin/matlab")
const RUN_MATLAB = get(ENV, "SEQGEN_RUN_MATLAB", "true") != "false"
const RUN_PYPULSEQ = get(ENV, "SEQGEN_RUN_PYPULSEQ", "true") != "false"

struct BenchRow
    framework::String
    case::String
    operation::String
    ntr::Int
    blocks::Int
    ns::Int
    seconds::Float64
    bytes::Int
end

function run_and_collect(cmd)
    env = copy(ENV)
    env["SEQGEN_NTR"] = NTR
    env["SEQGEN_NS"] = NS
    env["SEQGEN_OUTDIR"] = OUTDIR
    output = read(setenv(cmd, env), String)
    print(output)
    return parse_rows(output)
end

function parse_rows(output)
    rows = BenchRow[]
    for line in split(output, '\n')
        isempty(line) && continue
        startswith(line, "framework,") && continue
        fields = split(line, ',')
        length(fields) == 8 || continue
        push!(rows, BenchRow(
            fields[1],
            fields[2],
            fields[3],
            parse(Int, fields[4]),
            parse(Int, fields[5]),
            parse(Int, fields[6]),
            parse(Float64, fields[7]),
            parse(Int, fields[8]),
        ))
    end
    return rows
end

function add_totals(rows)
    out = copy(rows)
    for framework in unique(row.framework for row in rows), case_name in unique(row.case for row in rows)
        selected = [row for row in rows if row.framework == framework && row.case == case_name]
        length(selected) == 3 || continue
        push!(out, BenchRow(framework, case_name, "total", selected[1].ntr, selected[1].blocks, selected[1].ns, sum(row.seconds for row in selected), selected[1].bytes))
    end
    return out
end

function row_lookup(rows)
    return Dict((row.framework, row.case, row.operation) => row for row in rows)
end

function print_markdown_table(rows)
    rows = add_totals(rows)
    table = row_lookup(rows)
    println()
    println("| Case | Operation | Koma | PyPulseq | MATLAB dev | Koma vs PyPulseq | Koma vs MATLAB |")
    println("|---|---|---:|---:|---:|---:|---:|")
    for case_name in ("radial", "spiral"), operation in ("construct", "write_seq", "read_seq", "total")
        koma = table[("koma", case_name, operation)]
        py = get(table, ("pypulseq", case_name, operation), nothing)
        ml = get(table, ("matlab-dev", case_name, operation), nothing)
        py_s = isnothing(py) ? "NA" : @sprintf("%.3f s", py.seconds)
        ml_s = isnothing(ml) ? "NA" : @sprintf("%.3f s", ml.seconds)
        py_x = isnothing(py) ? "NA" : @sprintf("%.2fx", py.seconds / koma.seconds)
        ml_x = isnothing(ml) ? "NA" : @sprintf("%.2fx", ml.seconds / koma.seconds)
        println("| $case_name | `$operation` | $(@sprintf("%.3f s", koma.seconds)) | $py_s | $ml_s | $py_x | $ml_x |")
    end
end

mkpath(OUTDIR)
rows = BenchRow[]
koma_cmd = `$(Base.julia_cmd()) --project=$(joinpath(ROOT, "KomaMRIFiles", "test")) --threads=auto $(joinpath(HERE, "koma_seqgen.jl"))`
append!(rows, run_and_collect(koma_cmd))
if RUN_PYPULSEQ
    append!(rows, run_and_collect(`uv run --with pypulseq python $(joinpath(HERE, "pypulseq_seqgen.py"))`))
end
if RUN_MATLAB
    append!(rows, run_and_collect(`$MATLAB -batch $(string("run('", joinpath(HERE, "matlab_seqgen.m"), "')"))`))
end
print_markdown_table(rows)
