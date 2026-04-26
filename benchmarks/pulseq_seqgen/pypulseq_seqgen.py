import copy
import math
import os
import statistics
import time

import numpy as np
import pypulseq as pp

NTR = int(os.environ.get("SEQGEN_NTR", "100000"))
NS = int(os.environ.get("SEQGEN_NS", "128"))
REPEATS = int(os.environ.get("SEQGEN_REPEATS", "3"))
WARMUP_NTR = int(os.environ.get("SEQGEN_WARMUP_NTR", "1"))
OUTDIR = os.environ.get("SEQGEN_OUTDIR", os.path.join("/tmp", "koma_seqgen_bench"))
os.makedirs(OUTDIR, exist_ok=True)

system = pp.Opts(max_grad=40, grad_unit="mT/m", max_slew=170, slew_unit="T/m/s")
rf = pp.make_block_pulse(flip_angle=1e-6, duration=0.5e-3, system=system)
radial_gx = pp.make_trapezoid("x", amplitude=1e4, flat_time=1e-3, rise_time=10e-6, fall_time=10e-6, system=system)
spoiler_gx = pp.make_trapezoid("x", amplitude=-5e3, flat_time=0.5e-3, rise_time=10e-6, fall_time=10e-6, system=system)
radial_adc = pp.make_adc(64, dwell=10e-6, delay=0.0, system=system)

spiral_t = np.linspace(0, 1, NS)
spiral_amp = 2e2
spiral_gx = pp.make_arbitrary_grad("x", spiral_amp * np.sin(2 * np.pi * 4 * spiral_t) * spiral_t, first=0.0, last=0.0, system=system)
spiral_gy = pp.make_arbitrary_grad("y", spiral_amp * np.cos(2 * np.pi * 4 * spiral_t) * spiral_t, first=0.0, last=0.0, system=system)
spiral_adc = pp.make_adc(NS, dwell=10e-6, delay=0.0, system=system)


def phased(event, phase):
    out = copy.copy(event)
    out.phase_offset = phase
    return out


def radial_sequence(n=NTR):
    seq = pp.Sequence(system)
    for k in range(n):
        phase = math.pi if k % 2 else 0.0
        seq.add_block(phased(rf, phase), spoiler_gx)
        seq.add_block(*pp.rotate(radial_gx, phased(radial_adc, phase), angle=math.pi * k / n, axis="z", system=system))
        seq.add_block(*pp.rotate(spoiler_gx, angle=math.pi * k / n, axis="z", system=system))
    return seq


def spiral_sequence(n=NTR):
    seq = pp.Sequence(system)
    for k in range(n):
        phase = math.pi if k % 2 else 0.0
        seq.add_block(phased(rf, phase), spoiler_gx)
        seq.add_block(*pp.rotate(spiral_gx, spiral_gy, phased(spiral_adc, phase), angle=2 * math.pi * k / n, axis="z", system=system))
        seq.add_block(*pp.rotate(spoiler_gx, angle=2 * math.pi * k / n, axis="z", system=system))
    return seq


def run_case(name, build):
    warm = build(WARMUP_NTR)
    warmfile = os.path.join(OUTDIR, f"pypulseq_{name}_warmup.seq")
    warm.write(warmfile, create_signature=True, remove_duplicates=True, check_timing=True)
    pp.Sequence(system).read(warmfile, remove_duplicates=True)

    construct_times = []
    write_times = []
    read_times = []
    blocks = 0
    bytes_ = 0
    for rep in range(REPEATS):
        t0 = time.perf_counter()
        seq = build(NTR)
        construct_times.append(time.perf_counter() - t0)
        blocks = len(seq.block_events)
        filename = os.path.join(OUTDIR, f"pypulseq_{name}_{NTR}tr_{NS}ns_rep{rep + 1}.seq")
        if os.path.exists(filename):
            os.remove(filename)
        t0 = time.perf_counter()
        seq.write(filename, create_signature=True, remove_duplicates=True, check_timing=True)
        write_times.append(time.perf_counter() - t0)
        bytes_ = os.path.getsize(filename)
        seq = None
        t0 = time.perf_counter()
        seq2 = pp.Sequence(system)
        seq2.read(filename, remove_duplicates=True)
        read_times.append(time.perf_counter() - t0)
    t_construct = statistics.median(construct_times)
    t_write = statistics.median(write_times)
    t_read = statistics.median(read_times)
    print(f"pypulseq,{name},construct,{NTR},{blocks},{NS},{t_construct:.9f},{bytes_}", flush=True)
    print(f"pypulseq,{name},write_seq,{NTR},{blocks},{NS},{t_write:.9f},{bytes_}", flush=True)
    print(f"pypulseq,{name},read_seq,{NTR},{blocks},{NS},{t_read:.9f},{bytes_}", flush=True)


print("framework,case,operation,ntr,blocks,ns,seconds,bytes")
run_case("radial", radial_sequence)
run_case("spiral", spiral_sequence)
