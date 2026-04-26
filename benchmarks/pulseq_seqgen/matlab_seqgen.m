pulseq_path = getenv('PULSEQ_MATLAB_PATH');
if isempty(pulseq_path)
    candidates = { ...
        '/tmp/pulseq-dev/matlab', ...
        fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'pulseq', 'matlab') ...
    };
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'dir')
            pulseq_path = candidates{i};
            break
        end
    end
end
if isempty(pulseq_path) || ~exist(pulseq_path, 'dir')
    error('Set PULSEQ_MATLAB_PATH to the Pulseq MATLAB folder, e.g. /tmp/pulseq-dev/matlab');
end
addpath(genpath(pulseq_path));
warning('off', 'all');

NTR = str2double(getenv('SEQGEN_NTR'));
if isnan(NTR), NTR = 100000; end
NS = str2double(getenv('SEQGEN_NS'));
if isnan(NS), NS = 128; end
REPEATS = str2double(getenv('SEQGEN_REPEATS'));
if isnan(REPEATS), REPEATS = 3; end
WARMUP_NTR = str2double(getenv('SEQGEN_WARMUP_NTR'));
if isnan(WARMUP_NTR), WARMUP_NTR = 1; end
OUTDIR = getenv('SEQGEN_OUTDIR');
if isempty(OUTDIR), OUTDIR = fullfile(tempdir, 'koma_seqgen_bench'); end
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', 'MaxSlew', 170, 'SlewUnit', 'T/m/s');
rf = mr.makeBlockPulse(1e-6, 'Duration', 0.5e-3, 'system', sys);
radial_gx = mr.makeTrapezoid('x', 'Amplitude', 1e5, 'FlatTime', 1e-3, 'RiseTime', 10e-6, 'FallTime', 10e-6, 'system', sys);
spoiler_gx = mr.makeTrapezoid('x', 'Amplitude', -5e4, 'FlatTime', 0.5e-3, 'RiseTime', 10e-6, 'FallTime', 10e-6, 'system', sys);
radial_adc = mr.makeAdc(64, 'Dwell', 10e-6, 'Delay', 0, 'system', sys);

spiral_t = linspace(0, 1, NS);
spiral_amp = 2e3;
spiral_gx = mr.makeArbitraryGrad('x', spiral_amp * sin(2*pi*4*spiral_t) .* spiral_t, 'first', 0, 'last', 0, 'system', sys);
spiral_gy = mr.makeArbitraryGrad('y', spiral_amp * cos(2*pi*4*spiral_t) .* spiral_t, 'first', 0, 'last', 0, 'system', sys);
spiral_adc = mr.makeAdc(NS, 'Dwell', 10e-6, 'Delay', 0, 'system', sys);

fprintf('framework,case,operation,ntr,blocks,ns,seconds,bytes\n');
run_case('radial', @radial_sequence, NTR, NS, REPEATS, WARMUP_NTR, OUTDIR, sys, rf, radial_gx, spoiler_gx, radial_adc, spiral_gx, spiral_gy, spiral_adc);
run_case('spiral', @spiral_sequence, NTR, NS, REPEATS, WARMUP_NTR, OUTDIR, sys, rf, radial_gx, spoiler_gx, radial_adc, spiral_gx, spiral_gy, spiral_adc);

function seq = radial_sequence(NTR, sys, rf, radial_gx, spoiler_gx, radial_adc, ~, ~, ~)
    seq = mr.Sequence(sys);
    for k = 1:NTR
        phase = pi * mod(k - 1, 2);
        adc = radial_adc;
        adc.phaseOffset = phase;
        rf_k = rf;
        rf_k.phaseOffset = phase;
        seq.addBlock(rf_k, spoiler_gx);
        ev = mr.rotate('z', pi * (k - 1) / NTR, radial_gx, adc);
        seq.addBlock(ev{:});
        ev = mr.rotate('z', pi * (k - 1) / NTR, spoiler_gx);
        seq.addBlock(ev{:});
    end
end

function seq = spiral_sequence(NTR, sys, rf, ~, spoiler_gx, ~, spiral_gx, spiral_gy, spiral_adc)
    seq = mr.Sequence(sys);
    for k = 1:NTR
        phase = pi * mod(k - 1, 2);
        adc = spiral_adc;
        adc.phaseOffset = phase;
        rf_k = rf;
        rf_k.phaseOffset = phase;
        seq.addBlock(rf_k, spoiler_gx);
        ev = mr.rotate('z', 2*pi * (k - 1) / NTR, spiral_gx, spiral_gy, adc);
        seq.addBlock(ev{:});
        ev = mr.rotate('z', 2*pi * (k - 1) / NTR, spoiler_gx);
        seq.addBlock(ev{:});
    end
end

function run_case(name, build, NTR, NS, REPEATS, WARMUP_NTR, OUTDIR, sys, rf, radial_gx, spoiler_gx, radial_adc, spiral_gx, spiral_gy, spiral_adc)
    warm = build(WARMUP_NTR, sys, rf, radial_gx, spoiler_gx, radial_adc, spiral_gx, spiral_gy, spiral_adc);
    warmfile = fullfile(OUTDIR, sprintf('matlab_dev_%s_warmup.seq', name));
    warm.write(warmfile, true);
    warm2 = mr.Sequence(sys);
    warm2.read(warmfile);
    clear warm warm2;

    construct_times = zeros(REPEATS, 1);
    write_times = zeros(REPEATS, 1);
    read_times = zeros(REPEATS, 1);
    blocks = 0;
    bytes = 0;
    for rep = 1:REPEATS
        tic;
        seq = build(NTR, sys, rf, radial_gx, spoiler_gx, radial_adc, spiral_gx, spiral_gy, spiral_adc);
        construct_times(rep) = toc;
        blocks = length(seq.blockEvents);
        filename = fullfile(OUTDIR, sprintf('matlab_dev_%s_%dtr_%dns_rep%d.seq', name, NTR, NS, rep));
        if exist(filename, 'file'), delete(filename); end
        tic;
        seq.write(filename, true);
        write_times(rep) = toc;
        info = dir(filename);
        bytes = info.bytes;
        clear seq;
        tic;
        seq2 = mr.Sequence(sys);
        seq2.read(filename);
        read_times(rep) = toc;
        clear seq2;
    end
    t_construct = median(construct_times);
    t_write = median(write_times);
    t_read = median(read_times);
    fprintf('matlab-dev,%s,construct,%d,%d,%d,%.9f,%d\n', name, NTR, blocks, NS, t_construct, bytes);
    fprintf('matlab-dev,%s,write_seq,%d,%d,%d,%.9f,%d\n', name, NTR, blocks, NS, t_write, bytes);
    fprintf('matlab-dev,%s,read_seq,%d,%d,%d,%.9f,%d\n', name, NTR, blocks, NS, t_read, bytes);
end
