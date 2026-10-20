function generate_reference(pulseq_root)
% Regenerate with MATLAB and the official Pulseq checkout:
% addpath('test/test_files/gradient_spectrum'); generate_reference('/path/to/pulseq')
% The Julia tests read the saved results and do not require MATLAB or Pulseq.
addpath(fullfile(pulseq_root, 'matlab'));
fixture_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(fixture_dir)));
reference = struct();
reference.matlab_version = version;
[status, commit] = system(sprintf('git -C "%s" rev-parse HEAD', pulseq_root));
assert(status == 0, 'Cannot identify the Pulseq revision');
reference.pulseq_commit = strtrim(commit);

% Spiral exercises arbitrary gradients in a single unwindowed segment;
% EPI exercises 50% overlap, Hann windowing and RMS aggregation of segments.
for name = {'spiral', 'epi_ramp'}
    name = name{1};
    seq = mr.Sequence();
    seq.read(fullfile(repo_root, 'examples', '1.sequences', [name '.seq']));
    for fmax = [3000, 1200]
        [R, Rax, F] = seq.gradSpectrum([], fmax, false);
        reference.(name).(sprintf('f%d', fmax)) = struct('R', R, 'Rax', Rax, 'F', F);
    end
end

save(fullfile(fixture_dir, 'pulseq_reference.mat'), '-struct', 'reference', '-v7');
fprintf('Saved gradient spectrum references: MATLAB %s, Pulseq %s\n', ...
    reference.matlab_version, reference.pulseq_commit);
end
