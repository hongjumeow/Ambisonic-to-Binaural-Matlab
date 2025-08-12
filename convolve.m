%% Load Source to Convolve
[source, fs_src] = audioread('your-directory-to-source');

%% Convert source to mono if stereo
if size(source, 2) == 2
    source = mean(source, 2);
end

%% Resample source to brir's sampling rate
fs_brir = 32000;
if fs_src ~= fs_brir
    source = resample(source, fs_brir, fs_src);
end

%% Convolve with BRIRs
brir_files = dir('your-directory-to-converted-brirs\*.wav');
for i = 1:length(brir_files)
    % Load BRIR (stereo IR)
    [brir, fs_brir_check] = audioread(fullfile(brir_files(i).folder, brir_files(i).name));

    % Convolve with source signal
    left  = conv(source, brir(:,1), 'full');
    right = conv(source, brir(:,2), 'full');

    % Combine stereo result
    spatialized = [left, right];

    % Normalize
    spatialized = spatialized / max(abs(spatialized), [], 'all');

    % Save result
    [~, base, ~] = fileparts(brir_files(i).name);
    out_filename = fullfile('spatialized_output', [base, '_convolved.wav']);
    if ~exist('spatialized_output', 'dir')
        mkdir('spatialized_output');
    end
    audiowrite(out_filename, spatialized, fs_brir);
end