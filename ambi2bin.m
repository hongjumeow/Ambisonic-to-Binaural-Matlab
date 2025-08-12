addpath(genpath('your-directory-to\sfs-matlab'));
addpath(genpath('your-directory-to\Spherical-Harmonic-Transform'));
SFS_start();

%% Match Normalization
function Ysn3d = n3d_to_sn3d(Yn3d, order)
    Ysn3d = Yn3d;
    col = 1;
    for n = 0:order
        ncols = 2*n + 1;
        Ysn3d(:, col:col+ncols-1) = Ysn3d(:, col:col+ncols-1) / sqrt(2*n+1);
        col = col + ncols;
    end
end

%% Spherical Transform
% Define decoding grid
[az_grid, el_grid] = meshgrid(0:15:345, -30:30:60);  % az x el grid
az_list = az_grid(:);
el_list = el_grid(:);

% Convert to SH-compatible dirs
az_rad = deg2rad(az_list);
el_rad = deg2rad(el_list);
colat_rad = pi/2 - el_rad;
dirs = [wrapTo2Pi(az_rad), colat_rad];  % [1368 x 2]

az_deg_list = rad2deg(dirs(:,1));                 % [0, 360) azimuth in degrees
el_deg_list = rad2deg(pi/2 - dirs(:,2));          % convert colatitude to elevation

% SH order and basis
order = 8;
Y = getSH(order, dirs, 'real'); % Spherical Harmonic basis for all dirs

Y = n3d_to_sn3d(Y, order);

%% HRTF loading and preprocessing
sofa_filename = 'your-directory-to/RIEC_hrir_subject_001.sofa';
hrtf = sofaread(sofa_filename);
fs_hrtf = hrtf.SamplingRate;

az_cipic = deg2rad(hrtf.SourcePosition(:,1));   % azimuth in radians
el_cipic = deg2rad(hrtf.SourcePosition(:,2));   % elevation in radians
colat_cipic = pi/2 - el_cipic;
dirs_hrtf = [wrapTo2Pi(az_cipic), colat_cipic];  % for SH matching

%% Target directions for SARITA Upsample
% Not many more than # HOA
az_step = 10;
el_step = 10;
az_tgt = 0:az_step:350;           % [0, 350) deg
el_tgt = -30:el_step:60;          % [-30, 80] deg

[az_tgt_grid, el_tgt_grid] = meshgrid(az_tgt, el_tgt);
az_tgt_list = az_tgt_grid(:);
el_tgt_list = el_tgt_grid(:);

% Optional: collapse azimuth at the zenith (el = 90 deg), since azimuth is undefined there.
is_zenith = abs(el_tgt_list - 90) < 1e-9;
az_tgt_list(is_zenith) = 0;

% Convert to SH-compatible [az, colat] in radians
az_tgt_rad  = deg2rad(az_tgt_list);
el_tgt_rad  = deg2rad(el_tgt_list);
colat_tgt   = pi/2 - el_tgt_rad;
target_dirs = [wrapTo2Pi(az_tgt_rad), colat_tgt];   % size: [252 x 2]

% % Or Lebedev Sphere
% leb = getLebedevSphere(194); % Not many more than (N+1)^2 = 25
% [target_az, target_el, ~] = cart2sph(leb.x, leb.y, leb.z);
% target_dirs = [wrapTo2Pi(target_az), pi/2 - target_el];

%% Get all .h5 ambisonic files
data_root = 'your-directory-to-ambisonic-h5files';
h5_files = dir(fullfile(data_root, '**', '*.h5'));
fs_hoa = 32000;

%% Convolve Ambisonic with HRTF
for f_idx = 1:length(h5_files) % for all HOA RIRs
    h5_file = fullfile(h5_files(f_idx).folder, h5_files(f_idx).name);
    [~, h5_name, ~] = fileparts(h5_file);  % without extension

    % Output directory
    out_dir = fullfile(h5_files(f_idx).folder, h5_name);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    % Load HOA IR
    try
        uuid_list = h5info(h5_file, '/spatial_ir').Datasets;
        if isempty(uuid_list)
            warning('No spatial_ir datasets in %s', h5_file);
            continue;
        end
        uuid_key = uuid_list(1).Name;
        hoa_ir = h5read(h5_file, ['/spatial_ir/', uuid_key]);
    catch err
        warning('Failed to load IR from %s: %s', h5_file, err.message);
        continue;
    end

    % Convert: HOA(81) -> Spherical Domain(96)
    drirs = Y * hoa_ir.'; 
    
    % Upsample with SARITA: Spherical Domain(96) -> Target Domain (432)
    [drirs_upsampled, DRTFs_ups] = Sarita_upsampling( ...
        drirs, ...
        dirs, ...
        target_dirs, ...
        0.0875, ...
        'fs', fs_hoa, ...
        'f_transit', 700);
    
    M = size(drirs_upsampled, 1);
    L_hrtf = size(hrtf.Numerator, 3);
    L_resampled = ceil(L_hrtf * fs_hoa / fs_hrtf);

    h_L = zeros(M, L_resampled);
    h_R = zeros(M, L_resampled);

    for i = 1:M
        az1 = target_dirs(i,1); colat1 = target_dirs(i,2);
        az2 = dirs_hrtf(:,1); colat2 = dirs_hrtf(:,2);

        az_deg_target = rad2deg(az1);
        el_deg_target = rad2deg(pi/2 - colat1);  % elevation
        
        az_deg_cipic = hrtf.SourcePosition(:,1);
        el_deg_cipic = hrtf.SourcePosition(:,2);
        
        tol = 1e-3;
        idx = find(abs(az_deg_cipic - az_deg_target) < tol & abs(el_deg_cipic - el_deg_target) < tol, 1);
        
        if isempty(idx)
            % fallback to angular distance
            disp('idx of proper HRTF direction does not exist, using nearest neighbor')
            az2 = deg2rad(hrtf.SourcePosition(:,1));
            el2 = deg2rad(hrtf.SourcePosition(:,2));
            colat2 = pi/2 - el2;
            dist = acos(sin(colat1).*sin(colat2) + cos(colat1).*cos(colat2).*cos(az1 - az2));
            [~, idx] = min(dist);
        end
        h_L_44k = squeeze(hrtf.Numerator(idx, 1, :));
        h_R_44k = squeeze(hrtf.Numerator(idx, 2, :));
        h_L(i,:) = resample(h_L_44k, fs_hoa, fs_hrtf);
        h_R(i,:) = resample(h_R_44k, fs_hoa, fs_hrtf);
    end

    % Convolve with HRTF
    T = size(drirs_upsampled, 2);
    L = size(h_L, 2);
    N = T + L - 1;
    brirs = zeros(M, 2, N);

    for i = 1:M
        % Convolve
        brir_L = conv(drirs_upsampled(i,:), h_L(i,:), 'full');
        brir_R = conv(drirs_upsampled(i,:), h_R(i,:), 'full');
        brir = [brir_L; brir_R];
    
        % Normalize
        max_val = max(abs(brir), [], 'all');
        if max_val > 1
            brir = brir / max_val;
        end
    
        % File naming
        az_deg = rad2deg(target_dirs(i,1));
        el_deg = rad2deg(pi/2 - target_dirs(i,2));
        filename = sprintf('az%03d_el%+03d.wav', round(az_deg), round(el_deg));
        fullpath = fullfile(out_dir, filename);
    
        % Save
        audiowrite(fullpath, brir.', fs_hoa);  % [2 x N] → transpose for [N x 2]
    end
    
    fprintf('Finished processing: %s → %d BRIRs saved to %s\n', ...
        h5_name, M, out_dir);
end
