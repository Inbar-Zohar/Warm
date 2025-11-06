function [meas, prm, flg] = extract_raw_data_from_MFLI_folder2(folderPath, n_files, verbose)
% Handle optional inputs
if ~exist('verbose', 'var') || isempty(verbose); verbose = true; end
if ~exist('n_files', 'var'); n_files = []; end

% Handle empty folderPath
if nargin < 1 || isempty(folderPath)
    folderPath = uigetdir();
    if folderPath == 0
        error('Folder selection cancelled.');
    end
end

% For sim:

%%

if ~exist('SR3') 
    instrreset; 
    SR3 = SR86X('192.168.2.19'); 
    fopen(SR3.Ins); 
end

DorB = struct;
DorB.Xe129_LIA = GenericDevice(SR3);

Y_raw = SR3.ShowY;
Y = get_act_LIA_signal(Y_raw, 2, DorB.Xe129_LIA.apply('readConfiguration'));


%%
% Modular definition
isotopes = {'Xe_131', 'Xe_129'}; % Can be extended to any number of isotopes
demods_idxs = [1, 2]; % Indices for each isotope
pid_idxs = [1, 2]; % PID indices for each isotope
sample_fields = {'x', 'y', 'frequency'};
output_fields = {'Xraw', 'Yraw', 'drive_freq'};

% Validate input lengths
if length(demods_idxs) ~= length(isotopes) || length(pid_idxs) ~= length(isotopes)
    error('Length of isotopes, demods_idxs, and pid_idxs must match.');
end

% Create scalar parameter mapping
scalarMap = {};
for i = 1:length(isotopes)
    iso = isotopes{i};
    iso_field = strrep(iso, '_', '');  % Remove underscore for field name
    demod_idx = demods_idxs(i);
    pid_idx = pid_idxs(i);
    
    % LIA configuration
    scalarMap{end+1, 1} = ['WP.' iso_field '.LIAconfig.time_const'];
    scalarMap{end, 2} = @(data) data.dev32531.demods(demod_idx).timeconstant.value;
    
    scalarMap{end+1, 1} = ['WP.' iso_field '.LIAconfig.filter_order'];
    scalarMap{end, 2} = @(data) data.dev32531.demods(demod_idx).order.value;
    
    % PID configuration (conditional)
    scalarMap{end+1, 1} = ['WP.' iso_field '.PIDconfig.kP'];
    scalarMap{end, 2} = @(data) 2*pi*data.dev32531.pids(pid_idx).p.value*180/pi;
    scalarMap{end, 3} = @(data) isfield(data.dev32531, 'pids') && length(data.dev32531.pids) >= pid_idx;
    
    scalarMap{end+1, 1} = ['WP.' iso_field '.PIDconfig.kI'];
    scalarMap{end, 2} = @(data) 2*pi*data.dev32531.pids(pid_idx).i.value*180/pi;
    scalarMap{end, 3} = @(data) isfield(data.dev32531, 'pids') && length(data.dev32531.pids) >= pid_idx;
end

% Create extractMap as a cell array [fieldName, accessorFunction]
extractMap = {};
for i = 1:length(isotopes)
    iso = isotopes{i};
    idx = demods_idxs(i);
    for j = 1:length(sample_fields)
        sfield = sample_fields{j};
        ofield = output_fields{j};
        key = [iso '_' ofield];
        extractMap{end+1, 1} = key;
        extractMap{end, 2} = @(data) data.dev32531.demods(idx).sample.(sfield);
    end
    % Add time_clock for each isotope
    extractMap{end+1, 1} = [iso '_time_clock'];
    extractMap{end, 2} = @(data) data.dev32531.demods(idx).sample.timestamp;
end

% Get list of .mat files (no sorting, assuming natural order is correct)
files = dir(fullfile(folderPath, '*.mat'));
if isempty(files)
    error('No .mat files found in the specified folder.');
end
if ~isempty(n_files) && n_files < length(files); files = files(1:n_files); end

% Use modification time of the first file
firstFileTime = datetime(files(1).date);

% Initialize meas struct with DateTime_var and DateNumber
meas = struct();
meas.DateTime_var = firstFileTime;
meas.DateNumber = datenum(firstFileTime);

% Initialize prm struct
prm = struct();
prm.WP = struct();

% Initialize extracted fields
for k = 1:size(extractMap, 1)
    meas.(extractMap{k, 1}) = [];
end

% Loop through each file, check consistency within each isotope, and concatenate
record_lengths = zeros(length(files), length(isotopes)); % Store lengths for printout
for i = 1:length(files)
    filePath = fullfile(folderPath, files(i).name);
    data = load(filePath);
    
    % Extract scalar parameters from first file only
    if i == 1
        for k = 1:size(scalarMap, 1)
            fieldPath = scalarMap{k, 1};
            extractor = scalarMap{k, 2};
            
            % Check condition if exists (for PID parameters)
            if size(scalarMap, 2) >= 3 && ~isempty(scalarMap{k, 3})
                condition = scalarMap{k, 3};
                if ~condition(data)
                    continue;  % Skip if condition not met
                end
            end
            
            % Extract value
            value = extractor(data);
            
            % Assign value to prm struct
            eval(['prm.' fieldPath ' = ' num2str(value) ';']);
        end
    end
    
    % Check consistency of record lengths within each isotope
    for j = 1:length(isotopes)
        iso = isotopes{j};
        iso_fields = extractMap(startsWith(extractMap(:, 1), iso), :);
        lengths = zeros(1, size(iso_fields, 1));
        for k = 1:size(iso_fields, 1)
            val = iso_fields{k, 2}(data);
            lengths(k) = length(val(:));
        end
        if any(lengths ~= lengths(1))
            error('Inconsistent record lengths within isotope %s in file %s. Found lengths: %s', ...
                  iso, files(i).name, mat2str(lengths));
        end
        time_clock_idx = strcmp(iso_fields(:, 1), [iso '_time_clock']);
        record_lengths(i, j) = lengths(time_clock_idx);
    end
    
    % Extract and concatenate
    for k = 1:size(extractMap, 1)
        fieldName = extractMap{k, 1};
        val = extractMap{k, 2}(data);
        meas.(fieldName) = [meas.(fieldName); val(:)];
    end
end

% Print record lengths in one row per file
if verbose
    for i = 1:length(files)
        fprintf('File %s:', files(i).name);
        for j = 1:length(isotopes)
            fprintf(' %s time_clock length = %d', isotopes{j}, record_lengths(i, j));
            if j < length(isotopes)
                fprintf(',');
            end
        end
        fprintf('\n');
    end
end

% Align time_clock vectors to the isotope with the earliest start time
time_clocks = cellfun(@(iso) meas.([iso '_time_clock']), isotopes, 'UniformOutput', false);
first_times = cellfun(@(tc) double(tc(1)), time_clocks);
[~, ref_idx] = max(first_times); % Isotope with latest start time
ref_iso = isotopes{ref_idx};
ref_time_clock = meas.([ref_iso '_time_clock']);

% Determine the common time range (latest start to earliest end)
last_times = cellfun(@(tc) double(tc(end)), time_clocks);
latest_end = min(last_times); % Earliest end time 

% Trim each isotope to align with the reference's start and the common end
for i = 1:length(isotopes)
    iso = isotopes{i};
    time_clock = meas.([iso '_time_clock']);
    
    % Skip if this is the reference isotope
    if i == ref_idx
        continue;
    end
    
    % Check and trim start
    if ~isempty(time_clock) && double(time_clock(1)) ~= double(ref_time_clock(1))
        start_idx = find(double(time_clock) == double(ref_time_clock(1)), 1);
        if isempty(start_idx)
            error('Cannot align time_clock for %s: No matching start value with %s.', iso, ref_iso);
        end
        meas.([iso '_Xraw']) = meas.([iso '_Xraw'])(start_idx:end);
        meas.([iso '_Yraw']) = meas.([iso '_Yraw'])(start_idx:end);
        meas.([iso '_drive_freq']) = meas.([iso '_drive_freq'])(start_idx:end);
        meas.([iso '_time_clock']) = meas.([iso '_time_clock'])(start_idx:end);
    end
    
    % Trim end to match the latest common end time
    time_clock = meas.([iso '_time_clock']);
    if double(time_clock(end)) > latest_end
        end_idx = find(double(time_clock) <= latest_end, 1, 'last');
        if isempty(end_idx)
            error('Cannot align time_clock for %s: No matching end value within %s range.', iso, ref_iso);
        end
        meas.([iso '_Xraw']) = meas.([iso '_Xraw'])(1:end_idx);
        meas.([iso '_Yraw']) = meas.([iso '_Yraw'])(1:end_idx);
        meas.([iso '_drive_freq']) = meas.([iso '_drive_freq'])(1:end_idx);
        meas.([iso '_time_clock']) = meas.([iso '_time_clock'])(1:end_idx);
    end
end

% Trim reference isotope to match the latest common end time
if double(ref_time_clock(end)) > latest_end
    end_idx = find(double(ref_time_clock) <= latest_end, 1, 'last');
    if isempty(end_idx)
        error('Cannot align reference time_clock for %s: No matching end value.', ref_iso);
    end
    meas.([ref_iso '_Xraw']) = meas.([ref_iso '_Xraw'])(1:end_idx);
    meas.([ref_iso '_Yraw']) = meas.([ref_iso '_Yraw'])(1:end_idx);
    meas.([ref_iso '_drive_freq']) = meas.([ref_iso '_drive_freq'])(1:end_idx);
    meas.([ref_iso '_time_clock']) = meas.([ref_iso '_time_clock'])(1:end_idx);
end

% Verify all time_clock vectors are identical
ref_time_clock = meas.([ref_iso '_time_clock']);
for i = 1:length(isotopes)
    iso = isotopes{i};
    if ~isequal(double(meas.([iso '_time_clock'])), double(ref_time_clock))
        error('time_clock vectors for %s and %s are not identical after alignment.', iso, ref_iso);
    end
end

% Store only one time_clock and remove individual isotope time_clocks
meas.time_clock = meas.([ref_iso '_time_clock']);
for i = 1:length(isotopes)
    meas = rmfield(meas, [isotopes{i} '_time_clock']);
end

% Add calculated fields for each isotope
for i = 1:length(isotopes)
    iso = isotopes{i};
    prefix = [iso '_'];
    X = meas.([prefix 'Xraw']);
    Y = meas.([prefix 'Yraw']);
    
    R = sqrt(X.^2 + Y.^2);
    Theta = atan2(Y, X) * 180 / pi;  % In degrees
    meanTheta = mean(Theta);
    dY = R .* sind(Theta - meanTheta);
    
    drive_freq = meas.([prefix 'drive_freq']);
    Omega_deg_hr = (drive_freq - 0*mean(drive_freq)) * 360 * 3600;
    iso_num = iso(end-2:end);
    
    meas.([prefix 'R']) = R;
    meas.([prefix 'Theta']) = Theta;
    meas.([prefix 'Y']) = dY;
    meas.(['Omega_' iso_num '_deg_hr']) = Omega_deg_hr;
end

% Add Omega and Bz analysis (only for Xe_131 and Xe_129)
if ~all(ismember({'Xe_131', 'Xe_129'}, isotopes))
    error('Omega and Bz analysis requires both Xe_131 and Xe_129 to be present.');
end
ga = (-7.44e7) / (2*pi) / (1e4);  % [Hz/G] Xe129 gyromagnetic ratio
gb = (+2.2056e7) / (2*pi) / (1e4);  % [Hz/G] Xe131 gyromagnetic ratio
R = abs(ga / gb);

Omega_deg_hr = meas.Omega_129_deg_hr / (1+R) - meas.Omega_131_deg_hr * R / (1+R);
meas.Omega_deg_hr_mean = mean(Omega_deg_hr);
meas.Omega_deg_hr = Omega_deg_hr - mean(Omega_deg_hr); 

Bz = (meas.Omega_129_deg_hr + meas.Omega_131_deg_hr) / (abs(ga) + gb) / (360 * 3600);  % [G]
meas.Bz_mean = mean(Bz);
meas.Bz = Bz - mean(Bz);

% Add time_s and fs
time_s = double(meas.time_clock) / 60e6;
meas.time_s = time_s - time_s(1);
fs = 1 / diff(time_s(1:2));
meas.fs = fs;

flg.plot_individual_xenons = 1;
