%% ACCESS THE RAW DATA FOLDER
raw_folder = fullfile('..', 'raw_data');
save_folder = '.';
all_data = [];

%% Stage mapping
stage_map = containers.Map(...
    {'4', '3', '2','1', '0', '5'}, ...
    {'N1', 'N2', 'N3','N3', 'REM', 'Wake'});

%% Filter parameters
low_cutoff = 0.5;
high_cutoff = 24.9;

%% Loop through all 10 subjects
for i = 1:10
    name = sprintf('R%d', i);
    edfFilename = fullfile(raw_folder, [name '.edf']);
    xmlFilename = fullfile(raw_folder, [name '.xml']);

    try
        [hdr, record] = edfread(edfFilename);
        [events, stages, epochLength, annotation] = readXML(xmlFilename);

    catch ME
        warning('Error reading files for %s: %s', name, ME.message);
        continue;
    end

    fs = 50;  % Hz
    ch_left = find(strcmp(hdr.label, 'EOGL'));
    ch_right = find(strcmp(hdr.label, 'EOGR'));

    if isempty(ch_left) || isempty(ch_right)
        warning('%s: EOGL or EOGR channel not found.', name);
        continue;
    end

    % Raw signals
    EOGL_raw = record(ch_left, :);
    EOGR_raw = record(ch_right, :);

    % Remove trailing zeros (zero-padding)
    last_valid_EOGL = find(EOGL_raw ~= 0, 1, 'last');
    last_valid_EOGR = find(EOGR_raw ~= 0, 1, 'last');
    last_valid = min(last_valid_EOGL, last_valid_EOGR);
    EOGL_raw = EOGL_raw(1:last_valid);
    EOGR_raw = EOGR_raw(1:last_valid);

    % Bandpass filter design
    %% Method 4: Eliptic (5-24.9 Hz) 
    [b, a] = ellip(4, 0.1, 60, [low_cutoff, high_cutoff]/(fs/2), 'bandpass');
    EOGL_filt = filtfilt(b, a, double(EOGL_raw));
    EOGR_filt = filtfilt(b, a, double(EOGR_raw));

    % Combine channels AFTER filtering
    combined_filt = (EOGL_filt + EOGR_filt) / 2;

    samples_per_epoch = fs * epochLength;  % = 1500
    num_epochs = floor(length(combined_filt) / samples_per_epoch) - 1 ;
    dominantStages = getDominantEpochStage(stages, fs, 30);

    subj_data = struct();
    subj_data.name = name;
    subj_data.epochs = cell(num_epochs, 1);

    for e = 1:num_epochs
        idx_start = (e - 1) * samples_per_epoch + 1;
        idx_end = e * samples_per_epoch;

        signal = combined_filt(idx_start:idx_end);

        if e <= length(dominantStages)
            label_raw = dominantStages(e);
            if isKey(stage_map, num2str(label_raw))
                label = stage_map(num2str(label_raw));
            end
        end

        subj_data.epochs{e} = struct( ...
            'signal', signal, ...
            'stage', label ...
        );
    end

    all_data = [all_data; subj_data]; %#ok<AGROW>
    fprintf('Processed subject: %s\n', name);
end


%% Save processed data
save(fullfile(save_folder, 'all_subjects_EOG_epochs.mat'), 'all_data', '-v7.3');
fprintf('Saved all processed data to: %s\n', fullfile(save_folder, 'all_subjects_EOG_epochs.mat'));

%% Create a CSV with Patient Number, Epoch, and Stage
% Initialize a cell array to store the data for the CSV
csvData = {};

% Loop through all subjects
for i = 1:length(all_data)
    patientID = i;  % Patient number (or use a unique identifier from the subject name)
    epochs = all_data(i).epochs;
    
    % Loop through all epochs for the current subject
    for e = 1:length(epochs)
        stage = epochs{e}.stage;  % Get the stage for the current epoch
        % Append patient ID, epoch number, and stage to the CSV data
        csvData = [csvData; {patientID, e, stage}];
    end
end

% Convert the cell array to a table for easy export
csvTable = cell2table(csvData, 'VariableNames', {'PatientID', 'Epoch', 'Stage'});

% Save the table as a CSV file
writetable(csvTable, 'Patient_Epochs_Stages.csv');
disp('CSV file created: Patient_Epochs_Stages.csv');


%% Summary by subject
fprintf('\nSummary of epochs per subject:\n');
for i = 1:length(all_data)
    numEpochs = length(all_data(i).epochs);
    fprintf('%s: %d epochs\n', all_data(i).name, numEpochs);
end
