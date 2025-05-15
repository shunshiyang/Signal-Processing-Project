% --- generate_EOG_features_table.m ---
% Load the dataset
data = load(fullfile('..', 'processed_data', 'all_subjects_EOG_epochs.mat'));
all_data = data.all_data;

% Sampling frequency and epoch duration
fs = 50;  % Hz
epochLength = 30;  % seconds

% Initialize final results table
finalTable = table();

% Loop over all patients
for p = 1:length(all_data)
    patient = all_data(p);
    patientID = p;
    epochs = patient.epochs;
    numEpochs = length(epochs);
   
    for e = 1:numEpochs
        epochData = epochs{e};

        if isfield(epochData, 'signal') && isfield(epochData, 'stage')
            signal = epochData.signal;
            stage = epochData.stage;
            
            % Use external function for feature extraction
            featureRow = extractEOGFeatures(signal, fs, epochLength);

            % Convert featureRow struct to table
            featureTable = struct2table(featureRow);

            % Add metadata as new columns
            featureTable.PatientID = patientID;
            featureTable.Epoch = e;
            featureTable.Stage = string(stage); 

            % Append to finalTable
            finalTable = [finalTable; featureTable];
        else
            warning('Missing EOGL/EOGR/stage for patient %d epoch %d. Skipping.', p, e);
        end
    end
    fprintf('Patient %d: Features extracted\n', patientID);
end

% Reorder for clarity
finalTable = movevars(finalTable, {'PatientID', 'Epoch', 'Stage'}, 'Before', 1);

% Save table
save('EOG_Features_Table.mat', 'finalTable');

% Load the table
data = load('EOG_Features_Table.mat');

% Extract the table (assuming the first variable is the table)
varNames = fieldnames(data);
T = data.(varNames{1});  % Adjust if the table variable has a known name

% Identify metadata columns
metaVars = {'PatientID', 'Epoch', 'Stage'};

% Get all feature names (excluding metadata)
featureNames = setdiff(T.Properties.VariableNames, metaVars);

% Define the subset of features you want to KEEP
featuresToKeep = {'Mean','Std','Variance','RMS','ZeroCrossings','SlopeSignChanges','WaveformLength','IEMG',...
                  'MAV','SSI','Kurtosis','Skewness','TotalPower','DeltaPower','ThetaPower','AlphaPower',...
                  'BetaPower','AlphaThetaRatio','BetaAlphaRatio','SpectralEntropy','SEF95','SEF50','SpectralCentroid',...
                  'SpectralFlatness','SpectralRollOff','HjorthActivity','HjorthMobility','HjorthComplexity','MovementDensity','BlinkRate','SEM_Rate','REM_Rate'};  

% Combine with metadata
finalVars = [metaVars, featuresToKeep];

% Create the new table
T_new = T(:, finalVars);

% Save as new .mat file and .csv
save('Reduced_EOG_Features_Table.mat', 'T_new');
writetable(T_new, 'Reduced_EOG_Features_Table.csv');

fprintf('\nReduced table saved with %d rows × %d columns\n', size(T_new,1), size(T_new,2));


% Load the reduced features table
load('Reduced_EOG_Features_Table.mat');

% Count number of epochs per sleep stage
[stage_counts, stage_names] = groupcounts(T_new.Stage);
total_epochs = sum(stage_counts);
percentages = 100 * stage_counts / total_epochs;

% Create figure with proper sizing
figure('Color', 'white', 'Position', [100 100 800 600]);
h = bar(stage_counts, 'FaceColor', [0.2 0.4 0.6], 'EdgeColor', 'none');

% Set y-axis limits with 25% padding
y_max = max(stage_counts) * 1.25;
ylim([0 y_max]);

% Add TOTAL COUNTS above bars (bold black)
for i = 1:length(stage_counts)
    text(i, stage_counts(i) + 0.03*y_max, ... % Positioned just above bar
        sprintf('%d', stage_counts(i)), ... % Just the number
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold', ...
        'FontSize', 11, ...
        'Color', 'k'); % Black color
end

% Add PERCENTAGES inside bars (bold white)
for i = 1:length(stage_counts)
    text(i, stage_counts(i)/2, ...
        sprintf('%.1f%%', percentages(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', 'white', ...
        'FontWeight', 'bold', ...
        'FontSize', 11);
end

% Formatting
title('Sleep Stage Distribution', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Sleep Stage', 'FontSize', 12);
ylabel('Number of Epochs', 'FontSize', 12);
xticklabels(stage_names);
set(gca, 'FontSize', 11, 'GridAlpha', 0.3, 'Layer', 'top');





