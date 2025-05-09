%% Clear all variables and prompt
% EEG Signal processing pipeline
% Preprocessing and feature extraction for EEG only

%% Clear prompt and all variable
clc, clf;
close all;
clear all;

%% Setup
params = configureParams();
%% Create structure to store processed data
procEdfData = struct();

%% Process each patient's data
patientIDs = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10'};
% patientIDs = {'R1', 'R2'};
for k = 1:length(patientIDs)
    patient = patientIDs{k};
    
    % Read EDF and XML files
    edfFilename = [patient '.edf'];
    xmlFilename = [patient '.xml'];
    
    fprintf('Processing patient %s\n', patient);
    
    try

        [hdr, record, stages, events, epochLength] = loadData(k, params.dataFolder);
        
        % Find EEG channel indices
        eegChannels = find(contains(lower(hdr.label), 'eeg'));
        
        if isempty(eegChannels)
            warning('No EEG channels found for patient %s', patient);
            continue;
        end
        
        % Process each EEG channel
        for ch = 1:length(eegChannels)
            channelIdx = eegChannels(ch);
            channelLabel = hdr.label{channelIdx};
            
            fprintf('Processing EEG channel: %s\n', channelLabel);
            
            % Get sampling rate using the correct field (samples instead of frequency)
            if isfield(hdr, 'samples')
                Fs = hdr.samples(channelIdx);
                fprintf('Using sampling rate: %d Hz from hdr.samples\n', Fs);
            else
                warning('Cannot find sampling rate information. Using default value of 256 Hz');
                Fs = 256;
            end

            sleepStages_perEpoch = processSleepStages(stages, Fs, 30);   
            
            % Get EEG signal
            eegSignal = record(channelIdx, :);
            
            % Divide signal into 30-second epochs
            epochSamples = Fs * 30; % 30 seconds
            numEpochs = floor(length(stages) / 30);
            
            fprintf('Signal length: %d samples, creating %d epochs of %d samples each\n', ...
                    length(stages), numEpochs, 30);
            
            % Initialize arrays for storing processed data
            rawEpochs = zeros(numEpochs, epochSamples);
            ppEpochs = zeros(numEpochs, epochSamples);
            features = zeros(numEpochs, 12); % 12 features per epoch
            
            % Process each epoch
            for e = 1:numEpochs
                % Extract epoch
                startIdx = (e-1) * epochSamples + 1;
                endIdx = e * epochSamples;
                
                if endIdx <= length(eegSignal)
                    epochData = eegSignal(startIdx:endIdx);
                    rawEpochs(e, :) = epochData;
                    
                    % Preprocess epoch
                    ppEpochs(e, :) = preProEEG(epochData, Fs);
                    
                    % Extract features
                    if e < numEpochs
                        nextEpochData = eegSignal(endIdx+1:min(endIdx+epochSamples, length(eegSignal)));
                        
                        % Ensure both epochs have the same length
                        if length(nextEpochData) < epochSamples
                            % Pad with zeros if needed
                            nextEpochData = [nextEpochData, zeros(1, epochSamples - length(nextEpochData))];
                        end
                        
                        features(e, :) = featExtEEG([ppEpochs(e, :); nextEpochData], Fs);
                    else
                        % For the last epoch, duplicate it
                        features(e, :) = featExtEEG([ppEpochs(e, :); ppEpochs(e, :)], Fs);
                    end
                end
            end
            
            % Ensure stages has the right length
            if length(stages) >= numEpochs
                epochStages = stages(1:numEpochs);
            else
                warning('Number of stages (%d) is less than number of epochs (%d). Padding with zeros.', ...
                        length(stages), numEpochs);
                epochStages = [stages, zeros(1, numEpochs - length(stages))];
            end
            
            % Store processed data in structure
            procEdfData.(patient).(channelLabel).raw = rawEpochs;
            procEdfData.(patient).(channelLabel).preprocessed = ppEpochs;
            procEdfData.(patient).(channelLabel).features = features;
            procEdfData.(patient).(channelLabel).fs = Fs;
            procEdfData.(patient).(channelLabel).stages = epochStages;
            fprintf('features_Matrix: %d row x %d colum\n', size(features, 1), size(features, 2));
            fprintf('sleepStages_perEpoch: %d\n', length(sleepStages_perEpoch));
            featureTable = createFeatureTable(features, sleepStages_perEpoch);
            saveFeatures(featureTable, k, channelLabel, params.FeatureFolder);
            
            fprintf('Completed processing for channel %s\n', channelLabel);
        end
        
    catch ME
        warning('Error processing patient %s: %s\n%s', patient, ME.message, getReport(ME));
        continue;
    end
end

%% Save processed data
save('eegProcessedData.mat', 'procEdfData', '-v7.3');
fprintf('All processed data saved to eegProcessedData.mat\n');

%% EEG preprocessing function
function out = preProEEG(signal, Fs)
    % Apply bandpass filter for EEG preprocessing
    try
        out = bandpass(signal, [0.5 30], Fs);
    catch ME
        % In case bandpass fails (e.g., Signal Processing Toolbox not available)
        warning('Bandpass filter failed: %s. Using simple filtering instead.', ME.message);
        
        % Simple high-pass filtering (remove DC component)
        signal = signal - mean(signal);
        
        % Simple low-pass filtering (moving average)
        windowSize = round(Fs / 30); % Roughly corresponds to 30 Hz cutoff
        if windowSize > 1
            b = ones(1, windowSize) / windowSize;
            out = filter(b, 1, signal);
        else
            out = signal;
        end
    end
end

%% EEG feature extraction function
function out = featExtEEG(signal, Fs)
    try
        % Extract spectral features
        % Calculate power in different frequency bands using bandpower function
        
        % Delta (0.5-4 Hz)
        deltaFeature = bandpower(signal(1,:), Fs, [0.5 4]);
        
        % Theta (4-8 Hz)
        thetaFeature = bandpower(signal(1,:), Fs, [4 8]);
        
        % Alpha (8-12 Hz)
        alphaFeature = bandpower(signal(1,:), Fs, [8 12]);
        
        % Beta (12-30 Hz)
        betaFeature = bandpower(signal(1,:), Fs, [12 30]);
        
    catch ME
        % If bandpower fails, calculate power manually
        warning('bandpower function failed: %s. Using manual power calculation.', ME.message);
        
        % Calculate power manually using FFT
        L = length(signal(1,:));
        NFFT = 2^nextpow2(L);
        Y = fft(signal(1,:), NFFT)/L;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        
        % Find indices for each frequency band
        deltaIdx = (f >= 0.5) & (f <= 4);
        thetaIdx = (f > 4) & (f <= 8);
        alphaIdx = (f > 8) & (f <= 12);
        betaIdx = (f > 12) & (f <= 30);
        
        % Calculate power in each band
        P = abs(Y(1:NFFT/2+1)).^2;
        deltaFeature = sum(P(deltaIdx));
        thetaFeature = sum(P(thetaIdx));
        alphaFeature = sum(P(alphaIdx));
        betaFeature = sum(P(betaIdx));
    end
    
    % Calculate temporal features
    meanFeature = mean(signal(1,:));
    stdFeature = std(signal(1,:));
    varianceFeature = var(signal(1,:));
    
    try
        kurtosisFeature = kurtosis(signal(1,:));
        skewnessFeature = skewness(signal(1,:));
    catch ME
        % If kurtosis/skewness functions are not available
        warning('kurtosis/skewness functions failed: %s. Using manual calculation.', ME.message);
        
        % Manual kurtosis calculation
        x = signal(1,:) - mean(signal(1,:));
        n = length(x);
        m4 = sum(x.^4)/n;
        m2 = sum(x.^2)/n;
        kurtosisFeature = m4/(m2^2) - 3;
        
        % Manual skewness calculation
        m3 = sum(x.^3)/n;
        skewnessFeature = m3/(m2^(3/2));
    end
    
    try
        % Calculate spectral edge frequency
        [pxx, f] = pwelch(signal(1,:), [], [], [], Fs);
        totalPower = sum(pxx);
        cumulativePower = cumsum(pxx) / totalPower;
        sef95 = f(find(cumulativePower >= 0.95, 1, 'first'));
    catch ME
        % If pwelch fails
        warning('pwelch function failed: %s. Using simple spectral analysis.', ME.message);
        
        % Simple spectral edge frequency calculation
        L = length(signal(1,:));
        NFFT = 2^nextpow2(L);
        Y = fft(signal(1,:), NFFT)/L;
        f = Fs/2*linspace(0,1,NFFT/2+1);
        P = abs(Y(1:NFFT/2+1)).^2;
        
        cumP = cumsum(P) / sum(P);
        sef95 = f(find(cumP >= 0.95, 1, 'first'));
    end
    
    % Visual features (K-complex, Sleep spindles)
    kComplexCount = countKComplex(signal(1,:), Fs);
    sleepSpindleCount = countSleepSpindles(signal(1,:), Fs);
    
    % Combine all features
    out = [deltaFeature, thetaFeature, alphaFeature, betaFeature, ...
           meanFeature, stdFeature, varianceFeature, kurtosisFeature, ...
           skewnessFeature, sef95, kComplexCount, sleepSpindleCount];
end

%% K-complex detection function
function count = countKComplex(signal, Fs)
    try
        % Filter signal to enhance K-complexes
        filteredSignal = bandpass(signal, [0.5 4], Fs);
    catch ME
        % Simple filter if bandpass fails
        filteredSignal = signal - mean(signal);
        windowSize = round(Fs / 4); % Roughly corresponds to 4 Hz cutoff
        if windowSize > 1
            b = ones(1, windowSize) / windowSize;
            filteredSignal = filter(b, 1, filteredSignal);
        end
    end
    
    % Find potential K-complexes based on amplitude
    threshold = 2 * std(filteredSignal);
    
    try
        [~, locs] = findpeaks(-filteredSignal, 'MinPeakHeight', threshold);
    catch ME
        % Manual peak detection if findpeaks fails
        warning('findpeaks function failed: %s. Using manual peak detection.', ME.message);
        
        % Simple peak detection
        isPeak = zeros(size(filteredSignal));
        for i = 2:length(filteredSignal)-1
            if -filteredSignal(i) > threshold && ...
               -filteredSignal(i) > -filteredSignal(i-1) && ...
               -filteredSignal(i) > -filteredSignal(i+1)
                isPeak(i) = 1;
            end
        end
        locs = find(isPeak);
    end
    
    % Count potential K-complexes
    count = length(locs);
end

%% Sleep spindle detection function
function count = countSleepSpindles(signal, Fs)
    try
        % Filter signal to isolate spindle frequency band
        filteredSignal = bandpass(signal, [12 16], Fs);
    catch ME
        % Simple filter if bandpass fails
        filteredSignal = signal - mean(signal);
        b = fir1(round(Fs/4), [12 16]/(Fs/2));
        if ~isempty(b)
            filteredSignal = filter(b, 1, filteredSignal);
        end
    end
    
    % Calculate envelope
    try
        envelope = abs(hilbert(filteredSignal));
    catch ME
        % Simple envelope calculation if hilbert fails
        warning('hilbert function failed: %s. Using manual envelope.', ME.message);
        envelope = sqrt(filteredSignal.^2);
    end
    
    % Find potential spindles based on amplitude
    threshold = 2 * std(envelope);
    
    try
        [~, locs] = findpeaks(envelope, 'MinPeakHeight', threshold);
    catch ME
        % Manual peak detection if findpeaks fails
        warning('findpeaks function failed: %s. Using manual peak detection.', ME.message);
        
        % Simple peak detection
        isPeak = zeros(size(envelope));
        for i = 2:length(envelope)-1
            if envelope(i) > threshold && ...
               envelope(i) > envelope(i-1) && ...
               envelope(i) > envelope(i+1)
                isPeak(i) = 1;
            end
        end
        locs = find(isPeak);
    end
    
    % Count potential spindles
    count = length(locs);
end
