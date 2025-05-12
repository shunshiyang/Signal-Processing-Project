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
            features = zeros(numEpochs, 32); % 32 features per epoch
            
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
            procEdfData.(patient).(channelLabel).Fs = Fs;
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
        delta = bandpower(signal(1,:), Fs, [0.5 4]);
        
        % Theta (4-8 Hz)
        theta = bandpower(signal(1,:), Fs, [4 8]);
        
        % Alpha (8-12 Hz)
        alpha = bandpower(signal(1,:), Fs, [8 12]);
        
        % Beta (12-30 Hz)
        beta = bandpower(signal(1,:), Fs, [12 30]);
        
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
        delta = sum(P(deltaIdx));
        theta = sum(P(thetaIdx));
        alpha = sum(P(alphaIdx));
        beta = sum(P(betaIdx));
    end
    %another frequency domain femeanatures
    alpha_power = calculate_alpha_wave_power(signal, Fs);
    alpha_ratio = calculate_alpha_wave_ratio(signal, Fs);
    wave_percentages = calculate_wave_percentages(signal, Fs);
    wave_delta=wave_percentages.delta;
    wave_theta=wave_percentages.theta;
    wave_alpha=wave_percentages.alpha;
    wave_beta=wave_percentages.beta;
    centroid = calcentroid(signal, Fs);
    flatness = calculate_spectral_flatness(signal, Fs);
    % Calculate temporal features
    meanFeature = mean(signal(1,:));
    stdFeature = std(signal(1,:));
    variance = var(signal(1,:));
    ssi = calculate_energy_ssi(signal);
    max_val = calculate_maximum(signal);
    max_peak = calculate_maximum_peak(signal);
    lrssv = calculate_log_root_sum_sequential_variation(signal);
    fd = calculate_first_difference(signal);
    sd = calculate_second_difference(signal);
    zcr = calculate_zero_crossing_rate(signal);
    
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
    %Check
    % disp(['delta: ', mat2str(size(delta))]);
    % disp(['theta: ', mat2str(size(theta))]);
    % disp(['alpha: ', mat2str(size(alpha))]);
    % disp(['beta: ', mat2str(size(beta))]);
    % disp(['alpha_power: ', mat2str(size(alpha_power))]);
    % disp(['alpha_ratio: ', mat2str(size(alpha_ratio))]);
    % disp(['wave_percentages.delta: ', mat2str(size(wave_percentages.delta))]);
    % disp(['wave_percentages.theta: ', mat2str(size(wave_percentages.theta))]);
    % disp(['wave_percentages.alpha: ', mat2str(size(wave_percentages.alpha))]);
    % disp(['wave_percentages.beta: ', mat2str(size(wave_percentages.beta))]);
    % disp(['wave_percentages.gamma: ', mat2str(size(wave_percentages.gamma))]);
    % disp(['wave_delta: ', mat2str(size(wave_delta))]);
    % disp(['wave_theta: ', mat2str(size(wave_theta))]);
    % disp(['wave_alpha: ', mat2str(size(wave_alpha))]);
    % disp(['wave_beta: ', mat2str(size(wave_beta))]);
    % disp(['centroid: ', mat2str(size(centroid))]);
    % disp(['flatness: ', mat2str(size(flatness))]);
    % disp(['sef95: ', mat2str(size(sef95))]);
    % disp(['meanFeature: ', mat2str(size(meanFeature))]);
    % disp(['stdFeature: ', mat2str(size(stdFeature))]);
    % disp(['variance: ', mat2str(size(variance))]);
    % disp(['kurtosisFeature: ', mat2str(size(kurtosisFeature))]);
    % disp(['skewnessFeature: ', mat2str(size(skewnessFeature))]);
    % disp(['ssi: ', mat2str(size(ssi))]);
    % disp(['max_val: ', mat2str(size(max_val))]);
    % disp(['max_peak: ', mat2str(size(max_peak))]);
    % disp(['lrssv: ', mat2str(size(lrssv))]);
    % disp(['fd: ', mat2str(size(fd))]);
    % disp(['sd: ', mat2str(size(sd))]);
    % disp(['zcr: ', mat2str(size(zcr))]);
    % disp(['kComplexCount: ', mat2str(size(kComplexCount))]);
    % disp(['sleepSpindleCount: ', mat2str(size(sleepSpindleCount))]);
    % Combine all features
    out = zeros(1, 32);
    out(1) = delta;
    out(2) = theta;
    out(3) = alpha;
    out(4) = beta;
    out(5) = alpha_power;
    out(6) = alpha_ratio;
    out(7) = wave_percentages.delta;
    out(8) = wave_percentages.theta;
    out(9) = wave_percentages.alpha;
    out(10) = wave_percentages.beta;
    out(11) = wave_percentages.gamma;
    out(12) = wave_delta;
    out(13) = wave_theta;
    out(14) = wave_alpha;
    out(15) = wave_beta;
    out(16) = centroid;
    out(17) = flatness;
    out(18) = sef95;
    out(19) = meanFeature;
    out(20) = stdFeature;
    out(21) = variance;
    out(22) = kurtosisFeature;
    out(23) = skewnessFeature;
    out(24) = ssi;
    out(25) = max_val;
    out(26) = max_peak;
    out(27) = lrssv;
    out(28) = fd;
    out(29) = sd;
    out(30) = zcr;
    out(31) = kComplexCount;
    out(32) = sleepSpindleCount;
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
%% Centroid function
function centroid = calcentroid(signal, Fs)
    % Ensure signal is a one-dimensional vector
    if size(signal, 1) > 1
        signal = signal(1, :);  % If multi-channel, use only the first channel
    else
        signal = signal(1, :);  % Ensure using the first row
    end
    
    % Check signal length
    if length(signal) < 8
        warning('Signal too short for pwelch (< 8 samples). Using FFT directly.');
        % For short signals, use simple FFT method
        N = length(signal);
        if N < 2
            centroid = 0;  % Signal too short, return 0
            return;
        end
        
        % Use simple FFT
        nfft = max(256, 2^nextpow2(N));  % Ensure sufficient frequency resolution
        Y = fft(signal, nfft);
        P = abs(Y/N).^2;
        P = P(1:floor(nfft/2)+1);  % Single-sided spectrum
        f = Fs/2 * linspace(0, 1, floor(nfft/2)+1);
        
        % Normalize power spectrum
        if sum(P) > 0
            P = P / sum(P);
        else
            centroid = 0;
            return;
        end
        
        % Calculate spectral centroid
        centroid = sum(f .* P);
    else
        % For signals long enough, use pwelch
        try
            % Use pwelch with explicitly specified window size
            winSize = min(128, floor(length(signal)/2));  % Ensure window not larger than half signal length
            [pxx, f] = pwelch(signal, hamming(winSize), [], [], Fs);
            
            % Normalize power spectrum
            if sum(pxx) > 0
                pxx = pxx / sum(pxx);
            else
                centroid = 0;
                return;
            end
            
            % Calculate spectral centroid
            centroid = sum(f .* pxx);
        catch ME
            warning('pwelch failed in calcentroid: %s. Setting to 0.', ME.message);
            centroid = 0;
        end
    end
end

%% flatness
function flatness = calculate_spectral_flatness(signal, Fs)
    % Ensure signal is a one-dimensional vector
    if size(signal, 1) > 1
        signal = signal(1, :);  % If multi-channel, use only the first channel
    else
        signal = signal(1, :);  % Ensure using the first row
    end
    
    % Check signal length
    if length(signal) < 8
        warning('Signal too short for pwelch (< 8 samples). Using FFT directly.');
        % For short signals, use simple FFT method
        N = length(signal);
        if N < 2
            flatness = 0;  % Signal too short, return 0
            return;
        end
        
        % Use simple FFT
        nfft = max(256, 2^nextpow2(N));  % Ensure sufficient frequency resolution
        Y = fft(signal, nfft);
        spectrum = abs(Y(2:floor(nfft/2)+1));  % Exclude DC component
        
        % Ensure no zero values (would cause geometric mean to be zero)
        spectrum = max(spectrum, eps);
        
        if length(spectrum) < 2
            flatness = 0;
            return;
        end
        
        % Calculate geometric mean
        geo_mean = exp(mean(log(spectrum)));
        
        % Calculate arithmetic mean
        arith_mean = mean(spectrum);
        
        % Calculate spectral flatness
        if arith_mean > 0
            flatness = geo_mean / arith_mean;
        else
            flatness = 0;
        end
    else
        % For signals long enough, use pwelch
        try
            % Use pwelch with explicitly specified window size
            winSize = min(128, floor(length(signal)/2));  % Ensure window not larger than half signal length
            [pxx, ~] = pwelch(signal, hamming(winSize), [], [], Fs);
            
            % Exclude possible zero values
            pxx = max(pxx, eps);
            
            % Calculate geometric mean
            geo_mean = exp(mean(log(pxx)));
            
            % Calculate arithmetic mean
            arith_mean = mean(pxx);
            
            % Calculate spectral flatness
            if arith_mean > 0
                flatness = geo_mean / arith_mean;
            else
                flatness = 0;
            end
        catch ME
            warning('pwelch failed in calculate_spectral_flatness: %s. Setting to 0.', ME.message);
            flatness = 0;
        end
    end
end

%% edge_freq
function edge_freq=calculate_edge_freq(signal,Fs)
    if nargin < 3
        percentage = 95; 
    end

    signal = signal(:);
    
    N = length(signal);
    nfft = 2^nextpow2(N);
    
    [pxx, f] = pwelch(signal(1,:), [], [], [], Fs);
    
    cum_sum = cumsum(pxx);

    threshold = cum_sum(end) * (percentage / 100);

    idx = find(cum_sum >= threshold, 1, 'first');

    if isempty(idx)
        edge_freq = f(end);
    else
        edge_freq = f(idx);
    end
end

%% ssI

function ssi = calculate_energy_ssi(signal)
    signal = signal(:);
    ssi = sum(signal.^2);
end

%% max_val

function max_val = calculate_maximum(signal)
    signal = signal(:);
    max_val = max(signal);
end

%% max_peak

function max_peak = calculate_maximum_peak(signal)
    signal = signal(:);
    max_peak = max(abs(signal));
end

%% lrssv

function lrssv = calculate_log_root_sum_sequential_variation(signal)
    signal = signal(:);
    difFs = diff(signal);
    sum_squares = sum(difFs.^2);
    
    if sum_squares > 0
        lrssv = log10(sqrt(sum_squares));
    else
        lrssv = 0;
    end
end

%% fd

function fd = calculate_first_difference(signal)
    signal = signal(:);
    difFs = abs(diff(signal));
    fd = mean(difFs);
end

%% sd

function sd = calculate_second_difference(signal)
    signal = signal(:);
    second_difFs = diff(diff(signal));
    sd = mean(abs(second_difFs));
end

%% zcr
function zcr = calculate_zero_crossing_rate(signal)
    signal = signal(:);
    N = length(signal);
    sign_changes = abs(diff(sign(signal)));
    zero_crossings = sum(sign_changes > 0);
    zcr = zero_crossings / (N - 1);
end

%% alpha_power
function alpha_power = calculate_alpha_wave_power(signal, Fs)
    % Ensure signal is a one-dimensional vector
    if size(signal, 1) > 1
        signal = signal(1, :);  % If multi-channel, use only the first channel
    else
        signal = signal(1, :);  % Ensure using the first row
    end
    
    % Define Alpha wave band range (8-12 Hz)
    alpha_low = 8;
    alpha_high = 12;
    
    % Check signal length
    if length(signal) < 8
        warning('Signal too short for pwelch (< 8 samples). Using FFT directly.');
        % For short signals, use simple FFT method
        N = length(signal);
        if N < 2
            alpha_power = 0;  % Signal too short, return 0
            return;
        end
        
        % Use simple FFT
        nfft = max(256, 2^nextpow2(N));  % Ensure sufficient frequency resolution
        Y = fft(signal, nfft);
        P = abs(Y/N).^2;
        P = P(1:floor(nfft/2)+1);  % Single-sided spectrum
        f = Fs/2 * linspace(0, 1, floor(nfft/2)+1);
        
        % Find indices for Alpha band frequencies
        alpha_idx = (f >= alpha_low) & (f <= alpha_high);
        
        % Calculate Alpha band power
        if any(alpha_idx)
            alpha_power = sum(P(alpha_idx));
        else
            alpha_power = 0;
        end
    else
        % For signals long enough, use pwelch
        try
            % Use pwelch with explicitly specified window size
            winSize = min(128, floor(length(signal)/2));  % Ensure window not larger than half signal length
            [pxx, f] = pwelch(signal, hamming(winSize), [], [], Fs);
            
            % Find indices for Alpha band frequencies
            alpha_idx = (f >= alpha_low) & (f <= alpha_high);
            
            % Calculate Alpha band power
            if any(alpha_idx)
                alpha_power = sum(pxx(alpha_idx));
            else
                alpha_power = 0;
            end
        catch ME
            warning('pwelch failed in calculate_alpha_wave_power: %s. Setting to 0.', ME.message);
            alpha_power = 0;
        end
    end
end

%% alpha_ratio

function alpha_ratio = calculate_alpha_wave_ratio(signal, Fs)
    % Ensure signal is a one-dimensional vector
    if size(signal, 1) > 1
        signal = signal(1, :);  % If multi-channel, use only the first channel
    else
        signal = signal(1, :);  % Ensure using the first row
    end
    
    % Define Alpha wave band range (8-12 Hz)
    alpha_low = 8;
    alpha_high = 12;
    
    % Check signal length
    if length(signal) < 8
        warning('Signal too short for pwelch (< 8 samples). Using FFT directly.');
        % For short signals, use simple FFT method
        N = length(signal);
        if N < 2
            alpha_ratio = 0;  % Signal too short, return 0
            return;
        end
        
        % Use simple FFT
        nfft = max(256, 2^nextpow2(N));  % Ensure sufficient frequency resolution
        Y = fft(signal, nfft);
        P = abs(Y/N).^2;
        P = P(1:floor(nfft/2)+1);  % Single-sided spectrum
        f = Fs/2 * linspace(0, 1, floor(nfft/2)+1);
        
        % Find indices for Alpha band frequencies
        alpha_idx = (f >= alpha_low) & (f <= alpha_high);
        
        % Calculate Alpha band power
        alpha_power = sum(P(alpha_idx));
        
        % Calculate total power
        total_power = sum(P);
        
        % Calculate Alpha wave ratio
        if total_power > 0
            alpha_ratio = alpha_power / total_power;
        else
            alpha_ratio = 0;
        end
    else
        % For signals long enough, use pwelch
        try
            % Use pwelch with explicitly specified window size
            winSize = min(128, floor(length(signal)/2));  % Ensure window not larger than half signal length
            [pxx, f] = pwelch(signal, hamming(winSize), [], [], Fs);
            
            % Find indices for Alpha band frequencies
            alpha_idx = (f >= alpha_low) & (f <= alpha_high);
            
            % Calculate Alpha band power
            alpha_power = sum(pxx(alpha_idx));
            
            % Calculate total power
            total_power = sum(pxx);
            
            % Calculate Alpha wave ratio
            if total_power > 0
                alpha_ratio = alpha_power / total_power;
            else
                alpha_ratio = 0;
            end
        catch ME
            warning('pwelch failed in calculate_alpha_wave_ratio: %s. Setting to 0.', ME.message);
            alpha_ratio = 0;
        end
    end
end

%% wave_percentages

function wave_percentages = calculate_wave_percentages(signal, Fs)
    % Ensure signal is a one-dimensional vector
    if size(signal, 1) > 1
        signal = signal(1, :);  % If multi-channel, use only the first channel
    else
        signal = signal(1, :);  % Ensure using the first row
    end
    
    % Define frequency band ranges
    delta_range = [0.5, 4];   % Delta: 0.5-4 Hz
    theta_range = [4, 8];     % Theta: 4-8 Hz
    alpha_range = [8, 12];    % Alpha: 8-12 Hz
    beta_range = [12, 30];    % Beta: 12-30 Hz
    gamma_range = [30, Fs/2]; % Gamma: 30+ Hz (up to Nyquist frequency)
    
    % Create structure to store percentages for each band
    wave_percentages = struct('delta', 0, 'theta', 0, 'alpha', 0, 'beta', 0, 'gamma', 0);
    
    % Check signal length
    if length(signal) < 8
        warning('Signal too short for pwelch (< 8 samples). Using FFT directly.');
        % For short signals, use simple FFT method
        N = length(signal);
        if N < 2
            return;  % Signal too short, return default structure
        end
        
        % Use simple FFT
        nfft = max(256, 2^nextpow2(N));  % Ensure sufficient frequency resolution
        Y = fft(signal, nfft);
        P = abs(Y/N).^2;
        P = P(1:floor(nfft/2)+1);  % Single-sided spectrum
        f = Fs/2 * linspace(0, 1, floor(nfft/2)+1);
        
        % Find indices for each frequency band
        delta_idx = (f >= delta_range(1)) & (f <= delta_range(2));
        theta_idx = (f > theta_range(1)) & (f <= theta_range(2));
        alpha_idx = (f > alpha_range(1)) & (f <= alpha_range(2));
        beta_idx = (f > beta_range(1)) & (f <= beta_range(2));
        gamma_idx = (f > gamma_range(1)) & (f <= gamma_range(2));
        
        % Calculate power in each band
        delta_power = sum(P(delta_idx));
        theta_power = sum(P(theta_idx));
        alpha_power = sum(P(alpha_idx));
        beta_power = sum(P(beta_idx));
        gamma_power = sum(P(gamma_idx));
        
        % Calculate total power
        total_power = delta_power + theta_power + alpha_power + beta_power + gamma_power;
        
        % Calculate percentages for each band
        if total_power > 0
            wave_percentages.delta = delta_power / total_power * 100;
            wave_percentages.theta = theta_power / total_power * 100;
            wave_percentages.alpha = alpha_power / total_power * 100;
            wave_percentages.beta = beta_power / total_power * 100;
            wave_percentages.gamma = gamma_power / total_power * 100;
        end
    else
        % For signals long enough, use pwelch
        try
            % Use pwelch with explicitly specified window size
            winSize = min(128, floor(length(signal)/2));  % Ensure window not larger than half signal length
            [pxx, f] = pwelch(signal, hamming(winSize), [], [], Fs);
            
            % Find indices for each frequency band
            delta_idx = (f >= delta_range(1)) & (f <= delta_range(2));
            theta_idx = (f > theta_range(1)) & (f <= theta_range(2));
            alpha_idx = (f > alpha_range(1)) & (f <= alpha_range(2));
            beta_idx = (f > beta_range(1)) & (f <= beta_range(2));
            gamma_idx = (f > gamma_range(1)) & (f <= gamma_range(2));
            
            % Calculate power in each band
            delta_power = sum(pxx(delta_idx));
            theta_power = sum(pxx(theta_idx));
            alpha_power = sum(pxx(alpha_idx));
            beta_power = sum(pxx(beta_idx));
            gamma_power = sum(pxx(gamma_idx));
            
            % Calculate total power
            total_power = delta_power + theta_power + alpha_power + beta_power + gamma_power;
            
            % Calculate percentages for each band
            if total_power > 0
                wave_percentages.delta = delta_power / total_power * 100;
                wave_percentages.theta = theta_power / total_power * 100;
                wave_percentages.alpha = alpha_power / total_power * 100;
                wave_percentages.beta = beta_power / total_power * 100;
                wave_percentages.gamma = gamma_power / total_power * 100;
            end
        catch ME
            warning('pwelch failed in calculate_wave_percentages: %s. Using default values.', ME.message);
        end
    end
end