%% Load EDF and XML Files  
edfFilename = '../raw_data/R4.edf'; 
xmlFilename = '../raw_data/R4.xml'; 
[hdr, record] = edfread(edfFilename); 
[events, stages, epochLength, annotation] = readXML(xmlFilename); 

%% Identify EOG Channel 
eog_channel = find(strcmp(hdr.label, 'EOGL')); % Adjust if needed
eog_signal = record(eog_channel, :); 
fs = hdr.samples(eog_channel); % Original sampling frequency (50 Hz)

%% Method 1: Butterworth Bandpass 4º (5-24.9 Hz)
low_cutoff = 0.5;
high_cutoff = 24.5;
[b1, a1] = butter(4, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
eog_filtered1 = filtfilt(b1, a1, eog_signal);

%% Method 2: Butterworth Bandpass 8º (5-24.9 Hz)
[b2, a2] = butter(8, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
eog_filtered2 = filtfilt(b2, a2, eog_signal);

%% Method 3: Type I Chebyshev (5-24.9 Hz) 
[b3, a3] = cheby1(6, 0.1, [low_cutoff, high_cutoff] / (fs/2), 'bandpass');
eog_filtered3 = filtfilt(b3, a3, eog_signal);

%% Method 4: Eliptic (5-24.9 Hz) 
[b4, a4] = ellip(4, 0.1, 60, [low_cutoff, high_cutoff]/(fs/2), 'bandpass');
eog_filtered4 = filtfilt(b4, a4, eog_signal);

%% Select First 3 Epochs (0-90 seconds)
epoch_samples = fs * epochLength;
start_idx = 1; % First epoch
end_idx = 10 * epoch_samples; % Third epoch

% Extract segments for all versions
eog_raw = eog_signal(start_idx:end_idx);
eog_filt1 = eog_filtered1(start_idx:end_idx);
eog_filt2 = eog_filtered2(start_idx:end_idx);
eog_filt3 = eog_filtered3(start_idx:end_idx);
eog_filt4 = eog_filtered4(start_idx:end_idx);

% Time axis
t = (0:length(eog_raw)-1) / fs;

%% Compute FFTs
N = length(eog_raw);
f = (0:N-1) * (fs/N);

Y_raw = abs(fft(eog_raw)) / N;
Y_filt1 = abs(fft(eog_filt1)) / N;
Y_filt2 = abs(fft(eog_filt2)) / N;
Y_filt3 = abs(fft(eog_filt3)) / N;
Y_filt4 = abs(fft(eog_filt4)) / N;

%% Visualization over time
figure;

% Time Domain Plots
subplot(5, 1, 1);
plot(t, eog_raw);
title('Raw EOG (Epochs 1-3)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
grid on;

subplot(5, 1, 2);
plot(t, eog_filt1);
title('Butterworth Filtered 4º (0.5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
grid on;

hold on;
for epoch = 1:9
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(5, 1, 3);
plot(t, eog_filt2);
title('Butterworth Filtered 8º (0.5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
grid on;

hold on;
for epoch = 1:9
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(5, 1, 4);
plot(t, eog_filt3);
title('Type I Chebyshev (0.5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
grid on;

hold on;
for epoch = 1:9
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(5, 1, 5);
plot(t, eog_filt4);
title('Eliptic (0.5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
grid on;

hold on;
for epoch = 1:9
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

%% Visualization over frecuency
figure;
subplot(5, 1, 1);
plot(f(1:N/2), Y_raw(1:N/2));
title('FFT (Raw)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(5, 1, 2);
plot(f(1:N/2), Y_filt1(1:N/2));
title('FFT (Butterworth 4)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 30]); grid on;

subplot(5, 1, 3);
plot(f(1:N/2), Y_filt2(1:N/2));
title('FFT (Butterworth 8)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 30]); 
grid on

subplot(5, 1, 4);
plot(f(1:N/2), Y_filt3(1:N/2));
title('FFT (Chebyshev Tipo I)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 30]); grid on;

subplot(5, 1, 5);
plot(f(1:N/2), Y_filt4(1:N/2));
title('FFT (Eliptic)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 30]); grid on;

% Adjust figure
set(gcf, 'Position', [100, 100, 1200, 1200]);
sgtitle('EOG Signal Processing: Comparison of Filtering Methods');


%% Quantitative Filter Performance Analysis
% Calculate Artifact Rejection Ratio (ARR) and Signal-to-Artifact Improvement (SAI)

% Define frequency bands for artifacts
blink_band = [0, 5];    % Blink artifacts (0-5 Hz)
emg_band = [24.5, fs/2];  % EMG artifacts (>30 Hz, up to Nyquist)

% Initialize results table
filter_names = {'Butter4', 'Butter8', 'Cheby1', 'Ellip'};
results = table('Size', [4, 4], 'VariableTypes', {'string', 'double', 'double', 'double'}, ...
                'VariableNames', {'Filter', 'ARR_blink', 'ARR_emg', 'SAI'});
results.Filter = filter_names';

% Loop through each filtered signal
filtered_signals = {eog_filt1, eog_filt2, eog_filt3, eog_filt4};

for i = 1:4
    % Current filtered signal
    eog_filt = filtered_signals{i};
    
    % Compute residual (artifacts removed by the filter)
    residual = eog_raw - eog_filt;
    
    % Power spectral density of raw and residual signals
    [pxx_raw, f] = pwelch(eog_raw, [], [], [], fs);
    [pxx_res, ~] = pwelch(residual, [], [], [], fs);
    
    % Calculate power in artifact bands
    % Blink artifacts (0-5 Hz)
    blink_idx = (f >= blink_band(1)) & (f < blink_band(2));
    power_blink_raw = sum(pxx_raw(blink_idx));
    power_blink_res = sum(pxx_res(blink_idx));
    
    % EMG artifacts (>30 Hz)
    emg_idx = (f > emg_band(1)) & (f <= emg_band(2));
    power_emg_raw = sum(pxx_raw(emg_idx));
    power_emg_res = sum(pxx_res(emg_idx));
    
    % Artifact Rejection Ratio (ARR) in dB
    ARR_blink = 10 * log10(power_blink_raw / power_blink_res);
    ARR_emg = 10 * log10(power_emg_raw / power_emg_res);
        
    % Signal-to-Artifact Improvement (SAI)
    % SAR = 10*log10(signal_power / artifact_power)
    signal_power = var(eog_filt);
    artifact_power_raw = var(eog_raw - mean(eog_raw)); % Raw artifact power
    artifact_power_filt = var(residual);               % Residual artifact power
    
    SAR_raw = 10 * log10(signal_power / artifact_power_raw);
    SAR_filt = 10 * log10(signal_power / artifact_power_filt);
    SAI = SAR_filt - SAR_raw;
    
    % Store results
    results.ARR_blink(i) = ARR_blink;
    results.ARR_emg(i) = ARR_emg;
    results.SAI(i) = SAI;
end

%% Display Results
disp('Filter Performance Metrics:');
disp(results);

%% Plot Results
figure;
subplot(3, 1, 1);
b1 = bar(results.ARR_blink);
title('Artifact Rejection Ratio (Blink Band 0.5-5 Hz)');
ylabel('ARR (dB)');
set(gca, 'XTickLabel', filter_names);
ylim([0 15]); % Extend y-axis by 20% above max value
grid on;

for i = 1:length(b1.YData)
    text(i, b1.YData(i) + 0.3, sprintf('%.2f', b1.YData(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

subplot(3, 1, 2);
b2 = bar(results.ARR_emg);
title('Artifact Rejection Ratio (EMG Band 24.5-25 Hz)');
ylabel('ARR (dB)');
set(gca, 'XTickLabel', filter_names);
ylim([0 3]); % Extend y-axis by 20% above max value
grid on;

for i = 1:length(b2.YData)
    text(i, b2.YData(i) + 0.05, sprintf('%.2f', b2.YData(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

subplot(3, 1, 3);
b3 = bar(results.SAI);
title('Signal-to-Artifact Improvement (SAI)');
ylabel('Improvement (dB)');
set(gca, 'XTickLabel', filter_names);
ylim([0 15]); % Extend y-axis by 20% above max value
grid on;

for i = 1:length(b3.YData)
    text(i, b3.YData(i) + 0.3, sprintf('%.2f', b3.YData(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

sgtitle('Quantitative Filter Performance Comparison');



%% === Compute Peak Detection Parameters Automatically ===
% We'll use eog_filt1 here
signal_to_analyze = eog_filt4;

% Step 1: Get signal stats
min_val = min(signal_to_analyze);
max_val = max(signal_to_analyze);
range_val = max_val - min_val;
mean_val = mean(signal_to_analyze);
std_val = std(signal_to_analyze);

% Step 2: Set thresholds
min_peak_height = mean_val + 0.5 * std_val;   % Conservative threshold
min_peak_distance = 0.3 * fs;  % 0.3 seconds apart (adjustable)

% Optional: Plot signal stats
fprintf('Min: %.2f µV | Max: %.2f µV | Range: %.2f µV | Mean: %.2f | STD: %.2f\n', ...
    min_val, max_val, range_val, mean_val, std_val);
fprintf('Threshold set at: %.2f µV | Min peak distance: %.1f samples (%.2f s)\n', ...
    min_peak_height, min_peak_distance, min_peak_distance/fs);

%% === Detect Peaks ===
[peaks, locs] = findpeaks(signal_to_analyze, ...
    'MinPeakHeight', min_peak_height, ...
    'MinPeakDistance', round(min_peak_distance));

%% === Visualize Signal and Peaks ===
figure;
plot(t, signal_to_analyze, 'b'); hold on;
plot(t(locs), peaks, 'ro', 'MarkerFaceColor', 'r');
title('Peak Detection on Filtered EOG Signal (Butterworth 4º)');
xlabel('Time (s)'); ylabel('Amplitude (µV)');
grid on;
legend('Filtered Signal', 'Detected Peaks');

%% === Optional: Display Peak Statistics ===
fprintf('Number of detected peaks: %d\n', length(peaks));
fprintf('Minimum peak amplitude: %.2f µV\n', min(peaks));
fprintf('Mean peak amplitude: %.2f µV\n', mean(peaks));

% Optional histogram of peak amplitudes
figure;
histogram(peaks, 20);
title('Histogram of Detected Peak Amplitudes');
xlabel('Amplitude (µV)');
ylabel('Count');
grid on;