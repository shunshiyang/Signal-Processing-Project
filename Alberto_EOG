%% Load EDF and XML Files  
edfFilename = 'R4.edf'; 
xmlFilename = 'R4.xml'; 
[hdr, record] = edfread(edfFilename); 
[events, stages, epochLength, annotation] = readXML(xmlFilename); 

%% Identify EOG Channel 
eog_channel = find(strcmp(hdr.label, 'EOGL')); % Adjust if needed
eog_signal = record(eog_channel, :); 
fs = hdr.samples(eog_channel); % Original sampling frequency (50 Hz)

%% Method 1: Butterworth Bandpass 4º (5-24.9 Hz)
low_cutoff = 5;
high_cutoff = 24.9;
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

%% Method 5: Cascade (5-24.9 Hz) 
% highpass filter(5Hz)
[b_high, a_high] = butter(4, 5/(fs/2), 'high');
% lowpass filter (24.9Hz)
[b_low, a_low] = butter(4, 24.9/(fs/2), 'low');
eog_temp = filtfilt(b_high, a_high, eog_signal);
eog_filtered5 = filtfilt(b_low, a_low, eog_temp);

%% Select First 3 Epochs (0-90 seconds)
epoch_samples = fs * epochLength;
start_idx = 1; % First epoch
end_idx = 3 * epoch_samples; % Third epoch

% Extract segments for all versions
eog_raw = eog_signal(start_idx:end_idx);
eog_filt1 = eog_filtered1(start_idx:end_idx);
eog_filt2 = eog_filtered2(start_idx:end_idx);
eog_filt3 = eog_filtered3(start_idx:end_idx);
eog_filt4 = eog_filtered4(start_idx:end_idx);
eog_filt5 = eog_filtered5(start_idx:end_idx);

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
Y_filt5 = abs(fft(eog_filt5)) / N;

%% Visualization over time
figure;

% Time Domain Plots
subplot(6, 1, 1);
plot(t, eog_raw);
title('Raw EOG (Epochs 1-3)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

subplot(6, 1, 2);
plot(t, eog_filt1);
title('Butterworth Filtered 4º (5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

hold on;
for epoch = 1:2
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(6, 1, 3);
plot(t, eog_filt2);
title('Butterworth Filtered 8º (5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

hold on;
for epoch = 1:2
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(6, 1, 4);
plot(t, eog_filt3);
title('Type I Chebyshev (5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

hold on;
for epoch = 1:2
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(6, 1, 5);
plot(t, eog_filt4);
title('Eliptic (5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

hold on;
for epoch = 1:2
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

subplot(6, 1, 6);
plot(t, eog_filt5);
title('Cascade (5-24.9 Hz)');
xlabel('Time (s)'); ylabel('Amplitude (μV)');
xlim([0, 90]); grid on;

hold on;
for epoch = 1:2
    xline(epoch * 30, '--r', sprintf('Epoch %d', epoch));
end
hold off;

%% Visualization over frecuency
figure;
subplot(6, 1, 1);
plot(f(1:N/2), Y_raw(1:N/2));
title('FFT (Raw)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(6, 1, 2);
plot(f(1:N/2), Y_filt1(1:N/2));
title('FFT (Butterworth 4)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(6, 1, 3);
plot(f(1:N/2), Y_filt2(1:N/2));
title('FFT (Butterworth 8)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(6, 1, 4);
plot(f(1:N/2), Y_filt3(1:N/2));
title('FFT (Chebyshev Tipo I)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(6, 1, 5);
plot(f(1:N/2), Y_filt4(1:N/2));
title('FFT (Eliptic)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

subplot(6, 1, 6);
plot(f(1:N/2), Y_filt5(1:N/2));
title('FFT (Cascade)');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
xlim([0, 50]); grid on;

% Adjust figure
set(gcf, 'Position', [100, 100, 1200, 1200]);
sgtitle('EOG Signal Processing: Comparison of Filtering Methods');
