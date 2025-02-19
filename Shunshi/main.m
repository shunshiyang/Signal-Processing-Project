%% load edf and xml files
edfFilename = 'R4.edf';
xmlFilename = 'R4.xml';
[hdr, record] = edfread(edfFilename);
[events, stages, epochLength,annotation] = readXML(xmlFilename);

%%
numberOfEpochs = length(record(3,:)')/(30*hdr.samples(3))
%%
eeg_channel = find(strcmp(hdr.label, 'EEGsec'));
fs = hdr.samples(eeg_channel);
eeg_signal = record(eeg_channel, :);
%% Preprocess
processed_data = preprocess_eeg(eeg_signal, fs, 50);

%% Result analysis
%Choose on epoch to do the comparison
epochNumber = 1;
epochStart = (epochNumber * fs * 30) + 1;
epochEnd = epochStart + (30 * fs) - 1;
epoch_signal = eeg_signal(epochStart:epochEnd);
processed_epoch = processed_data(epochStart:epochEnd);

figure;
t = (0:length(epoch_signal)-1)/fs;

% Original signal
subplot(2,1,1);
plot(t, epoch_signal);
title('Original EEG Signal');
xlabel('Time (s)');
ylabel('Amplitude (μV)');

% Processed signal
subplot(2,1,2);
plot(t, processed_epoch);
title('Preprocessed EEG Signal');
xlabel('Time (s)');
ylabel('Amplitude (μV)');

figure;
% Power Spectrum of Original Signal
subplot(2,1,1);
[pxx_orig,f] = pwelch(epoch_signal,[],[],[],fs);
plot(f,10*log10(pxx_orig));
title('Power Spectrum of Original Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

% Power Spectrum of Processed Signal
subplot(2,1,2);
[pxx_proc,f] = pwelch(processed_epoch,[],[],[],fs);
plot(f,10*log10(pxx_proc));
title('Power Spectrum of Processed Signal');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

figure;
% Spectrogram of Original Signal
subplot(2,1,1);
spectrogram(epoch_signal,hamming(256),128,256,fs,'yaxis');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Original Signal');
colorbar;
ylabel(colorbar,'Power/Frequency (dB/Hz)');

% Spectrogram of Processed Signal
subplot(2,1,2);
spectrogram(processed_epoch,hamming(256),128,256,fs,'yaxis');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Processed Signal');
colorbar;
ylabel(colorbar,'Power/Frequency (dB/Hz)');