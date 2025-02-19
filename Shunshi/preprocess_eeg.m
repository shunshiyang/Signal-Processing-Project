function processed_data = preprocess_eeg(eeg_data, fs, powerline_freq)
    % Main preprocessing function for EEG signal
    processed_data = eeg_data;
    
    % 1. Remove baseline drift
    processed_data = remove_baseline_drift(processed_data, fs);
    
    % 2. Remove muscle noise
    processed_data = remove_muscle_noise(processed_data, fs);
    
    % 3. Remove powerline interference
    processed_data = remove_powerline_interference(processed_data, fs, powerline_freq);
end

function filtered_data = remove_baseline_drift(eeg_data, fs)
    % Remove baseline drift using a highpass filter (0.5 Hz)
    order = 4;
    cutoff = 0.5; 
    [b, a] = butter(order, cutoff/(fs/2), 'high');
    filtered_data = filtfilt(b, a, double(eeg_data));
end

function filtered_data = remove_muscle_noise(eeg_data, fs)
    % Remove high frequency muscle noise using a lowpass filter (35 Hz)
    order = 4;
    cutoff = 40;
    [b, a] = butter(order, cutoff/(fs/2), 'low');
    filtered_data = filtfilt(b, a, double(eeg_data));
end

function filtered_data = remove_powerline_interference(eeg_data, fs, powerline_freq)
    % Remove powerline interference using a narrow bandstop filter
    order = 2;
    bandwidth = 1;
    low_cutoff = (powerline_freq - bandwidth/2)/(fs/2);
    high_cutoff = (powerline_freq + bandwidth/2)/(fs/2);
    [b, a] = butter(order, [low_cutoff high_cutoff], 'stop');
    filtered_data = filtfilt(b, a, double(eeg_data));
end