function featureRow = extractEOGFeatures(epochSignal, fs, epochLength)
    featureRow = struct();
    try
        %% --- Non-Parametric Features (Already discussed) ---  27 features
        %TIME DOMAIN FEATURES 12 features
        
        featureRow.Mean = mean(epochSignal); %1
        featureRow.Std = std(epochSignal); %2
        featureRow.Variance = var(epochSignal); %3 
        featureRow.RMS = rms(epochSignal); %4
        featureRow.ZeroCrossings = sum(abs(diff(sign(epochSignal)))) / (2 * epochLength);%5
        featureRow.SlopeSignChanges = sum(diff(sign(diff(epochSignal))) ~= 0);%6
        featureRow.WaveformLength = sum(abs(diff(epochSignal)));%7
        featureRow.IEMG = sum(abs(epochSignal));%8
        featureRow.MAV = mean(abs(epochSignal));%9
        featureRow.SSI = sum(epochSignal.^2);%10
        featureRow.Kurtosis = kurtosis(epochSignal);%11 
        featureRow.Skewness = skewness(epochSignal);%12 

        % Welch PSD FREQUENCY-DOMAIN FEATURES 15 
        [pxx, f] = pwelch(epochSignal, [], [], [], fs); 
        featureRow.TotalPower = bandpower(pxx, f, [0.5 25], 'psd'); %13 
        featureRow.DeltaPower = bandpower(pxx, f, [0.5 4], 'psd');%14 
        featureRow.ThetaPower = bandpower(pxx, f, [4 8], 'psd');%15 
        featureRow.AlphaPower = bandpower(pxx, f, [8 13], 'psd'); %16 
        featureRow.BetaPower = bandpower(pxx, f, [13 25], 'psd'); %17 

        featureRow.AlphaThetaRatio = featureRow.AlphaPower / (featureRow.ThetaPower + eps); %18 
        featureRow.BetaAlphaRatio = featureRow.BetaPower / (featureRow.AlphaPower + eps); %19 
        
        %% --- Spectral Features ---
        psd_norm = pxx / sum(pxx);  
        featureRow.SpectralEntropy = -sum(psd_norm .* log2(psd_norm + eps)); %20

        cumulativePower = cumtrapz(f, pxx);
        totalPower = cumulativePower(end);%21 
        featureRow.SEF95 = f(find(cumulativePower >= 0.95 * totalPower, 1)); %22 
        featureRow.SEF50 = f(find(cumulativePower >= 0.50 * totalPower, 1)); %23 

        % Spectral Centroid
        featureRow.SpectralCentroid = sum(f .* pxx) / sum(pxx); %24 
    
        % Spectral Flatness 
        featureRow.SpectralFlatness = geomean(abs(pxx)) / mean(abs(pxx)); %25
    
        % Spectral Roll-Off (85%)
        threshold = 0.85 * sum(pxx);
        featureRow.SpectralRollOff = find(cumsum(pxx) >= threshold, 1, 'first');  %26 

        signalDiff1 = diff(epochSignal);
        signalDiff2 = diff(signalDiff1);
        var0 = var(epochSignal);
        var1 = var(signalDiff1);
        var2 = var(signalDiff2);
        featureRow.HjorthActivity = var0; %27
        featureRow.HjorthMobility = sqrt(var1 / var0); %28
        featureRow.HjorthComplexity = sqrt(var2 / var1) / featureRow.HjorthMobility; %29 


        %% --- Eye Movement Features (Blink Rate, REM, SEM, Movement Density) ---

        % Parameters
        timeInSeconds = epochLength;  % e.g., 30 seconds
        thresholdMovement = 3 * std(epochSignal);  % ~121.6 µV with your data
        
        % --- Movement Density ---
        movementEvents = sum(abs(diff(epochSignal)) > thresholdMovement);
        featureRow.MovementDensity = movementEvents / timeInSeconds;
        
        % --- Blink Detection ---
        [b,a] = ellip(4, 0.5, 40, [0.5 5]/(fs/2));  % 0.5–5 Hz filter for blink band
        blinkSignal = filtfilt(b, a, epochSignal);

        % Compute adaptive min peak height
        mean_val = mean(blinkSignal);
        std_val = std(blinkSignal);
        min_peak_heigh = mean_val + 0.5 * std_val;  % safer threshold
        min_peak_distance = 0.3 * fs;  % 300 ms minimum blink interval

        [blinkPeaks, blinkLocs] = findpeaks(blinkSignal, ...
            'MinPeakHeight', min_peak_heigh, ...
            'MinPeakDistance', round(min_peak_distance));
        featureRow.BlinkRate = length(blinkPeaks) / timeInSeconds;
        
        % --- Slow Eye Movements (SEM) ---
        [b_sem, a_sem] = ellip(4, 0.5, 40, [0.5 1.0]/(fs/2));  % 0.5–1.0 Hz
        semSignal = filtfilt(b_sem, a_sem, epochSignal);
        semSTD = std(semSignal);  % Use adaptive threshold
        [semPeaks, ~] = findpeaks(semSignal, 'MinPeakHeight', semSTD);
        featureRow.SEM_Rate = length(semPeaks) / timeInSeconds;
        
        % --- Rapid Eye Movements (REM) --- 
        timeInSeconds = epochLength;  % e.g., 30 seconds
        [b_rem, a_rem] = ellip(4, 0.5, 40, [1.0 10]/(fs/2));  % 1–10 Hz
        remSignal = filtfilt(b_rem, a_rem, epochSignal);
        remSignal = remSignal(:);
        remVelocity = abs([0; diff(remSignal)]);
        remThreshold = 2.5 * std(remVelocity);  % Slightly more sensitive than 3*std
        remPeaks = find(remVelocity > remThreshold);
        featureRow.REM_Rate = length(remPeaks) / timeInSeconds;
            
    catch ME
        
        
    end
end
