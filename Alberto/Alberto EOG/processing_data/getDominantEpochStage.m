function epochStages = getDominantEpochStage(stages, fs, epochLen)
    % Convert into dominant sleep stages per epoch
    samplesPerEpoch = epochLen; 
    numEpoch = floor(length(stages) / samplesPerEpoch);
    
    epochStages = zeros(1, numEpoch);
    for epoch = 1:numEpoch
        startIdx = (epoch-1)*samplesPerEpoch + 1;
        endIdx = epoch*samplesPerEpoch;
        epochWindow = stages(startIdx:endIdx);
        
        %The stage that appears the most (dominant) in a 30 second window becomes the epoch stage
        %  Resolve ties by picking lowest number (deeper sleep) 
        %so REM(0) is picked over Wake(5) for example 
        uniqueValues = unique(epochWindow);
        counts = histc(epochWindow, uniqueValues);
        maxCount = max(counts);
        dominantCandidates = uniqueValues(counts == maxCount);
        epochStages(epoch) = min(dominantCandidates);
    end
end