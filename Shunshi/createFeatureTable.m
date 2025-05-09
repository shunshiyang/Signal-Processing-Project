function featureTable = createFeatureTable(features, sleepStages_perEpoch)
    
    numEpochs = length(sleepStages_perEpoch);
    processedStages = processSleepStageLabels(sleepStages_perEpoch, numEpochs);
    
    
    featureTable = table( ...
        features(:, 1), ...  % deltaFeature
        features(:, 2), ...  % thetaFeature
        features(:, 3), ...  % alphaFeature
        features(:, 4), ...  % betaFeature
        features(:, 5), ...  % meanFeature
        features(:, 6), ...  % stdFeature
        features(:, 7), ...  % varianceFeature
        features(:, 8), ...  % kurtosisFeature
        features(:, 9), ...  % skewnessFeature
        features(:, 10), ... % sef95
        features(:, 11), ... % kComplexCount
        features(:, 12), ... % sleepSpindleCount
        categorical(processedStages(:)), ...  
        'VariableNames', {'Delta', 'Theta', 'Alpha', 'Beta', ...
                          'Mean', 'StdDev', 'Variance', 'Kurtosis', ...
                          'Skewness', 'SEF95', 'KComplexCount', 'SpindleCount', ...
                          'SleepStage'});
end

function processedStages = processSleepStageLabels(sleepStages_perEpoch, numEpochs)
    
    processedStages = sleepStages_perEpoch;
    
    for epoch = 1:numEpochs
        stageVal = str2double(string(sleepStages_perEpoch(epoch)));
        
        if isnan(stageVal)
            processedStages(epoch) = categorical("Unknown");
        else
            switch stageVal
                case 0
                    processedStages(epoch) = categorical("REM");
                case 2
                    processedStages(epoch) = categorical("N3");
                case 3
                    processedStages(epoch) = categorical("N2");
                case 4
                    processedStages(epoch) = categorical("N1");
                case 5
                    processedStages(epoch) = categorical("Wake");
                otherwise
                    processedStages(epoch) = categorical("Unknown");
            end
        end
    end
end