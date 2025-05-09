function sleepStages_perEpoch = processSleepStages(stages, Fs, epochLen)
    
    sleepStages_perEpoch = getDominantEpochStage(stages, Fs, epochLen);
    
    sleepStages_perEpoch = categorical(sleepStages_perEpoch);
end