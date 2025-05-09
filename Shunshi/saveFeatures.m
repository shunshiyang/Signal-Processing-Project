function saveFeatures(featureTable, subjectID, channelLabel, featureFolder)

    if ~exist(featureFolder, 'dir')
        mkdir(featureFolder);
    end
    
    channelLabel_clean = regexprep(channelLabel, '[^\w]', '_');

    outputFile = sprintf('%s/S%d_%s_Features.mat', featureFolder, subjectID, channelLabel_clean);

    save(outputFile, 'featureTable');
    fprintf('Features preserved to: %s\n', outputFile);
end