function [hdr, record, stages, events, epochLength] = loadData(subj, dataFolder)
    patientIDs = {'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10'};
    patientID = patientIDs{subj};
    edfFilename = [dataFolder patientID '.edf'];
    xmlFilename = [dataFolder patientID '.xml'];
    
    try
        [hdr, record] = edfread(edfFilename);
    catch
        warning('Couldnt read the document: %s', edfFilename);
        hdr = []; record = []; stages = []; events = []; epochLength = [];
        return;
    end
    
    try
        [events, stages, epochLength, annotation] = readXML(xmlFilename);
    catch
        warning('Couldnt read the document': %s', xmlFilename);
        hdr = []; record = []; stages = []; events = []; epochLength = [];
        return;
    end
    
    eeg_indices = find(contains(hdr.label, 'EEG', 'IgnoreCase', true));
    if isempty(eeg_indices)
        warning('No EEG channel');
        hdr = []; record = []; stages = []; events = []; epochLength = [];
        return;
    end
    
    stages(stages == 1) = 2;  
end