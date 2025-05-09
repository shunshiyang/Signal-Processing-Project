function params = configureParams()
    params = struct();
    
    params.dataFolder = 'Data/';
    params.FeatureFolder = 'Features/';
    
    params.subjectIDs = 1:10;
    params.numEpochsToPlot = inf; 
    
    params.toPlotFiltering = false;
    params.toPlotRR_QRS = false;
    params.toPlotRRAnalysis = false;
    params.toPlotResults = false;
    
    params.lowPassCutoff = 150;
    params.highPassCutoff = 0.67;
    params.QRS_LowCutoff = 10;
    params.QRS_HighCutoff = 25;
    
    params.vlf_band = [0.003, 0.04];
    params.lf_band = [0.04, 0.15];
    params.hf_band = [0.15, 0.4];
    params.Fs_interp = 4; 
end