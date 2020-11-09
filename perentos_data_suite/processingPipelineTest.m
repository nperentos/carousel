% this is a script that runs the first basic processing steps

%fileBase = 'NP50_2020-07-01_17-53-27';
%fileBase = 135
tic;
fileBase = 'NP51_2020-05-25_14-34-24';

%% **** MANUAL STEPS *****
    initialiseVideoAnalysis(fileBase);% this step is manual

%% someLFP and ADC
    convertDataCarousel(fileBase);
    % modify the dat file
    getGoodChannels(fileBase); 
    findCommonMode(fileBase);
    generateAuxVars(fileBase);
    % do you need to update spreadsheet??

%% realLFP?????

%% Units??????
    
%% VIDEO    
    
    runVideoAnalysis(fileBase);
    % user must check results here before next step
    % update spreadheet with scores of video analusis 0 super-bad 1 good 2
    % bad
    getVidVars(fileBase);
    
%% COMBINE BEHAVIORAL VARIABLLES FROM LFP AND VIDEO
    combineVidAuxVars(fileBase);
    
    
%%    
    elapsed_time = toc;        
    display(['the pipeline finished in ', num2str(elapsed_time), ' seconds']);