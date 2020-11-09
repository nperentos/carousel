function [session, behavior] = loadSession(fileBase)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [session, behavior] = loadSession(fileBase) loads variables 
% from processed data from a carousel session.
% session.mat
% behavior.mat
% rez.mat
% AND ANY OTHER ONES THAT WE WILL PRODUCE ALONG THE WAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PRE
    processedPath = getfullpath(fileBase);
    cd(processedPath);
    pth = fullfile(getfullpath(fileBase),[fileBase,'.lfp']);
    par = LoadXml([pth(1:end-4),'.xml']);
    SR = par.lfpSampleRate;
    %[~,settings,tScale] = getLFP(fileBase);
    s = dir([processedPath,fileBase,'.lfp']);
    tMax = s.bytes/(SR * par.nChannels*2);% 2 b/c data is 16 bit ints => 2 bytes per datapoint
    tScale = [1/SR:1/SR:tMax]; 
  

    
%% LOAD STUFF AND 'POINTERS'
    try
        load('session.mat');
        m = matfile('behavior.mat');
        behavior.data = m;
        behavior.name = m.name;
        behavior.info = m.info;
        %SPIKES???
        %OTHER???
        %CHECKS FOR EACH???
    catch
        display('note this function works with beavior.mat not interim_behavior...');
        error('didnt manage to load session properly - please check');
    end


    
%%    