function [data,settings,tScale] = getLFP(sessionFolder, dataPath,ch)
% takes folder name and loads the data
% works by default with NPs folder structure i.e. my data path is hardcoded
% can take other data paths if necessary in the optional dataPath input
% takes sessionFolder e.g. NP3_2018-04-11_19-37-06
% loads in workspace
%   lfp data into workspace as channels x timepoints
%   settings xml
%   sampling rate
%   xaxis time scale

defaultDataPath = '/storage2/perentos/data/recordings/';

%% ensure correct inputs

    % elseif nargin > 2
    %     error('too many input arguments -- max two')
if nargin > 0
    TF = isstrprop(sessionFolder(3:4),'digit');
    % if exist('processed','file'); pmod = '/processed'; else pmod = []; end
    if exist(fullfile(defaultDataPath, sessionFolder(1:4),'/',sessionFolder,'processed'),'file'); pmod = '/processed'; else pmod = []; end
    if sum(TF) == 2
        LFPpth = [defaultDataPath, sessionFolder(1:4),'/',sessionFolder,pmod,'/',sessionFolder,'.lfp'];
    elseif TF(1) == 1
        LFPpth = [defaultDataPath, sessionFolder(1:3),'/',sessionFolder,pmod,'/',sessionFolder,'.lfp'];
    else 
        disp('I was expecting a sessionFolder to start with NP (or other 2 initials characters followed by 1 or 2 digits');
        disp 'EXAMPLE:  NP3_2018-04-11_19-37-06'
        error('fix sessionFolder input variable and try again');e
    end
end

if nargin > 1 && ~isempty(dataPath)% alternative data path supplied
    defaultDataPath = dataPath;
    LFPpth = fullfile(defaultDataPath,sessionFolder,'processed',[sessionFolder,'.lfp']);
end

%%
disp 'loading binary ..lfp...';
if nargin == 3 % a subset of channels is requested
    data = load_binary(LFPpth,ch);
else
    data = load_binary(LFPpth);
end
disp 'convert to double...';
data = double(data);
settings = xml2struct([LFPpth(1:end-4),'.xml']); %grab sampling rate of LFP
SR = str2num(settings.parameters.fieldPotentials.lfpSamplingRate.Text);
settings = xml2struct([LFPpth(1:end-4),'.oe.xml']); % grab gain values
for i = 1:size(data,1)
    try
        gain(i) = str2num(settings.SETTINGS.SIGNALCHAIN{1,1}.PROCESSOR{1,1}.CHANNEL_INFO.CHANNEL{1,i}.Attributes.gain);
    catch
        gain(i) = str2num(settings.SETTINGS.SIGNALCHAIN.PROCESSOR{1,1}.CHANNEL_INFO.CHANNEL{1,i}.Attributes.gain);
    end         
end
data = data.*gain';
tScale = 1/SR:1/SR:size(data,2)/SR;
%settings = xml2struct([LFPpth(1:end-4),'.xml']); %grab sampling rate of LFP

    
