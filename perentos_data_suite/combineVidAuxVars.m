function combineVidAuxVars(fileBase,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combineVidAuxVars(fileBase) puts together the outputs from generateAuxVars
% (behavior_interim.mat) and getVidVars (all_video_results.mat)
% into a big structure called session

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('combining video and auxiliary variables into one matrix at lfp rate');

%% PRE

    options = {'ifPlot',1,'noVidFlag',0};
    options = inputparser(varargin,options); 

    processedPath = getfullpath(fileBase);
    cd(processedPath);
    pth = fullfile(getfullpath(fileBase),[fileBase,'.lfp']);
    par = LoadXml([pth(1:end-4),'.xml']);
    SR = par.lfpSampleRate;
    %[~,settings,tScale] = getLFP(fileBase);
    s = dir([processedPath,fileBase,'.lfp']);
    tMax = s.bytes/(SR * par.nChannels*2);% 2 b/c data is 16 bit ints => 2 bytes per datapoint
    tScale = [1/SR:1/SR:tMax]; 
    
    
    
%% IMPORT INTERIM DATA
    if exist(fullfile(processedPath,'all_video_results.mat')) == 2
        load(fullfile(processedPath,'all_video_results.mat'));
    else
        if exist(fullfile(getFullPath(fileBase),'video')) == 7
            error('there is a video folder but all video results.mat was not created. Please investigate');
        else
           noVidFlag = 1;
           % for now throw error till I decide what to do
           error('there is a video folder but all video results.mat was not created. Please investigate');
        end
    end                    
    load(fullfile(processedPath,'behavior_interim.mat'));
    load(fullfile(processedPath,'session.mat'));
    

    
%% UPSAMPLE THE VIDEO VARIABLES
    if ~isempty(videoFeatures)
        disp('upsampling video features to LFP rate... CANT WE AVOID THIS??');
        n1 = sum(session.helper.vidPulses); n2 = length(videoFeatures.data);
        % its possible that there are stray pulses in the ADC channel (due
        % to drift at onset or offset of pulse trains so we should eliminate them
        % we eliminate by finding pulses that are not 1/30s apart
        if n1 > n2
            tmp1 = find(session.helper.vidPulses ~= 0); % pulse positions
            tmp2 = diff(tmp1);
            tmp3 = find(tmp2 > 31 & tmp2 < 35); % fixed to 30 Hz unfortunately
            %special case for the last pulse
            if tmp3(1) == 1 & tmp1(tmp3(end)) ~= tmp1(end)% first pulse is good so add another pulse at the end to account for the diffing? This is untested however
                tmp3 = [1, tmp3+1];
            end
            tmp1 = tmp1(tmp3);
            mockPulses = session.helper.vidPulses; mockPulses(mockPulses ~= 0) = 0;
            mockPulses(tmp1) = 1;
            session.helper.vidPulses = mockPulses;
            clear tmp1 tmp2 tmp3
        end
            
        n1 = sum(session.helper.vidPulses); n2 = length(videoFeatures.data);
        n = abs(n1-n2);
        if n > 2; warning('the number of pulses and video-derived behavior data differ by more than two frames. Please investigate'); end
        % still, it is possible that carousel protocol was terminated
        % abnormaly. We will need special cases hardcoded if this was the
        % case
        % NP43_2019-12-08_18-33-23: OE recordng was terminated by accident
        % so there are more frames (~900 extra) than vid pulses
        clear tmp ztmp
        if strcmp(fileBase,'NP43_2019-12-08_18-33-23')
            % extra frames at the END to dropped
            videoFeatures.data = videoFeatures.data(1:sum(session.helper.vidPulses),:);
            tmp = interp1(tScale(session.helper.vidPulses==1),videoFeatures.data,tScale,'linear','extrap');
        elseif strcmp(fileBase,'NP53_2020-06-11_14-09-15')
            % extra frames at the START to drop (OE was not setup o accept
            % network event so did not start through matlab on time and was
            % started manualy instead?)
            videoFeatures.data = videoFeatures.data(end-sum(session.helper.vidPulses)+1:end,:);
            tmp = interp1(tScale(session.helper.vidPulses==1),videoFeatures.data,tScale,'linear','extrap');
        elseif strcmp(fileBase,'NP43_2019-12-11_16-08-24')
            % first we assume that n frames where collected before the
            % start of OE rec and so no vid pulses are available. So we rim
            % by n the videoFeatures
            videoFeatures.data(1:n,:) = [];
            % but there are loose pulses for the first 2.3 seconds (arduino
            % error?) for which we do have pulses but we do not want to
            % keep so we will trim those from both vidPulses as well as
            % videoFeatures
            pls = find(session.helper.vidPulses);
            IDX = find(diff(pls) > mode(diff(pls))+10);pls(IDX)  
            session.helper.vidPulses(1:pls(IDX)) = 0;
            videoFeatures.data(1:IDX,:) = [];            
            tmp = interp1(tScale(session.helper.vidPulses==1),videoFeatures.data,tScale,'linear','extrap');
        else
            tmp = interp1(tScale(session.helper.vidPulses==1),[zeros(n,size(videoFeatures.data,2)); videoFeatures.data],tScale,'linear','extrap');
        end
        ztmp = zscore(tmp);
        figure; imagesc(tScale,[],ztmp');caxis([prctile(ztmp(:),[1 99])]);        
        xlabel('time (s)'); ylabel('feature (#)'); title('z-scored behavioral feature activities');
        stampFig(fileBase);
        
    else
        display('no video features available!');
    end

    
    
%% MERGE BEHAVIORAL VARIABLES
    %     if ~isempty(videoFeatures)
    %         behavior.data = [behavior.data',tmp];
    %         behavior.name = [behavior.name;videoFeatures.name];
    %         behavior.info = videoFeatures.info;
    %     else
    %         behavior.data = [behavior.data'];
    %         behavior.name = [behavior.name];           
    %     end
    if ~isempty(videoFeatures)
        data = [behavior.data',tmp];
        name = [behavior.name;videoFeatures.name];
        info = videoFeatures.info;
    else
        data = [behavior.data'];
        name = [behavior.name];   
        info = 'no info here';
    end    


    
%% SAVE 
    fprintf(['saving auxiliary and video variables into behavior.mat at lfp rate...']);
    % save(fullfile(processedPath,'behavior.mat'),'behavior','-v7.3');
    save(fullfile(processedPath,'behavior.mat'),'data','name','info','-v7.3');
    % update session variablee
    save(fullfile(processedPath,'session.mat'),'session','-v7.3');
    print(fullfile(getfullpath(fileBase),'all_behavioral_variables.jpg'),'-djpeg');
    fprintf('DONE\n');
    
    
    
% %% GENERATE POSITION RESOLVED FIGURES FOR EACH BEHAVIORAL VARIABLE
%     disp('generating position resolved plots for behavioral variables...');
%     for i = [2 3 6 7 8 9 10 11 12 13 14 15]
%         carouselPlotVar(fileBase,pppNames{i},1);
%         pause(1); close;
%     end    
