function out = getMeanSpikePerPeriod(fileBase,clu_idx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getMeanSpikePerPeriod breaks a carousel recording into periods that
% represent single trials (breaks inbetween trial blocks are also treated
% as distint periods) and computes mean cluster spike waveform for each
% period. Its main use is to exclude clusters based on electrode drift as a
% function of time.
%
% INPUT:  fileBase: recording folder name
%         clu_idx:  optional to request subset of clusters - default is all
% OUTPUT: out.wvs: 3D matrix with period x waveamplitude x cell
%         out.cluID the cluster ID
%         out.dur duration of each of the periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    clu_idx = []; % process all cells
    ifSave = 1; % since all clusters will be processed save output in .mat    
end


%% GET ALL NECESSARY DATA INTO WORKSPACE AND EXTRACT 'PERIODS'
% import the kilosort data
    processedPath = getfullpath(fileBase);
    cd(processedPath);
    params.excludeNoise = 1;
    sp = loadKSdir(pwd,params);
    SR = sp.sample_rate;
    % check if sp.cids is sorted
    if ~issorted([sp.cids])
        disp('the cluster IDs are not sorted!');
        % only two variables to sort, cids and cgs
        I = sort(sp.cids);
        sp.cids = sp.cids(I);
        sp.cids = sp.cgs(I);
    end
    
% get the dat file and memory map it
    Res = round(sp.st * SR); % in samples
    Clu = sp.clu;
    chmap = readNPY(fullfile(pwd,'channel_map.npy')); % this is useful for
    FileName = dir(fullfile(pwd,'*.dat'));
    FileName = fullfile(FileName.folder,FileName.name);
    nSamples = FileLength(FileName)/2/sp.n_channels_dat;
    nChannels = sp.n_channels_dat;
    mmap = memmapfile(FileName, 'format',{'int16' [nChannels nSamples] 'x'},'offset',0,'repeat',1);
    
% load the necessary peripherals to construct events
    load peripherals; p = peripherals; clear peripherals;
    trial_end = p.events.atStart; % arriving at the eye position    
    load peripheralsPP; %ppp
    carouselSpeed = ppp(find(strcmp('carouselSpeed',pppNames) == 1),:);% find carouselSpeed from ppp
    tScale = ppp(find(strcmp('tScale',pppNames) == 1),:);% find carouselSpeed from ppp    
    first_start = find(carouselSpeed> 0.5,1,'first');

% build periods using awkward procedure
    % MISSING! specifically the very start, the start of the first trial, the end of the
    % middle trial (total trials/2 -1)and the very end of the rec
    % lets do a dirty fix for now (MUST FIX ARDUINO CODE!!!!)
    periods = [1; first_start; trial_end;];
    periods = [[periods] [periods(2:end);length(tScale)]];
    % we still need to inject the middle trial ending as well as start of
    % midbreak end (the two are the same but in different column)
    zper = zscore(diff(periods,[],2));
    for i = 1:6
        [pos,iN] = find(zper > i);
        if length(iN) == 3
            break
        end
    end
    
    % keep the middle occurance
    tr = pos(2);
    
    % find the last movement point of the carousel after tr and before tr+1
    carouselSpeedTmp = carouselSpeed(1:periods(tr,2));
    mid_trial_end = find(carouselSpeedTmp > 0.5,1,'last');
    % inject the missing trial in the middle
    periods2 = [periods(1:tr,:); [mid_trial_end mid_trial_end]; periods(tr+1:end,:)];   
    % now swap the stop times of tr and tr + 1
    periods2([tr:tr+1],2) = periods2([tr+1,tr],2);    

% plot to check that we are selecting the right periods    
    close all; figure; %subplot(211);
    plot(tScale,carouselSpeed );

    hold on;
    plot([tScale(unique(periods2)); tScale(unique(periods2))], repmat(ylim,length(unique(periods2)),1)','--r','linewidth',1);
    
    periods = periods2.*(SR/1000);
    
    
%% EXTRACT SPIKES FOR EACH CLUSTER AND EACH PERIOD
    if isempty(clu_idx)
        clu_idx = sp.cids;
    end
    
    parfor i = 1:length(clu_idx)
        display(['unit',num2str(i),'/',num2str(length(clu_idx))])
        myT = Res(find(Clu == clu_idx(i)));% lets fix the subset request later
        Lag = [-15:40];
        myT(myT<1-Lag(1) | myT>nSamples-Lag(end))=[]; % eliminate spikes to close to start/end of file
        nLag = length(Lag);
        nmyT = length(myT);
        myTLagged = reshape(bsxfun(@plus,myT, Lag),[],1); % idxs into dat file for each spike
        Amp = mmap.Data.x(:,myTLagged);
        Amp = reshape(Amp, nChannels, nmyT, nLag);
        Amp = Amp(chmap+1,:,:); % remove irrelevant channels
        mn = repmat(mean(Amp,3),[1, 1, size(Amp,3)]); % all the waveforms across whole electrode space
        mAmp = double(Amp)-mn; % mean corected spike waveforms   
       %mAmp is channel X n_of_spikes X timesamples
        mW = squeeze(mean(mAmp,2));
        % sdW = squeeze(std(mAmp,[],2));



        % find max amp waveform/channel
        % looks for non inverted waveforms too
        %[~,cluCh(i)] =max(max(-mW'));
        [tmp1Amp,tmp1Ch] =max(max(-mW'));
        [tmp2Amp,tmp2Ch] =max(max(mW'));

        [a,b] = max([tmp1Amp tmp2Amp]);
        if b == 1
            cluCh(i) = tmp1Ch;
        else
            cluCh(i) = tmp2Ch;
        end

        % now we know cluster-channel assignment, lets go back to the mean
        % corrected spikes and extract only those within periods of
        % interest
        for per = 1:length(periods)
            tmp_res = find(myT > periods(per,1) & myT < periods(per,2));
            if length(tmp_res) < 2
                out(i).wv(per,:) = nan(1,nLag);
            else
                out(i).wv(per,:) = mean(sq(mAmp(cluCh(i),tmp_res,:)),1);
            end
            out(i).dur(per) = diff(periods(per,:))/(SR);
            out(i).nSpk(per) = length(tmp_res);            
        end
        out(i).cids = clu_idx(i);        
    end
    for i = 1:length(clu_idx)
        % coefficient of variation for the waves' peak datum across periods
        x = min(out(i).wv,[],2);
        out(i).wvCV = nanstd(x)./abs(nanmean(x));
        % out(i).wvCV = mean(nanstd(out(i).wv(:,10:30))./abs(nanmean(out(i).wv(:,10:30))));
        a = findchangepts(x(~isnan(x)),'Statistic','linear');
        a = a - 1 + sum(isnan(x(1:a-1)));
        out(i).wvChPt = a-1; % changepoint
        if isempty(a); out(i).wvChPt = a-1; end
        [~,P,o,df] = ttest2(x(1:a-1),x(a:end)); %changepoint significance
        out(i).wvChPtPval = P; % changepoint
    end
    out(1).note = 'average cluster waveform as a function of trial. Includes delay periods too';

    
%% PLOT
    close all;
    ifPlot = 0;
    if ifPlot
        for i = 1:length(clu_idx)
            figure('pos',[1932         747        1596         800]);
            ax = subplot(1,5,[1 2 3]); 
            imAlpha=ones(size(out(i).wv));
            imAlpha(isnan(out(i).wv))=0;
            imagesc(out(i).wv,'AlphaData',imAlpha);
            set(gca,'color',0.5*[1 1 1]);
            colormap('summer'); 
            clb = colorbar; xlabel(clb,'ampl (\muV)');
            ylabel('periods (#)'); axis xy;
            ax_sz = get(ax,'pos');

            axes('pos',ax_sz);
            shadedErrorBar([],nanmean(out(i).wv,1),nanstd(out(i).wv),[],1); 
            hold on;
            plot([10 10],ylim,'--k')
            axis tight;
            axis off;


            subplot(1,5,[5]); 
            [haxes,hline1,hline2] = plotyy([1:1:length(periods)],out(i).nSpk./out(i).dur,...
                                           [1:1:length(periods)],out(i).dur);
            view([90 -90]); box off;
            set(haxes(:),'xlim',[1 length(periods)]); % ,'ycolor','k'
            set(haxes(2),'LineWidth',1); % ,'ycolor','r'
            ylabel(haxes(2),'duration (s)');
            ylabel(haxes(1),'rate (Hz)');
            set(hline2,'LineWidth',2,'Marker','o','linestyle','none');
            set(hline1,'LineWidth',2);% ,'Color',[0 0 1]
            xlabel('periods (#)'); 
            mtit(['cluster ID: ',(num2str(out(i).cluID ))]);
        end
    end



%% SAVE
    if ifSave
        cluWavePerTrial = out;
        save([pwd,'/cluster_waveform_by_trial.mat'],'cluWavePerTrial');
    end



