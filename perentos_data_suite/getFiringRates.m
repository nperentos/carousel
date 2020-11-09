function out = getFiringRates(fileBase,varargin)

% generates smoothed firing rate with window length win (default: 50ms)
% directly from kilosort output. 
% To process only a subset of units, provide cluster IDs 
% example call out = getFiringRates(fileBase,{'win',0.05,'units',[]})
% example call out = getFiringRates(fileBase,{'win',0.05,'units',cluIDs})
% example call out = getFiringRates(fileBase)


    ms = 1000;
    options = {'win',0.05,'units',[]};
    options = inputparser(varargin,options);
    win = options.win*ms; % 50 ms
    units = options.units;

% import the kilosort data
    processedPath = getfullpath(fileBase);
    cd(processedPath);
    sp = loadKSdir(pwd);
% load peripherals
    load peripherals; p = peripherals; clear peripherals;
% get LFP tScale for double checking
    tScale = carouselGetVar(fileBase,'tScale');

% ms time bins 
    L = ceil(sp.st(end)*ms); % round to the nearest ms - matches the lfp resolution
    bins = [0:win:L];
    nBins = length(bins)-1;

% all MUA and SUA
    if isempty(options.units)
        nC = length(unique(sp.cids(sp.cgs>0))); % n of cells
        cellList = unique(sp.cids(sp.cgs>0));
    else   
        nC = length(units);
        cellList = units;
    end

    rates = zeros(nBins,nC);

% convolution window     
    winLength = round(3*1000);% three second gaussian window
    gwin = gausswin(winLength,winLength/(win)); % standard deviation of 50 ms (win)

    for i = 1:nC % number of clusters
        res = sp.st(sp.clu == cellList(i)).*ms; % spike times for this cluster
        tmp = histcounts(res,bins);
        rates(:,i) = conv(tmp,gwin,'same')+eps;    
    end

%% plot 'report'
    close all;
    t = (bins +diff(bins(1:2)))/1000;
    figure; 
    poi = [618 1687 2299 4103];
    poi = repmat(poi,2,1);

    subplot(313);
    plot(tScale,normalize_array1(p.carouselSpeed)); axis tight; hold on;
    plot(tScale,-normalize_array1(p.runSpeed)); legend('carouselSpeed','runSpeed');
    plot(poi,repmat(ylim',1,length(poi)),'--k','linewidth',2);
    clm = flipud(colormap('gray'));
    subplot(3,1,[1 2]); 
    ir = zscore(rates)';
    imagesc(t,[],ir); hold on;
    plot(poi,repmat(ylim',1,length(poi)),'--k','linewidth',2);

    colormap(clm);
    colormap('jet');cax = caxis; caxis([cax(1),prctile(ir(:),95)]);
    linkaxes(get(gcf,'children'),'x');

%% OUTPUTS  
    out.r = rates;
    out.t = t;
