function getPlaceTuninng(filebase,flag,pth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getPlaceTuninng(filebase,flag)
% generates tuning curves for all SUs and MUAs as generated from Kilosort 2
% and curated in phy.
% The script is specific to the carousel setup
% NPerentos - 26/06/2019
% -----INPUTS-----
% filebase: main folder with data e.g. filebase = 'NP25_2019-06-02_merged';
% flag:     defines the task used 
%           flag = 1 for continuous unidirectional
%           flag = 2 for cued direction
%           flag = 3 for passive
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
if nargin < 2
    display('no task flag supplied - assuming continuous uni-direction');
    flag = 1;
    display('no data path provided assuming default:: ''/storage2/perentos/data/recordings/''');
    processedPath = getfullpath(filebase);
end

if nargin == 3
    processedPath = getfullpath(filebase,pth);
end
% bring in peripherals and spikes

cd(processedPath)
disp 'loading peripherals'; 
load('peripherals.mat');
cd KS;
sp = loadKSdir(pwd);
tit = {'MUA','SU'}; 
display('in kilosort neurons labelled with a ''1'' are multiunit clusters. Those labeled with 2 are putative single units');

%% UNIDIRECTIONAL PROTOCOL
if flag == 1
    for U = 1:2 % unit type (SU and MUA)
        
        if U == 1; display('processing multi unit clusters...'); end
        if U == 2; display('processing putative single unit clusters...'); end

        figure;
        gU = sp.cids(sp.cgs == U);    
        p = numSubplots(length(gU));

        % sub select data where carousel or animal is moving
        idxMv = find(carouselSpeed>1); % idxMv = find(runSpeed>1);
        tmpLoc = location(idxMv);
        locEdges = linspace(min(tmpLoc),max(tmpLoc),360/5); % 5 degree bins or 2*pi*25/(360/5)=2cm
        positionInd = discretize(tmpLoc,locEdges);
        tmptScale = tScale(idxMv);
        % cycle through cells
        for i = 1:length(gU)   
            idxSp = find(sp.clu == gU(i));
            spT = round(sp.st(idxSp).*1000);
            res = zeros(size(tScale));
            res(spT) = 1; % spikes array as a series of 1s and 0s for the whole tScale
            tmpSpikes = res(idxMv); % keep spikes for movement periods

            % average cell tuning for carousel movement (or run?) periods only
            meanCellTuning = accumarray(positionInd',tmpSpikes,[numel(locEdges)-1,1],@mean); 

            % now cell tuning split in trials    
            clear cellTuning;    
            for j = 1:length(rewards.atStart)      
                % find the spikes that fall between carousel start points (per trial split)
                if j == 1
                    from = 1;
                else
                    from = find(tmptScale >= tScale( rewards.atStart(j))  ,1);
                end
                if j == length(rewards.atStart)
                    to = tmptScale(end);
                else
                    to = find(tmptScale >  tScale( rewards.atStart(j+1)),1) - 1;
                end                
                cellTuning(:,j) = accumarray(positionInd(from:to)',tmpSpikes(from:to),[numel(locEdges)-1,1],@sum); 
            end


            % plot
            subplot(p(1),p(2),i);
            imagesc(locEdges,1:1:size(cellTuning,2),cellTuning');  axis xy;
            clm = flipud(colormap('gray'));
            hold on;
            colormap(clm);
            hold on;    
            plot(locEdges(1:end-1)+(locEdges(2)-locEdges(1)),...
                smooth((meanCellTuning./max(meanCellTuning)).*length(rewards.atStart),3));
            title(num2str(gU(i)));
        end
        mtit(tit{U});
    end
end


%% CUED DIRECTION
rewards.atStart = round(rewards.atStart); % /15.764
rewards.atACW = round(rewards.atACW); % /15.764
[idxLR,idxLR]=sort(direction);
if flag == 2
    for U = 1:2 % unit type (SU and MUA)
        %special case for IIT probe with autozeroing issue

        % is the direction variable available? % (should be if the protocol was bidirectional)        
        if ~exist('direction','var')            
            error('direction variable is missing - cannot proceed. Check if truly the  cued direction behavioral protocol was used in this session');
        end
        
        if U == 1; display('processing multi unit clusters...'); end
        if U == 2; display('processing putative single unit clusters...'); end

        figure;
        gU = sp.cids(sp.cgs == U);    
        p = numSubplots(length(gU));

        % sub select data where carousel or animal is moving
        idxMv = find(carouselSpeed>1); % idxMv = find(runSpeed>1);
        tmpLoc = location(idxMv);
        locEdges = linspace(min(tmpLoc),max(tmpLoc),360/5); % 5 degree bins or 2*pi*25/(360/5)=2cm
        positionInd = discretize(tmpLoc,locEdges);
        tmptScale = tScale(idxMv);

        % cycle through cells
        for i = 1:length(gU)   
            idxSp = find(sp.clu == gU(i));
            spT = round(sp.st(idxSp).*1000);
            res = zeros(size(tScale));
            res(spT) = 1; % spikes array as a series of 1s and 0s for the whole tScale
            tmpSpikes = res(idxMv); % keep spikes for movement periods

            % average cell tuning for carousel movement (or run?) periods only
            meanCellTuning = accumarray(positionInd',tmpSpikes,[numel(locEdges)-1,1],@mean); 

            % now cell tuning split in trials    
            clear cellTuning;    
            for j = 1:length(rewards.atStart)      
                % find the spikes that fall between carousel start points (per trial split)
                if j == 1
                    from = 1;
                else
                    from = find(tmptScale >= tScale( rewards.atStart(j))  ,1);
                end
                if j == length(rewards.atStart)
                    to = tmptScale(end);
                elseif j>1 & j~= length(rewards.atStart)
                    to = find(tmptScale >  tScale( rewards.atStart(j+1)),1) - 1;
                elseif j == 1
                    to = find(tmptScale >  tScale( rewards.atStart(j)),1);
                end
                display(['from ', num2str(from), ' to ',num2str(to)]);
                cellTuning(:,j) = accumarray(positionInd(from:to)',tmpSpikes(from:to),[numel(locEdges)-1,1],@sum);
                if strcmp(direction(j),'L')
                    cellTuning(:,j) = flipud(cellTuning(:,j));
                end
                
            end


            % plot
%             figure(10+10*U-1); 
            
            subplot(p(1),p(2),i);
%             imagesc( locEdges(2:end), [], cellTuning(:,strfind(direction,'R'))');  axis xy; hold on;
%             clm = flipud(colormap('gray'));
%             hold on;
%             colormap(clm);
%             figure(10+10*U+1);
%             subplot(p(1),p(2),i);
%             imagesc(locEdges(2:end), [], cellTuning(:,strfind(direction,'L'))');

            imagesc(locEdges(2:end), [], cellTuning(:,idxLR)');xlim([0 4400]);
            clm = flipud(colormap('gray'));
            hold on;
            colormap(clm);
            hold on;    
            plot(locEdges(1:end-1)+(locEdges(2)-locEdges(1)),...
                smooth((meanCellTuning./max(meanCellTuning)).*length(rewards.atStart),3));
            title(num2str(gU(i)));
        end
        mtit(tit{U});
    end
end


%% PASSIVE ROTATIONS 
if flag == 3
   display 'not implemented yet'; 
end
