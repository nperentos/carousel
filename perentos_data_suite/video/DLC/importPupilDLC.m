function [pupilDLC] = importPupilDLC(fileBase, startRow, endRow)

%% DISCOVER THE CSV FILE
    fullPath = fullfile(getFullPath(fileBase),'video');
    test = (dir(fullfile(fullPath,'*pupil*.csv')));
    
    if size(test,1) ~=1
        error('no csv or more than one found - cannot proceed');        
    else
        filename = [test.folder,'/',test.name];
    end

    load(fullfile(getfullpath(fileBase),'faceROIs.mat'));
    eye = ROI.pts(find(strcmp(ROI.names,'eye')),:);
    
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [BODYPARTS1,pupilBASETOP,pupilBASETOP1,pupilBASETOP2,pupilTIP,pupilTIP1,pupilTIP2,pupilBASEBOTTOM,pupilBASEBOTTOM1,pupilBASEBOTTOM2]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BODYPARTS1,pupilBASETOP,pupilBASETOP1,pupilBASETOP2,pupilTIP,pupilTIP1,pupilTIP2,pupilBASEBOTTOM,pupilBASEBOTTOM1,pupilBASEBOTTOM2]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [bodyparts1,pupilbasetop,pupilbasetop1,pupilbasetop2,pupiltip,pupiltip1,pupiltip2,pupilbasebottom,pupilbasebottom1,pupilbasebottom2] = importfile('top_cam_0_date_2019_06_05_time_20_13_39_v001DeepCut_resnet50_pupilDLCDec2shuffle1_950000.csv',4, 153970);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/12/11 10:02:42

%% Initialize variables.
    delimiter = ',';
    if nargin<=2
        startRow = 4;
        endRow = inf;
    end


%% Format for each line of text:
%   column2: double (%f)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
%	column7: double (%f)
%   column8: double (%f)
%	column9: double (%f)
%   column10: double (%f)
%	column11: double (%f)
%   column12: double (%f)
%	column13: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% GET ELLIPSE AREA AND CENTER
% X = pupilDLC(:,  [1 10 4 7]);
% Y = pupilDLC(:,1+[1 10 4 7]);
% for i = 1:5%size(X,1)
%     A(i,:) = fitellipse(X(i,:),Y(i,:)); % A = arrayfun(@(X,Y) fitellipse(X,Y),X(1,:),Y(1,:),'UniformOutput',false);
% end

%% Create output variable
pup = [dataArray{1:end-1}];
clearvars filename delimiter startRow formatSpec fileID dataArray ans;


likelihoods = pup(:,[3:3:24]);
otl = find(sum(likelihoods,2) ~= 8);
X = pup(:,  [1:3:24]);
Y = pup(:,[2:3:24]);
refl = pup(:,[25:26]);
clear A
for i = 1:size(X,1)
    %A(i,:) = fitellipse(X(i,:),Y(i,:));
    A(i,:) = circfit(X(i,:),Y(i,:));
    %A(i,:) = EllipseDirectFit([X(i,:)',Y(i,:)']);
end
% A = [1 1 1 1 180].*A;


idx = find(sum(likelihoods,2) <3);
[~,re] = sort(rand(size(idx)));
idx = idx(re);
%idx = 1:100;
vidfle = dir('side*.avi'); figure; clf;
if length(vidfle) == 1
    obj = VideoReader(vidfle.name);
    for i = 1:100%length(idx)
        obj.CurrentTime = idx(i);
        fr = obj.readFrame;
        imshow(fr); hold on;
        % dlc labels
        plot(pupilDLC(idx(i),1),pupilDLC(idx(i),2),'or');
        plot(pupilDLC(idx(i),4),pupilDLC(idx(i),5),'og');
        plot(pupilDLC(idx(i),7),pupilDLC(idx(i),8),'ob');
        plot(pupilDLC(idx(i),10),pupilDLC(idx(i),11),'om');
        xlim([eye(1) eye(1)+eye(3)]);
        ylim([eye(2) eye(2)+eye(4)]);
        % fitted circle
        viscircles([A(i,2), A(i,3)],A(i,1));
%         % fitted ellipse
%         Returned = A(i,:);
%         % Draw the returned ellipse
%         t = linspace(0,pi*2);
%         x = Returned(3) * cos(t);
%         y = Returned(4) * sin(t);
%         nx = x*cos(Returned(5))-y*sin(Returned(5)) + Returned(1); 
%         ny = x*sin(Returned(5))+y*cos(Returned(5)) + Returned(2);
%         hold on
%         plot(nx,ny,'r-')
        pause;clf;
    end
end



%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
varNames = {'pupil_XL','pupil_YL','pupil_PL',...
            'pupil_XR','pupil_YR','pupil_PR',...
            'pupil_XU','pupil_YU','pupil_PU',...
            'pupil_XD','pupil_YD','pupil_PD'};
        
        
%% quick vis
figure;
scatter(pupilDLC(:,1),pupilDLC(:,2),'r');
hold on;
scatter(pupilDLC(:,4),pupilDLC(:,5),'g');
scatter(pupilDLC(:,7),pupilDLC(:,8),'b');
scatter(pupilDLC(:,10),pupilDLC(:,11),'m');
        
        
        
        
        
        
        
        
        
        
        
        
        
% make sure that the order of labels is as expected
    mus = mean(T(:,1:6));
    [~, i_apex_x]=max(mus(1:3));
    [~, i_base_right_y]=max(mus(4:6));
    [~, i_base_left_y]=min(mus(4:6));    
    seq = [i_base_left_y i_apex_x i_base_right_y];
    T = T(:,[seq seq+3 seq+6]);
    
% center of mass or per point movement tolerance away from the mean
    px = 10;


% find unreliable frames based on DLC posteriors and center of mass calc.
    tmp = min(T(:,7:9),[],2);
    otl1 = find(tmp<0.995);    

% find unreliable frames based center of mass of three points
    cx = median(vc(T(:,1:3))); cy = median(vc(T(:,4:6)));% center of mass
    cxi = mean(T(:,1:3),2); cyi = mean(T(:,4:6),2); % COM per timepoint    
    otl2 = find(abs(cxi-cx)>px|abs(cyi-cy)>px); 

% a per pupil point deviation rather than COM deviation? It might be more robust
%     cx = median(T(:,1:3)); cy = median(T(:,4:6));
%     otl2 = [];
%     for i = 1:3
%         otl2 = [otl2; find(abs(T(:,i)-cx(i))>px)];
%     end
%     for i = 1:3
%         otl2 = [otl2; find(abs(T(:,i+3)-cy(i))>px)];
%     end
%     otl2 = unique(otl2);



    
    % all detected outliers
    otl = union(otl1,otl2);
    disp(['detected ',num2str(round(length(otl)/6)),' possible outliers']);
% mean corrected
    Tmr = T;
    mns = mean(T,1);
    for i = 1:6
        Tmr(otl,i) = mns(i);
    end    
% interpolate the nans
    Tint = T;
    Tint(otl,1:6) = nan; Tnans = Tint;
    t = [1:size(Tint,1)]';
    nanx = isnan(Tint(:,i));
    Tint(nanx,1:6) = interp1(t(~nanx),Tint(~nanx,1:6),t(nanx),'linear');

% summary plot
    close all;
    figure('pos',[1 1 1500 1100]);    
    hx = tight_subplot(6,1,0.03);
    for i = 1:6  
        axes(hx(i));
        jplot((Tint(:,i)),'linewidth',2); hold on;yl=ylim;
        jplot((Tmr(:,i)),'linewidth',1.5);    
        jplot((T(:,i)),'linewidth',1);
        %jplot((Tnans(:,1)),'*','linewidth',1.5);
        ylim(yl);        
%     chH = get(gca,'Children'); %uistack(hx(i),'top')
%     set(gca,'Children',[chH(end);chH(1:end-1)]);
    end
    legend('interpolated','mean corrected','original');%,'missing vals');% 'X1 original',
    %
    linkaxes(hx,'x','y');
    % find the middle of the interpolated values series and center aroun it
    tt = otl(round(end/3));
    xlim([tt-1000,tt+1000]);
    %xlim([73474.1322751323          73590.5343915344]);
    %xlim([ 53586       55586]);
    xlim([1,60*30]);
    xlabel('time (frames)');
    mtit('pupil position, X1,X2, X3, Y1 Y2 and Y3');
    print(fullfile(getfullpath(fileBase),'pupil_DLC_tracking.jpg'),'-djpeg');


%% GET PRINCIPAL COMPONENTS
    rmpath(genpath('/storage2/perentos/code/thirdParty/drtoolbox'));
    [coeff, score, latent] = pca(Tint(:,1:6));

%% EXTRACT A VECTOR OF pupil ORIENTATION    
for i = 1:2 % process both raw and interpolated traces
    clear L phi;
    if i == 2; v = T; elseif i == 1; v = Tint; end
    %  pupil base midpoint confirmed as long as labels 1 and 3 are indeed the pupil bases which should always be the xase
    m = [(v(:,1)+v(:,3))/2,  (v(:,4)+v(:,6))/2]; % pupil ase midpoint
    apx = [v(:,2),v(:,5)]; % pupil apex
    %L = diag(pdist2(m,apx)); % eucledian length from base midpoint to apex
    % the above susceptible to memory limitations - break it in pieces    
    lll=9;
    rs = diff(fix(linspace(0,length(m),lll+1)));
    C = mat2cell(m, rs, 2);
    D = mat2cell(apx, rs, 2);
    for jj = 1:lll
        Lcell{jj} = diag(pdist2(C{jj},D{jj}))';
    end
    L = cell2mat(Lcell)';
    phi = acos((apx(:,1)-m(:,1))./L).*sign(m(:,2)-apx(:,2)); % eucledian angle of said vector wrt to the x-axis of the 

    [Out, XBins, YBins, Pos] = hist2([rad2deg(phi),L],100,100);
    XBins = XBins(2:end) - mean(diff(XBins))/2;
    YBins = YBins(2:end) - mean(diff(YBins))/2;   
    
    if i == 2
        pupil.vector.raw.phi = phi;
        pupil.vector.raw.L   = L;
        pupil.vector.raw.baseMidPoint   = m;
        pupil.vector.raw.apex = apx;
        pupil.vector.raw.polarDistr = Out;
        pupil.vector.raw.polarXBins = XBins;
        pupil.vector.raw.polarYBins = YBins;
    elseif i == 1
        pupil.vector.int.phi = phi;
        pupil.vector.int.L   = L;
        pupil.vector.int.baseMidPoint = m;
        pupil.vector.int.apex   = apx;        
        pupil.vector.int.polarDistr = Out;
        pupil.vector.int.polarXBins = XBins;
        pupil.vector.int.polarYBins = YBins;
    end

end


%% SUMMARY PLOT
close all;
cd(fullfile(getFullPath(fileBase),'video'));
lst = dir('top*.avi');
if length(lst) == 1
    obj = VideoReader(lst.name);
end

figure('pos',[10 10 1100 800]); subplot(3,4,[1 2 5 6] );
for oo = 9%1:100
    obj.CurrentTime = oo-1;
    fr = obj.readFrame;
    imshow(fr); hold on;
    nPx = 100;
    xlim([cx-nPx cx+nPx]);
    ylim([cy-nPx cy+nPx]);
    plot(T(oo,1),T(oo,4),'xr','markersize',14,'linewidth',2);
    plot(T(oo,2),T(oo,5),'xg','markersize',14,'linewidth',2);
    plot(T(oo,3),T(oo,6),'xb','markersize',14,'linewidth',2);
    plot([m(oo,1) T(oo,2)],[m(oo,2) v(oo,5)],'w');
    plot(xlim,[m(oo,2) m(oo,2)])
    title(['frame ',int2str(oo),',  ',num2str(rad2deg(pupil.vector.raw.phi(oo))),char(176)]);
    view(-90,90);
    %pause; clf;
end


subplot(3,4,[3 4]);%subplot(3,4,[3 4 7 8]);
[h,c] = polarPcolor(pupil.vector.raw.polarYBins,pupil.vector.raw.polarXBins,pupil.vector.raw.polarDistr','Nspokes',4,'Ncircles',2,'Rscale','log');
view(180, -90);
c.YLabel.String = 'raw signal';
subplot(3,4,[7 8]);
[h,c] = polarPcolor(pupil.vector.int.polarYBins,pupil.vector.int.polarXBins,pupil.vector.int.polarDistr','Nspokes',4,'Ncircles',2,'Rscale','log');
view(180, -90)
c.YLabel.String = 'interp. signal';

subplot(3,4,[9 10]);
histogram(rad2deg(pupil.vector.raw.phi),[-60:2:60]);%[N1,X1] = 
set ( gca, 'xdir', 'reverse' );
hold on;
histogram(rad2deg(pupil.vector.int.phi),[-60:2:60]); %[N2,X2] = 
set ( gca, 'xdir', 'reverse' )
xlabel(['pupil direction (',char(176),')']); 
ylabel(['count (',char(35),')']); 
legend('raw','interpolated');
axis tight;hold on;

subplot(3,4,[11 12]);
histogram(pupil.vector.raw.L,[0:0.5:40]); hold on;
histogram(pupil.vector.int.L,[0:0.5:40]);
xlabel(['pupil length (pixels)']); ylabel(['count (',char(35),')']); axis tight;
xlim([0 40]);legend('raw','interpolated');
ForAllLabels('fontsize', 10, 'fontweight','normal');
ForAllSubplots('set(gca,''TickDir'',''out'',''box'',''off'',''fontsize'',10)');
print(fullfile(getfullpath(fileBase),'pupil_angle.jpg'),'-djpeg');

%% ASSEMBLE REMAINING VARIABLES
    pupil.dataRaw = T;
    pupil.dataInterp = Tint;
    pupil.dataMr= Tmr;
    pupil.varNames = varNames;
    pupil.otl1 = otl1;
    pupil.otl2 = otl2;
    pupil.otl = otl;
    note.px = px;
    pupil.PCs.coeff = coeff;
    pupil.PCs.score = score;
    pupil.PCs.latent= latent;  
    pupil.note = {'dataMr: mean corrected','dataIntert:interpolated',...
                 'otl1:outliers based on DLC posterior',...
                 'otl2:outliers based on center of mass',...
                 'otl:union of otl1 & otl2',...
                 'px:center of mass displacement allowed (in pixels)',...
                 'PCs:principal compoents'};

  

%% SAVE
    disp('saving pupil data...');
    save(fullfile(getfullpath(fileBase),['pupilTracking.mat']),'pupil');  
    disp(['pupil tracking data saved in ', fullfile(getfullpath(fileBase)),'pupilTracking.mat']);
    disp('DONE!');  






%     for i = 1:6 
%         nanx = isnan(Tint(:,i));
%         Tint(nanx,i) = interp1(t(~nanx),Tint(~nanx,i),t(nanx),'pchip');
%     end










% [Out1, XBins, YBins, Pos] = hist2([rad2deg(phi(1:30000)),L(1:30000)],100,100);
% [Out2, XBins, YBins, Pos] = hist2([rad2deg(phi(30000:end)),L(30000:end)],100,100);
% 
% figure;
% subplot(121);
% polarPcolor(YBins',XBins',Out1);
% subplot(122);
% polarPcolor(YBins',XBins',Out2);





%% i copied this from N_scrap and I place it here temporarily till i return to this analysis
%% next
figure;

% for i = 1:length(pup)
%     hold off;
%     scatter(pup(i,[2:3:24]),pup(i,[3:3:24]));
%     hold on;
%     scatter(pup(i,[26]),pup(i,[27]),'r');
%     xlim([0 100]);ylim([0 100]); set(gca,'YDir','reverse');
%     drawnow;
% end


% data = pup(:,[[2:3:24],[3:3:24]]);
data = pup(:,[3:3:24]);
[sig,mu,mah,outliers] = robustcov(data);
idx1 = find(outliers == 1);

idx2 = find(sum(data,2) ~= 8);
idx = idx2;
%figure;
for i = 1:length(idx)
    hold off;
    scatter(pup(idx(i),[1:3:24]),pup(idx(i),[2:3:24]));
    %hold on;
    %scatter(pup(idx(i),[26]),pup(idx(i),[27]),'r');
    xlim([0 100]);ylim([0 100]); set(gca,'YDir','reverse');
    title(num2str(idx(i)));
    drawnow;
    % pause;
end
% this is an outlier frame and we will test how to find the goodness of fit
% to it using ellipse or circle fitting
otl = find(sum(likelihoods,2) ~= 8);
for i = 20:length(otl)    
    Xc = X(otl(i),:)';
    Yc = Y(otl(i),:)';
    clf;
    plot(Xc,Yc,'xr'); hold on; % the dlc data points
    % polygons COM and area as most robust estimation of area and centroid
    [ geom, iner, cpmo ] = polygeom(Xc,Yc); %area and centroid of said polygon
    plot(geom(2),geom(3),'ok');
    % fit circle (used for residuals)
    out = circfit(Xc,Yc); % fit a circle to the data points
    viscircles([out(2),out(3)],out(1)); 
    resid(i,1) = sum(abs(sqrt((Xc-out(2)).^2 + (Yc-out(3)).^2)-out(1))); % circular fit residual
    % fi ellipse (used for residuals)
    [params] = fitellipse(Xc,Yc);
    Returned = round(params.*[1 1 1 1 180]);
    t = linspace(0,pi*2,40);
    x = params(3) * cos(t);
    y = params(4) * sin(t);
    nx = x*cos(params(5))-y*sin(params(5)) + params(1); 
    ny = x*sin(params(5))+y*cos(params(5)) + params(2);
    % since we would need to figure out the rotation of the fitted ellipse and
    % the corresponding radius at that angle to be able to compute the residual
    % on an individual point's basis, we can take a shortcut and instead
    % compute all the distances of the each point in question to those define
    % by the ellipse(nx,ny). The minimum of this set will be the equivalent to
    % the residual as defined above for circfit
    resid = 0;
    for pt = 1:8
        resid_(pt) = min(sqrt((Xc(pt)-nx).^2 + (Yc(pt)-ny).^2));
    end
    resid(i,2) = sum(resid_);
    
    hold on;
    plot(nx,ny,'b-');
    
    text(geom(2),geom(3),num2str(resid(i,2)));
    title(num2str(sum(likelihoods(otl(i),:))/8));
    xlim([0 100]);ylim([0 100]);
    axis equal;
    drawnow; pause(0.5*(resid(i,2)/20));
end













