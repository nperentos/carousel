function [nose] = importDLC(fileBase, startRow, endRow)

%% DISCOVER THE CSV FILE
    fullPath = fullfile(getFullPath(fileBase),'video');
    test = (dir(fullfile(fullPath,'top*.csv')));
    
    if size(test,1) ~=1
        error('no csv or more than one found - cannot proceed');        
    else
        filename = [test.folder,'/',test.name];
    end


%IMPORTFILE Import numeric data from a text file as column vectors.
%   [BODYPARTS1,NOSEBASETOP,NOSEBASETOP1,NOSEBASETOP2,NOSETIP,NOSETIP1,NOSETIP2,NOSEBASEBOTTOM,NOSEBASEBOTTOM1,NOSEBASEBOTTOM2]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [BODYPARTS1,NOSEBASETOP,NOSEBASETOP1,NOSEBASETOP2,NOSETIP,NOSETIP1,NOSETIP2,NOSEBASEBOTTOM,NOSEBASEBOTTOM1,NOSEBASEBOTTOM2]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [bodyparts1,nosebasetop,nosebasetop1,nosebasetop2,nosetip,nosetip1,nosetip2,nosebasebottom,nosebasebottom1,nosebasebottom2] = importfile('top_cam_0_date_2019_06_05_time_20_13_39_v001DeepCut_resnet50_NoseDLCDec2shuffle1_950000.csv',4, 153970);
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
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
% For more information, see the TEXTSCAN documentation.
    formatSpec = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';


%% Open the text file.
    fileID = fopen(filename,'r');


%% Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end

    
%% Close the text file.
    fclose(fileID);


%% Post processing for unimportable data.
    % No unimportable data rules were applied during the import, so no post
    % processing code is included. To generate code which works for
    % unimportable data, select unimportable cells in a file and regenerate the
    % script.

    
%% Allocate imported array to column variable names
    %T.n = dataArray{:, 1};
    % T.X = [dataArray{:, 2}, dataArray{:, 5}, dataArray{:, 8}]';
    % T.Y = [dataArray{:, 3}, dataArray{:, 6}, dataArray{:, 9}]';
    % T.P = [dataArray{:, 4}, dataArray{:, 7}, dataArray{:, 10}]';
    % original T
    varNames = {'nose_X1','nose_X2','nose_X3','nose_Y1','nose_Y2','nose_Y3','nose_p1','nose_p2','nose_p3'};
    T = zeros(size(dataArray{:, 2},1),9);
    T(:,1) = dataArray{:, 2};
    T(:,2) = dataArray{:, 5};
    T(:,3) = dataArray{:, 8};
    T(:,4) = dataArray{:, 3};
    T(:,5) = dataArray{:, 6};
    T(:,6) = dataArray{:, 9};
    T(:,7) = dataArray{:, 4};
    T(:,8) = dataArray{:, 7};
    T(:,9) = dataArray{:, 10};

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

% a per nose point deviation rather than COM deviation? It might be more robust
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
    mtit('nose position, X1,X2, X3, Y1 Y2 and Y3');
    print(fullfile(getfullpath(fileBase),'nose_DLC_tracking.jpg'),'-djpeg');


%% GET PRINCIPAL COMPONENTS
    rmpath(genpath('/storage2/perentos/code/thirdParty/drtoolbox'));
    [coeff, score, latent] = pca(Tint(:,1:6));

%% EXTRACT A VECTOR OF NOSE ORIENTATION    
for i = 1:2 % process both raw and interpolated traces
    clear L phi;
    if i == 2; v = T; elseif i == 1; v = Tint; end
    %  nose base midpoint confirmed as long as labels 1 and 3 are indeed the nose bases which should always be the xase
    m = [(v(:,1)+v(:,3))/2,  (v(:,4)+v(:,6))/2]; % nose ase midpoint
    apx = [v(:,2),v(:,5)]; % nose apex
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
        nose.vector.raw.phi = phi;
        nose.vector.raw.L   = L;
        nose.vector.raw.baseMidPoint   = m;
        nose.vector.raw.apex = apx;
        nose.vector.raw.polarDistr = Out;
        nose.vector.raw.polarXBins = XBins;
        nose.vector.raw.polarYBins = YBins;
    elseif i == 1
        nose.vector.int.phi = phi;
        nose.vector.int.L   = L;
        nose.vector.int.baseMidPoint = m;
        nose.vector.int.apex   = apx;        
        nose.vector.int.polarDistr = Out;
        nose.vector.int.polarXBins = XBins;
        nose.vector.int.polarYBins = YBins;
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
    title(['frame ',int2str(oo),',  ',num2str(rad2deg(nose.vector.raw.phi(oo))),char(176)]);
    view(-90,90);
    %pause; clf;
end


subplot(3,4,[3 4]);%subplot(3,4,[3 4 7 8]);
[h,c] = polarPcolor(nose.vector.raw.polarYBins,nose.vector.raw.polarXBins,nose.vector.raw.polarDistr','Nspokes',4,'Ncircles',2,'Rscale','log');
view(180, -90);
c.YLabel.String = 'raw signal';
subplot(3,4,[7 8]);
[h,c] = polarPcolor(nose.vector.int.polarYBins,nose.vector.int.polarXBins,nose.vector.int.polarDistr','Nspokes',4,'Ncircles',2,'Rscale','log');
view(180, -90)
c.YLabel.String = 'interp. signal';

subplot(3,4,[9 10]);
histogram(rad2deg(nose.vector.raw.phi),[-60:2:60]);%[N1,X1] = 
set ( gca, 'xdir', 'reverse' );
hold on;
histogram(rad2deg(nose.vector.int.phi),[-60:2:60]); %[N2,X2] = 
set ( gca, 'xdir', 'reverse' )
xlabel(['nose direction (',char(176),')']); 
ylabel(['count (',char(35),')']); 
legend('raw','interpolated');
axis tight;hold on;

subplot(3,4,[11 12]);
histogram(nose.vector.raw.L,[0:0.5:40]); hold on;
histogram(nose.vector.int.L,[0:0.5:40]);
xlabel(['nose length (pixels)']); ylabel(['count (',char(35),')']); axis tight;
xlim([0 40]);legend('raw','interpolated');
ForAllLabels('fontsize', 10, 'fontweight','normal');
ForAllSubplots('set(gca,''TickDir'',''out'',''box'',''off'',''fontsize'',10)');
print(fullfile(getfullpath(fileBase),'nose_angle.jpg'),'-djpeg');

%% ASSEMBLE REMAINING VARIABLES
    nose.dataRaw = T;
    nose.dataInterp = Tint;
    nose.dataMr= Tmr;
    nose.varNames = varNames;
    nose.otl1 = otl1;
    nose.otl2 = otl2;
    nose.otl = otl;
    note.px = px;
    nose.PCs.coeff = coeff;
    nose.PCs.score = score;
    nose.PCs.latent= latent;  
    nose.note = {'dataMr: mean corrected','dataIntert:interpolated',...
                 'otl1:outliers based on DLC posterior',...
                 'otl2:outliers based on center of mass',...
                 'otl:union of otl1 & otl2',...
                 'px:center of mass displacement allowed (in pixels)',...
                 'PCs:principal compoents'};

  

%% SAVE
    disp('saving nose data...');
    save(fullfile(getfullpath(fileBase),['noseTracking.mat']),'nose');  
    disp(['Nose tracking data saved in ', fullfile(getfullpath(fileBase)),'noseTracking.mat']);
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



















