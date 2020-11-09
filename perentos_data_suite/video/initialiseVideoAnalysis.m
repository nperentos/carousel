function initialiseVideoAnalysis(fileBase)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function initialiseVideoAnalysis promps the user to define ROIs for the
% video analysis which are a prerequisite to all the video analysis except
% those that use DeepLabCut.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin == 0
    error('Cannot proceed without a fileBase - please provide');
end


%% look for videos
    
    fullPath = fullfile(getFullPath(fileBase),'video');
    
    % side view
    test1 = (dir(fullfile(fullPath,'side*.avi')));
    if ~size(test1,1)
        warning('there is no side*.avi video file .. skipping session');
        return;
    elseif size(test1,1)>1
        error('multiple side*.avi files were found - please investigate');
    elseif size(test1,1) == 1
        disp 'found a side*.avi file; continuing with data processing...'
        flag = [1];
        fileName = [test1.folder,'/',test1.name];
        obj1 = VideoReader(fileName);
        %[pth,fle,ext] = fileparts(fileName);
    end
    
    % top view
    test2 = (dir(fullfile(fullPath,'top*.avi')));
    if ~size(test1,1)
        warning('there is no top*.avi file .. skipping session');
        return;
    elseif size(test2,1)>1
        error('multiple top*.avi files were found - please investigate');
    elseif size(test2,1) == 1
        disp 'found a top*.avi file; continuing with data processing...'
        flag = [flag 1];
        fileName = [test2.folder,'/',test2.name];
        obj2 = VideoReader(fileName);
        %[pth,fle,ext] = fileparts(fileName);
    end
    
    
%% initialise the ROIs for later processing

    if flag(1) 
        trackPupilMorpho(fileBase,1,0);
        subFrameSVD(fileBase);
    end
    if flag(2)
        trackWhiskersGabor(fileBase);
    end
    
%% a figure showing the ROIs for summary visualisation
    % load the ROIs
    figure;

    face = load(fullfile(getfullpath(fileBase),'faceROIs'));
    whisker = load(fullfile(getfullpath(fileBase),'whiskerROIs'));

    % plot first frame for side
    subplot(121);
    imagesc(rgb2gray(im2double(obj1.readFrame))); hold on; axis equal; colormap gray;
    % add the side ROIs:
    for  i = 1:length(face.ROI.names)
        rectangle('Position',face.ROI.pts(i,:),'EdgeColor','r','LineWidth',2,'LineStyle','--');
    end



    % plot first frame for top
    subplot(122);
    imagesc(rgb2gray(im2double(obj2.readFrame))); hold on; axis equal; colormap gray;
    % add the top ROIs:
    for  i = 1:length(whisker.ROI.names)
        rectangle('Position',whisker.ROI.pts(i,:),'EdgeColor','r','LineWidth',1.5,'LineStyle','--');end

    print(gcf,fullfile(getfullpath(fileBase),'ROI_definitions.jpg'),'-djpeg','-r600');
