function [h] = stampFig(fileBase,h,other)

if nargin < 2
    h = gcf;
end
if nargin == 3
    if ~iscell(other)
        error('input variable other must be of type ''cell''');
    end
end

T= getCarouselDataBase;
i = find(strcmp(T.session,fileBase));
h = axes('pos',[0 0 1 1]);
text(0.01,0.98,[fileBase],'fontsize',9,'interpreter','none');
hold on;
text(0.01,0.95,['db_entry=',num2str(i)],'fontsize',9,'interpreter','none');

if exist('other') == 1
    for i = 1:length(other)
        text(0.01,0.92,other{i},'fontsize',9,'interpreter','none'); 
    end
end

xlim([0 1]); ylim([0 1]);
axis off;
