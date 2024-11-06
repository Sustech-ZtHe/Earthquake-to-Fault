function Scatter3_E2F(overevent,points,TOTAL,name,style,l)

% points=f_sort1
% overevent=hypotout
% TOTAL=NEvents

if nargin==6 && style==0 
    if ismember(points(:,1)<0,1)
        scatter3(abs(points(:, 1)), points(:, 2), points(:, 3), 5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
        ax.XTickLabel=strcat(ax.XTickLabel, '°W')
        box on;grid off;grid minor
    else
        scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
        ax.XTickLabel=strcat(ax.XTickLabel, '°E')
        box on;grid off;grid minor
    end
    ax.YTickLabel=strcat(ax.YTickLabel, '°N');
    pos = get(ax, 'Position');view(0,90);
    annotation('textbox', [pos(1)+pos(3)*0.45, pos(2)+pos(4)*0.1, 0.6*pos(3), 0.1*pos(4)], 'LineStyle','none',...
                'String', ['Fault lines=',num2str(length(unique(points(:,14)))),' Events=',num2str(size(points,1)),' Remaining=',num2str(TOTAL-size(points,1)) ], ...
                'FontAngle','normal','FontSize', 12, 'Color', 'black'); 
    annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 1*pos(3), 0.85*pos(4)], 'LineStyle','none',...
            'String', name, ...
            'FontAngle','normal','FontSize', 12, 'Color', 'black');ax.FontSize=13;
end
if nargin==5 && style==1 
    if min(points(:,1))<0
        scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),10, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
        ax.XTickLabel=strcat(ax.XTickLabel, '°W')
        box on;grid off;grid off
    else
        scatter3(points(:, 1), points(:, 2), points(:, 3),10, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
        ax.XTickLabel=strcat(ax.XTickLabel, '°E')
        box on;grid off;grid off
    end
    ax.YTickLabel=strcat(ax.YTickLabel, '°N');
    pos = get(ax, 'Position');view(0,90)
    if name~=char(10)
    annotation('textbox', [pos(1)+pos(3)*0.45, pos(2)+pos(4)*0.1, 0.8*pos(3), 0.01*pos(4)], 'LineStyle','none',...
            'String', ['Fault Lines=',num2str(length(unique(points(:,14)))),char(10),'Events=',num2str(size(points,1)),'(',num2str(round(size(points,1)/TOTAL,2)*100),'%)'], ...
            'FontSize', 15, 'Color', 'black'); 
    annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 0.8*pos(3), 0.85*pos(4)], 'LineStyle','none',...
                'String', name, ...
                'FontSize', 15, 'Color', 'black');
    end
end
if nargin==2 
    if size(points,2)==6
        if min(points(:,1))<0
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),5, points(:,4:6),"filled");
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
            ax.XTickLabel=strcat(ax.XTickLabel, '°W')
            box on;grid off;grid minor
        else
            scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,4:6),"filled");
            ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
            ax.XTickLabel=strcat(ax.XTickLabel, '°E')
            box on;grid off;grid minor
        end
        ax.YTickLabel=strcat(ax.YTickLabel, '°N');
    else
        if min(points(:,1))<0
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
            ax.XTickLabel=strcat(ax.XTickLabel, '°W')
            box on;grid off;grid minor
        else
            scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
            ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12);
            ax.XTickLabel=strcat(ax.XTickLabel, '°E')
            box on;grid off;grid minor
        end
        ax.YTickLabel=strcat(ax.YTickLabel, '°N');
    end
end


d1=(max(overevent(:,9))-min(overevent(:,9)))/50;
d2=(max(overevent(:,8))-min(overevent(:,8)))/50;
if min(points(:,1))<0
    xlim([abs(max(overevent(:,9))+d1) abs(min(overevent(:,9))-d1)]);
    ylim([min(overevent(:,8))-d2 max(overevent(:,8))+d2]); 
else
    xlim([min(overevent(:,9))-d1 max(overevent(:,9))+d1]);
    ylim([min(overevent(:,8))-d2 max(overevent(:,8))+d2]);
end











