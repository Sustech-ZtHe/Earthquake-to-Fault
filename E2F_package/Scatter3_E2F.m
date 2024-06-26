function Scatter3_E2F(overevent,points,TOTAL,name,style,l)
% overevent=hypotout
% points=f_sort0
% TOTAL=NEvents
% name=['Unaccept']
% style=1
%style==0 [0 0.2 0.6] power(exp(points(:, 4)),1/5)+1
if nargin==6 && style==0 
    if ismember(points(:,1)<0,1)
        scatter3(abs(points(:, 1)), points(:, 2), points(:, 3), 5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','Xdir','reverse');box on;grid off;
    else
        scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse');box on;grid off;
    end
ax.XMinorTick="on";ax.YMinorTick="on";ax.ZMinorTick="on";
ax.TickLength=[0.02 0.1];grid minor;
pos = get(ax, 'Position');view(0,90);
annotation('textbox', [pos(1)+pos(3)*0.45, pos(2)+pos(4)*0.1, 0.6*pos(3), 0.1*pos(4)], 'LineStyle','none',...
            'String', ['Fault lines=',num2str(length(unique(points(:,14)))),' Events=',num2str(size(points,1)),' Remaining=',num2str(TOTAL-size(points,1)) ], ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black'); 
annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 1*pos(3), 0.85*pos(4)], 'LineStyle','none',...
            'String', name, ...
            'FontAngle','normal','FontSize', 15, 'Color', 'black');ax.FontSize=13;
% Toc2me 0.45
% Grant/Okamc 0.45
% Turkey 
end

% power(exp(points(:, 4)),1/5)

if nargin==5 && style==1 
    if ismember(points(:,1)<0,1)
        scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse','Xdir','reverse');box on;grid off;
        xticks([round(abs(max(points(:,1))),1):round((abs(min(points(:,1)))-abs(max(points(:,1))))/2,1):round(abs(min(points(:,1))),1)])
    else
        scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
        ax=gca;set(ax,'Zdir','reverse');box on;grid off;
    end
    ax.XMinorTick="on";ax.YMinorTick="on";ax.ZMinorTick="on";
    ax.TickLength=[0.02 0.1];grid minor;ax.FontSize=15;
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
        if ismember(points(:,1)<0,1)
%             f=convhull(points(:, 1:3))
%             patch('vertices',[abs(points(:, 1)), points(:, 2), points(:, 3)], ...
%                 'faces',f,'facecolor','red','facealpha',0.3,'linestyle','none')
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),5, points(:,4:6),"filled");
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse');box on;grid off;
        else
%             f=convhull(points(:, 1:3))
%             patch('vertices',points(:, 1:3), ...
%                 'faces',f,'facecolor','red','facealpha',0.3,'linestyle','none')
            scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,4:6),"filled");
            ax=gca;set(ax,'Zdir','reverse');box on;grid off;
        end
    else
        if ismember(points(:,1)<0,1)
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse');box on;grid off;
        else
            scatter3(points(:, 1), points(:, 2), points(:, 3),5, points(:,11:13),"filled");
            ax=gca;set(ax,'Zdir','reverse');box on;grid off;
        end
    end
    ax.XMinorTick="on";ax.YMinorTick="on";ax.ZMinorTick="on";
    ax.TickLength=[0.02 0.1];grid minor;ax.FontSize=10;
end

d1=(max(overevent(:,9))-min(overevent(:,9)))/50;
d2=(max(overevent(:,8))-min(overevent(:,8)))/50;
if ismember(points(:,1)<0,1)
    xlim([abs(max(overevent(:,9))+d1) abs(min(overevent(:,9))-d1)]);
    ylim([min(overevent(:,8))-d2 max(overevent(:,8))+d2]); 
else
    xlim([min(overevent(:,9))-d1 max(overevent(:,9))+d1]);
    ylim([min(overevent(:,8))-d2 max(overevent(:,8))+d2]);
end


if nargin==6 || nargin==5
    if ismember(points(:,1)<0,1)
        ax.XTickLabel=strcat(ax.XTickLabel, '°W')
    else
        ax.XTickLabel=strcat(ax.XTickLabel, '°E')
    end
    ax.YTickLabel=strcat(ax.YTickLabel, '°N');
end


% if nargin==6 && (l==2 || l==3 || l==4)
%     set(ax,'XTickLabel',[],'YTickLabel',[]);
% end
% if nargin==6 && (l==6 || l==7 || l==8)
%     set(ax,'YTickLabel',[]);
% end









