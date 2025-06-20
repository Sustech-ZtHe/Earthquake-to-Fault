function Scatter3_E2F(hypo,points,name,style)
if ~isempty(points)
    psize=8;
    if nargin==4 && style==0 
        if ismember(points(:,1)<0,1)
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3), psize, points(:,11:13),"filled");view(0,90);
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
            box on;grid off;grid minor
        else
            scatter3(points(:, 1), points(:, 2), points(:, 3),psize, points(:,11:13),"filled");view(0,90);
            ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
            box on;grid off;grid minor
        end
        pos = get(ax, 'Position');
        annotation('textbox', [pos(1)+pos(3)*0.4, pos(2)+pos(4)*0.05, 0.8*pos(3), 0.1*pos(4)], 'LineStyle','none',...
                    'String', ['Fault candidates=',num2str(length(unique(points(:,14)))),newline,'Clustered events=',num2str(size(points,1)),newline,'Remaining events=',num2str(size(hypo,1)-size(points,1)) ], ...
                    'FontAngle','normal','FontSize', 10, 'Color', 'black'); 
        annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 1*pos(3), 0.85*pos(4)], 'LineStyle','none',...
                'String', name, ...
                'FontAngle','normal','FontSize', 12, 'Color', 'black');
    end
    if nargin==4 && style==1 
        if min(points(:,1))<0
            scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),psize, points(:,11:13),"filled");view(0,90)
            ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',15,'color',[0.9 0.9 0.9]);
            box on;grid off;grid off
        else
            scatter3(points(:, 1), points(:, 2), points(:, 3),psize, points(:,11:13),"filled");view(0,90)
            ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',15,'color',[0.9 0.9 0.9]);
            box on;grid off;grid off
        end
        pos = get(ax, 'Position');
        if name~=char(10)
        annotation('textbox', [pos(1)+pos(3)*0.01, pos(2)+pos(4)*0.08, 1.5*pos(3), 0.01*pos(4)], 'LineStyle','none',...
                'String', ['Fault candidates=',num2str(length(unique(points(:,14)))),newline,'Events used in Hough transform=',num2str(size(points,1))], ...
                'FontSize', 15, 'Color', 'black'); 
        annotation('textbox', [pos(1)+pos(3)*0.01, pos(2)+pos(4)*0.14, 0.4*pos(3), 0.85*pos(4)], 'LineStyle','none',...
                    'String', name, ...
                    'FontSize', 20, 'Color', 'black');
        end
    end
    if nargin==2 
        if size(points,2)==6
            if min(points(:,1))<0
                scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),psize, points(:,4:6),"filled");view(0,90);
                ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
                box on;grid off;grid minor
            else
                scatter3(points(:, 1), points(:, 2), points(:, 3),psize, points(:,4:6),"filled");view(0,90);
                ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
                box on;grid off;grid minor
            end
        else
            if min(points(:,1))<0
                scatter3(abs(points(:, 1)), points(:, 2), points(:, 3),psize, points(:,11:13),"filled");view(0,90);
                ax=gca;set(ax,'Zdir','reverse','Xdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
                box on;grid off;grid minor
            else
                scatter3(points(:, 1), points(:, 2), points(:, 3),psize, points(:,11:13),"filled");view(0,90);
                ax=gca;set(ax,'Zdir','reverse','YMinorTick','on','XMinorTick','on','ZMinorTick','on','FontSize',12,'color',[0.9 0.9 0.9]);
                box on;grid off;grid minor
            end
        end
    end
    d1=(max(hypo(:,9))-min(hypo(:,9)))/50;
    d2=(max(hypo(:,8))-min(hypo(:,8)))/50;
    if min(points(:,1))<0
        xlim([abs(max(hypo(:,9))+d1) abs(min(hypo(:,9))-d1)]);
        ax.XTickLabel = string(ax.XTick) + "°W";
    else
        xlim([min(hypo(:,9))-d1 max(hypo(:,9))+d1]);
        ax.XTickLabel = string(ax.XTick) + "°E";
    end
    if min(points(:,2))>0
        ylim([min(hypo(:,8))-d2 max(hypo(:,8))+d2]);
        ax.YTickLabel = string(ax.YTick) + "°N";
    else
        ylim([abs(max(hypo(:,8))+d2) abs(min(hypo(:,8))-d2)]);
        ax.YTickLabel = string(ax.YTick) + "°S";
    end
end










