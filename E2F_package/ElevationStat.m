function ElevationStat(DipA,maxrlim)
    polarhistogram(DipA.*pi/180, (0:10:90).*pi/180, ...
        'LineWidth', 1.5, 'FaceColor', [70 130 180]./255);
    thetalim([0, 90]);
    ax = gca;
    
    % 设置极坐标轴属性
    ax.GridLineStyle = "--";
    ax.LineWidth = 2;
    ax.GridColor = [0 0 0];
    ax.FontSize = 15;
    
    % 关键修正：启用半径轴并设置刻度
    ax.RAxis.Visible = "on";           % 开启半径轴
    if nargin==2
        rlim([0 maxrlim]);                  
    elseif nargin==1
        uprlim=max(histcounts(DipA.*pi/180, (0:10:90).*pi/180))+5;
        ax.RTick = 0:ceil(uprlim/5):uprlim; 
        rlim([0 uprlim]); 
    end
end