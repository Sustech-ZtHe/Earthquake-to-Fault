function AzimuthStat(Azimuth,maxrlim)
    polarhistogram(Azimuth.*pi/180, 36,...
        'LineWidth',1.5,'FaceColor',[220 20 60]./255);;
    ax=gca;
    ax.ThetaDir = 'clockwise'; % 顺时针方向
    ax.ThetaZeroLocation = 'top'; % 让0度位于顶部
    ax.GridLineStyle="--";
    ax.LineWidth=2;
    ax.GridColor=[0 0 0];
    ax.FontSize=12;
    ax.RAxis.Visible="on";
    if nargin==2
        rlim([0 maxrlim]);                  
    elseif nargin==1
        uprlim=max(histcounts(Azimuth.*pi/180,18));
        ax.RTick = 0:ceil(uprlim/5):uprlim; 
        rlim([0 uprlim]); 
    end
end
