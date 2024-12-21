function ElevationStat(DipA)
    polarhistogram(DipA.*pi/180, (0:10:90).*pi/180,...
        'LineWidth',1.5,'FaceColor',[70 130 180]./255);
    thetalim([0,90]);
    uprlim=ceil(1.2*max(histcounts(DipA.*pi/180,(0:10:90).*pi/180)));
    % thetaticklabels({'0','30','60','90'});
    ax=gca;
    % ax.ThetaDir = 'clockwise'; % 顺时针方向
    % ax.ThetaZeroLocation = 'top'; % 让0度位于顶部
    ax.GridLineStyle="--";ax.LineWidth=2;ax.GridColor=[0 0 0];ax.FontSize=15;ax.RAxis.Visible="off";
    rlim([0 uprlim]);
end