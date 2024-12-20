function AzimuthStat(Azimuth)
    polarhistogram(Azimuth.*pi/180, 36,...
        'LineWidth',1.5,'FaceColor',[220 20 60]./255);
    uprlim=ceil(1.2*max(histcounts(Azimuth.*pi/180,36)));
    % thetaticklabels({'90', '60', '30', '0', '-30', '-60', '-90',...
    %     '-120', '-150','180', '150', '120'});
    ax=gca;
    ax.ThetaDir = 'clockwise'; % 顺时针方向
    ax.ThetaZeroLocation = 'top'; % 让0度位于顶部
    ax.GridLineStyle="--";ax.LineWidth=2;ax.GridColor=[0 0 0];ax.FontSize=12;ax.RAxis.Visible="off";
    rlim(ax,[0 uprlim]);
end
