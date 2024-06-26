function AzimuthStat(Azimuth)
    polarhistogram(Azimuth.*pi/180, 36,...
        'LineWidth',1,'FaceColor',[220 20 60]./255);
    uprlim=round(1.2*(max(histcounts(Azimuth.*pi/180,36))));
    thetaticklabels({'E', 'ENE', 'NE', 'NNE', 'N', 'NNW', 'NW',...
        'WNW', 'W', 'WSW', 'SW', 'SSW','S','SSE','SE','ESE'})
    ax=gca;
    ax.GridLineStyle="--";ax.LineWidth=2;ax.GridColor=[0 0 0];ax.FontSize=15;ax.RAxis.Visible="off"
    rlim([0 uprlim]);
end