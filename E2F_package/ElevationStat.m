function ElevationStat(DipA)
    polarhistogram(DipA.*pi/180, (0:10:90).*pi/180,...
        'LineWidth',1,'FaceColor',[70 130 180]./255);
    thetalim([0,90])
    uprlim=max(histcounts(DipA.*pi/180,(0:10:90).*pi/180));
    thetaticklabels({'E', 'ENE', 'NNE', 'N', 'NNW',...
        'WNW', 'W', 'WSW', 'SSW','S','SSE','ESE'})
    ax=gca;
    ax.GridLineStyle="--";ax.LineWidth=2;ax.GridColor=[0 0 0];ax.FontSize=15;
    rlim([0 uprlim]);
end