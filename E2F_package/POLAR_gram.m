function POLAR_gram(angle1,style,angle2)
if nargin==2 && style==0
polarhistogram(angle1.*pi/180, (-90:10:90).*pi/180,...
    'LineWidth',1,'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6]);
thetalim([-90 90]);
uprlim=max(histcounts(angle1,-90:10:90));
thetaticklabels({'S','SSE','ESE','E', 'ENE', 'NNE', 'N', 'NNW', ...
    'WNW', 'W', 'WSW',  'SSW',})
set(gca,'LineWidth',1,'FontSize',12,'RLim',[0 uprlim]);
elseif nargin==2 && style==1
polarhistogram(angle1.*pi/180, (-90:10:90).*pi/180,...
    'LineWidth',1,'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0]);
thetalim([-90 90]);
uprlim=max(histcounts(angle1,-90:10:90));
thetaticklabels({'S','SSE','ESE','E', 'ENE', 'NNE', 'N', 'NNW', ...
    'WNW', 'W', 'WSW',  'SSW',})
set(gca,'LineWidth',1,'FontSize',12,'RLim',[0 uprlim]);
elseif nargin==3
polarhistogram(angle1.*pi/180, (-90:10:90).*pi/180,...
    'LineWidth',1,'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6]);
hold on;
polarhistogram(angle2.*pi/180, (-90:10:90).*pi/180,...
    'LineWidth',1,'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0]);
thetalim([-90 90])
angle=[angle1;angle2];
uprlim=max(histcounts(angle,-90:10:90));
thetaticklabels({'S','SSE','ESE','E', 'ENE', 'NNE', 'N', 'NNW', ...
    'WNW', 'W', 'WSW',  'SSW',})
set(gca,'LineWidth',1,'FontSize',12,'RLim',[0 uprlim],'GridColor',[0 0 0]);
end
