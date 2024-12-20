function POLAR_gram(angle1,style,angle2)
if nargin==2 && style==0 && ~isempty(angle1)
    polarhistogram([angle1,angle1-180].*pi/180, deg2rad(0:10:360),...
        'LineWidth',1,'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6]);
    pax = gca;
    pax.ThetaZeroLocation = 'top'; % 将 0° 移动到顶部（北）
    pax.ThetaDir = 'clockwise';    % 顺时针方向
    set(pax,'LineWidth',1.5,'FontSize',12);
elseif nargin==2 && style==1 && ~isempty(angle1)
    polarhistogram([angle1,angle1-180].*pi/180, deg2rad(0:10:360),...
        'LineWidth',1,'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0]);
    pax = gca;
    pax.ThetaZeroLocation = 'top'; % 将 0° 移动到顶部（北）
    pax.ThetaDir = 'clockwise';    % 顺时针方向
    set(pax,'LineWidth',1.5,'FontSize',12);
elseif nargin==3 && (~isempty(angle1) || ~isempty(angle2))
    polarhistogram([angle1,angle1-180].*pi/180, deg2rad(0:10:360),...
        'LineWidth',1,'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6]);
    hold on;
    polarhistogram([angle2,angle2-180].*pi/180, deg2rad(0:10:360),...
        'LineWidth',1,'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0]);
    angle=[angle1;angle2];
    pax = gca;
    pax.ThetaZeroLocation = 'top'; % 将 0° 移动到顶部（北）
    pax.ThetaDir = 'clockwise';    % 顺时针方向
    set(pax,'LineWidth',1.5,'FontSize',12);
end
