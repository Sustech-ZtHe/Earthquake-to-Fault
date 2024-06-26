function [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,el,linea,Vol,elliptic_points] = plotcov_3d(C1, m1,sp_km,weight,ratio,draw_style,s,loca)
% Visualise a 3x3 covariance matrix by drawing ellipsoids at 1, 2 and 3 STD.

% Input arguments:
%  C         3x3 Covariance matrix.
%  MU        Optional 1x3 array defining the centre of the ellipsoids.
%            The default value is [0,0,0].
%  VARARGIN  Any keyword arguments can be passed to `ellipsoid`.
%
% Output arguments:
%  H         3x1 vector of plot handles. One per ellipsoid. In order, they
%            are the handles to the ellipsoids at 1, 2 and 3 STD.
%
if nargin<4 || isempty(m1), m1=[0 0 0]; end
    % Find sorted eigenvectors and eigenvalues for C.
    [V1,D1]  = eig(C1);
    [~,ix] = sort(diag(D1), 'descend');
    D1      = D1(ix,ix);
    V1      = V1(:,ix);
    
    % Define the scales at which to draw the ellipsoids.
    %stds   = [1 2 3]; % 1, 2 and 3 standard deviations.
    %   stds   = [1 2]; % 1, 2 standard deviations.
    stds=[2 3];
%     conf  =normcdf(stds);
    conf   = 2 * normcdf(stds) - 1;
    scale  = chi2inv(conf, 2);  
    % Set up a unit sphere. Will be scaled to ellipsoids.
    [X,Y,Z] = sphere(16);
    e = [X(:), Y(:), Z(:)];

    MM=[mean(weight(:,1))/ratio(1,1),1,1;1,mean(weight(:,2))/ratio(2,1),1;1,1,mean(weight(:,3))/ratio(3,1)];
    % WEIGHT=[weight(:,1)/c(1,1),weight(:,2)/c(2,1),weight(:,3)/c(3,1)];
    % wc=cov(WEIGHT);
    wc=cov(weight);
    [V2,D2]  = eig(wc);
    VV = V1 * sqrt(D1 * scale(1).*[1 1 1;1 1 1;1 1 1]);
%         VV    = V1 * sqrt(D1 * scale(1).*sqrt(MM));
    VV=VV*V2;   %.*cov(weight);
    %     VV=VV*orth(V2);

    ee = bsxfun(@plus, e*(VV'), m1);
    elliptic_points=ee;

%% PRI_AXES_SIZE
%     [V,D]  = eig(cov(elliptic_points));
% %     elliptic_parameter=Ellipsoid3D_Fitting_LLS_SVD(ee');
%     pri_axes_size=sqrt(D);
%     es=sort([pri_axes_size(1,1),pri_axes_size(2,2),pri_axes_size(3,3)]);
%     c=es(1);b=es(2);a=es(3);
%     ec1=abs(sqrt(a^2-b^2)/a);
%     ec2=abs(sqrt(a^2-c^2)/a);
%     ecc=(ec1+ec2)/2;
% %     N=a*b*c/(0.1*111*0.1*111*cos(meanlatitude*pi/180)*0.1)*100
%     N=a*b*c/(s^2*0.1*0.1*0.1)*150;

[coeff, ~, latent, ~, ~]=pca(elliptic_points(:,1:3));
% [az,el,~]=cart2sph(mean(elliptic_points(:,1)),mean(elliptic_points(:,2)),mean(elliptic_points(:,3)));
% az=az*180/pi;
% el=el*180/pi;

hypotxy = hypot(mean(elliptic_points(:,1))-loca(1),mean(elliptic_points(:,2))-loca(2));
r = hypot(hypotxy,mean(elliptic_points(:,3))-loca(3));
% el = atan2(mean(elliptic_points(:,3))-loca(3),hypotxy)*180/pi;
% az = atan2(mean(elliptic_points(:,2))-loca(2),mean(elliptic_points(:,1))-loca(1))*180/pi;

p1_u=2*sqrt(latent(1))*coeff(:, 1)'+m1;
p2_u=2*sqrt(latent(2))*coeff(:, 2)'+m1;
p3_u=2*sqrt(latent(3))*coeff(:, 3)'+m1;
p1_d=m1-2*sqrt(latent(1))*coeff(:, 1)';
p2_d=m1-2*sqrt(latent(2))*coeff(:, 2)';
p3_d=m1-2*sqrt(latent(3))*coeff(:, 3)';

% p1_u=[0.7079 1.968 3.734]
% p1_d=[0.812 1.956 3.68]
% p2_u=[0.296 1.019 3.148]
% p2_d=[0.336 0.992 3.171]
el = atan2(abs(p2_u(3)-p2_d(3)),pdist2([p2_u(1),p2_u(2)],[p2_d(1),p2_d(2)]))*180/pi;
az = atan2(abs(p1_u(2)-p1_d(2)),p1_u(1)-p1_d(1))*180/pi;
% atan2(1,1)*180/pi

a=sqrt(latent(1));
b=sqrt(latent(2));
c=sqrt(latent(3));
S=4*pi*a^2+2*pi*a*b;

% sphe=4*pi*(a*b*c)^(2/3)/S;
% linea=1-sphe;

linea=(b-c)/b;
Vol=a*b*c/(s^2*0.1*0.1*0.1)*100;
 
% % figure;scatter3(EP_km(:,1),EP_km(:,2),EP_km(:,3));hold on
% % figure;
% % scatter3(sp_km(:,1),sp_km(:,2),sp_km(:,3));hold on
% scatter3(elliptic_points(:,1),elliptic_points(:,2),elliptic_points(:,3))
% hold on
% scatter3(m1(1,1),m1(1,2),m1(1,3),100,'k','filled','o');
% scatter3(p1_u(1,1),p1_u(1,2),p1_u(1,3),100,'green','filled','o');
% scatter3(p2_u(1,1),p2_u(1,2),p2_u(1,3),100,'red','filled','o');
% scatter3(p3_u(1,1),p3_u(1,2),p3_u(1,3),100,'blue','filled','o');
% 
% scatter3(p1_d(1,1),p1_d(1,2),p1_d(1,3),100,'green','filled','o');
% scatter3(p2_d(1,1),p2_d(1,2),p2_d(1,3),100,'red','filled','o');
% scatter3(p3_d(1,1),p3_d(1,2),p3_d(1,3),100,'blue','filled','o');


% axis equal;view(0,90)
% scatter3(p1_d(1,1),p1_d(1,2),p1_d(1,3),100,'magenta','filled','o');
% scatter3(p3(1,1),p3(1,2),p3(1,3),100,'blue','filled','o');
% pdist2([mean(elliptic_points(:,1)),mean(elliptic_points(:,2)),mean(elliptic_points(:,3))],[0 0 0])
% cart2pol(coeff(1,:),coeff(2,:),coeff(3,:))*180/pi
% [az,el,r]=cart2sph(1,1,sqrt(2))
% el*180/pi
% [V,D]=eig(cov(elliptic_points(:,1:3)))
% [~,ix] = sort(diag(D), 'descend');
%     D      = D(ix,ix);
%     V      = V(:,ix);

% loca=median(hypo_out_points(:,1:3));
% hypotxy = hypot(mean(elliptic_points(:,1))-loca(1),mean(elliptic_points(:,2))-loca(2));
% r = hypot(hypotxy,mean(elliptic_points(:,3))-loca(3));
% el = atan2(mean(elliptic_points(:,3))-loca(3),hypotxy)*180/pi
% az = atan2(mean(elliptic_points(:,2))-loca(2),mean(elliptic_points(:,1))-loca(1))*180/pi
% figure;scatter3(elliptic_points(:,1),elliptic_points(:,2),elliptic_points(:,3));hold on
% scatter3(loca(1),loca(2),loca(3),100,'red','filled');
% 
% EP=KMtoDEG(elliptic_points,la0,lo0);
% loca=KMtoDEG(loca,la0,lo0);
% hypotxy = hypot(mean(EP(:,1))-loca(1),mean(EP(:,2))-loca(2))
% elev = atan2(mean(EP(:,3))-loca(3),hypotxy)*180/pi
% az = atan2(mean(EP(:,2))-loca(2),mean(EP(:,1))-loca(1))*180/pi
% r = hypot(hypotxy,mean(EP(:,3))-loca(3));
% % [az,el,~]=cart2sph(mean(EP(:,1)),mean(EP(:,2)),mean(EP(:,3)));
% % az=az*180/pi;
% % el=el*180/pi;
% figure;scatter3(EP(:,1),EP(:,2),EP(:,3));hold on
% scatter3(loca(1),loca(2),loca(3),100,'red','filled');
% xlim([-98 -97.6]);
% ylim([36.8 36.86]);
% zlim([3 10]);
% set(gca,'ZDir','reverse')
% atan2(1,-1)*180/pi
% axis equal
% atan2(0.02,0.1)

%%%%%%%%%Draw%%%%%%%%%%
if draw_style==1
    % Line styles (and handles) for the ellipsoids.
    styles = {
    {'LineStyle','--','LineWidth',2};    % 1 STD
    {'LineStyle',':','LineWidth',2,'color',[unique(sp_km(:,11)) unique(sp_km(:,12)) unique(sp_km(:,13))]}; % 2 STD
    };
    h = zeros(numel(stds), 1);
    washold = ishold;
    MM=[mean(weight(:,1))/ratio(1,1),1,1;1,mean(weight(:,2))/ratio(2,1),1;1,1,mean(weight(:,3))/ratio(3,1)];
    % WEIGHT=[weight(:,1)/c(1,1),weight(:,2)/c(2,1),weight(:,3)/c(3,1)];
    % wc=cov(WEIGHT);
    wc=cov(weight);
    [V2,D2]  = eig(wc);
    %   [~,ix] = sort(diag(D2), 'descend');
    %   D2      = D2(ix,ix);
    %   V2      = V2(:,ix);
    % Set up the i-th scaled ellipsoid. 
    for i = 1:numel(stds) 
        % Set up the i-th scaled ellipsoid.
        %     VV    = V1 * sqrt(D1 * scale(i));
        VV    = V1 * sqrt(D1 * scale(i).*sqrt(MM));
        %     VV=VV*V2;   %.*cov(weight);
        %     VV=VV*power(V2,3);
        VV=VV*V2;
        ee    = bsxfun(@plus, e*(VV'), m1);
        ex    = reshape(ee(:,1), size(X));
        ey    = reshape(ee(:,2), size(Y));
        ez    = reshape(ee(:,3), size(Z));
        elliptic_points=ee;
%         elliptic_parameter=Ellipsoid3D_Fitting_LLS_SVD(ee');

        % Draw with different properties for different ellipsoids.
        % Also, pick the colour from the first ellipsoid for aesthetics.    
%         if i==1 & isreal(ee)
%             h(i) = surf(ex,ey,ez,styles{i}{:});
%             hold on;
%             shading interp
%             %       'FaceColor','r',
%             colormap("gray")
%             set(h(i),'EdgeColor', [unique(subS(:,11)) unique(subS(:,12)) unique(subS(:,13))]);
%         else
%             h(i) = plot3(ee(:,1),ee(:,2),ee(:,3),styles{i}{:});
%         end  
    end
      
    plot3(ee(:,1),ee(:,2),ee(:,3),styles{i}{:});

    if ~washold
        hold off;
    end
end
%%%%%%%%%END%%%%%%%%%%

end
% p1_u=[0.87 1.92 3.34]
% p1_d=[0.69 1.59 3.31]
% p2_u=[0.779 1.76 3.29]
% p2_d=[0.786 1.75 3.36]

% p1_u=[0.824 1.597 3.33]
% p1_d=[0.704 1.37 3.35]
% p2_u=[0.769 1.481 3.29]
% p2_d=[0.759 1.49 3.38]

% p1_u=[1.0633 1.848 3.443]
% p1_d=[0.725 1.269 3.327]
% p2_u=[0.904 1.564 3.33]
% p2_d=[0.879 1.544 3.432]

% p1_u=[0.36 1.0 3.11]
% p1_d=[0.27 0.988 3.19]
% p2_u=[0.296 1.019 3.148]
% p2_d=[0.336 0.992 3.171]

% p1_u=[0.818 2.32 3.73]
% p1_d=[0.725 2.1 3.72]
% 
% p1_u=[0.714 2.02 3.74]
% p1_d=[0.842 2.18 3.7]

% p1_u=[0.7079 1.968 3.734]
% p1_d=[0.812 1.956 3.68]