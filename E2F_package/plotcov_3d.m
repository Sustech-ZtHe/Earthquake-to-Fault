function [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,el,linea,Vol,elliptic_points] = plotcov_3d(C1, m1,sp_km,weight,ratio,draw_style,s,loca)
if nargin<4 || isempty(m1), m1=[0 0 0]; end
    % Find sorted eigenvectors and eigenvalues for C.
    [V1,D1]  = eig(C1);
    [~,ix] = sort(diag(D1), 'descend');
    D1      = D1(ix,ix);
    V1      = V1(:,ix);

    stds=[2 3];
    conf   = 2 * normcdf(stds) - 1;
    scale  = chi2inv(conf, 2);  
    % Set up a unit sphere. Will be scaled to ellipsoids.
    [X,Y,Z] = sphere(16);
    e = [X(:), Y(:), Z(:)];

    MM=[mean(weight(:,1))/ratio(1,1),1,1;1,mean(weight(:,2))/ratio(2,1),1;1,1,mean(weight(:,3))/ratio(3,1)];
    wc=cov(weight);
    [V2,D2]  = eig(wc);
    VV = V1 * sqrt(D1 * scale(1).*[1 1 1;1 1 1;1 1 1]);
    VV=VV*V2;   %.*cov(weight);

    ee = bsxfun(@plus, e*(VV'), m1);
    elliptic_points=ee;
    [coeff, ~, latent, ~, ~]=pca(elliptic_points(:,1:3));
    hypotxy = hypot(mean(elliptic_points(:,1))-loca(1),mean(elliptic_points(:,2))-loca(2));
    r = hypot(hypotxy,mean(elliptic_points(:,3))-loca(3));
    p1_u=2*sqrt(latent(1))*coeff(:, 1)'+m1;
    p2_u=2*sqrt(latent(2))*coeff(:, 2)'+m1;
    p3_u=2*sqrt(latent(3))*coeff(:, 3)'+m1;
    p1_d=m1-2*sqrt(latent(1))*coeff(:, 1)';
    p2_d=m1-2*sqrt(latent(2))*coeff(:, 2)';
    p3_d=m1-2*sqrt(latent(3))*coeff(:, 3)';
    
    el = atan2(abs(p2_u(3)-p2_d(3)),pdist2([p2_u(1),p2_u(2)],[p2_d(1),p2_d(2)]))*180/pi;
    az = atan2(abs(p1_u(2)-p1_d(2)),p1_u(1)-p1_d(1))*180/pi;
    
    a=sqrt(latent(1));
    b=sqrt(latent(2));
    c=sqrt(latent(3));
    S=4*pi*a^2+2*pi*a*b;
    
    linea=(b-c)/b;
    Vol=a*b*c/(s^2*0.1*0.1*0.1)*100;

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
    end
      
    plot3(ee(:,1),ee(:,2),ee(:,3),styles{i}{:});

    if ~washold
        hold off;
    end
end
%%%%%%%%%END%%%%%%%%%%

end
