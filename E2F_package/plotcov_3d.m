function [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,strike,dip,linea,ellipsoid_points] = plotcov_3d(m1,sp_km)
%% Long Mid Short Up-Long Bo-Long Up-Mid Bo-Mid Up-Short Bo-Short strike dip ... ...
conf = 0.9;
scale_el = sqrt(chi2inv(conf, 3)); 
[coeff, ~, latent, ~, ~]=pca(sp_km(:,1:3));
p1_u=scale_el*sqrt(latent(1))*coeff(:, 1)'+m1;
p2_u=scale_el*sqrt(latent(2))*coeff(:, 2)'+m1;
p3_u=scale_el*sqrt(latent(3))*coeff(:, 3)'+m1;
p1_d=m1-scale_el*sqrt(latent(1))*coeff(:, 1)';
p2_d=m1-scale_el*sqrt(latent(2))*coeff(:, 2)';
p3_d=m1-scale_el*sqrt(latent(3))*coeff(:, 3)';
V_a = p1_u - p1_d; 
V_b = p2_u - p2_d; 
V_c = p3_u - p3_d; 
V_ah = [V_a(1), V_a(2), 0]; 
strike = atan2(V_ah(1), V_ah(2)); 
strike = strike * 180 / pi; 
if strike>180
    strike=strike-180;
elseif strike<0
    strike=strike+180;
end
n = cross(V_a, V_c); 
dip = atan2(norm(cross(n, [1, 0, 0])), dot(n, [1, 0, 0]));
dip = dip * 180 / pi; 
if dip>90
    dip=180-dip;
end
V_ap = V_a - dot(V_a, n) / dot(n, n) * n; 
V_bp = [V_b(1), V_b(2), 0]; 
a=scale_el*sqrt(latent(1));
b=scale_el*sqrt(latent(2));
c=scale_el*sqrt(latent(3));
linea=1-b/a*c/a;

theta = linspace(0, 2*pi, 20);
phi = linspace(-pi/2, pi/2, 10);
[Theta, Phi] = meshgrid(theta, phi);
X = m1(1) + a * coeff(1,1) * cos(Phi) .* cos(Theta) ...
        + b * coeff(1,2) * cos(Phi) .* sin(Theta) ...
        + c * coeff(1,3) * sin(Phi);
Y = m1(2) + a * coeff(2,1) * cos(Phi) .* cos(Theta) ...
        + b * coeff(2,2) * cos(Phi) .* sin(Theta) ...
        + c * coeff(2,3) * sin(Phi);
Z = m1(3) + a * coeff(3,1) * cos(Phi) .* cos(Theta) ...
        + b * coeff(3,2) * cos(Phi) .* sin(Theta) ...
        + c * coeff(3,3) * sin(Phi);
X = X(:);
Y = Y(:);
Z = Z(:);
ellipsoid_points = [X, Y, Z];
end
