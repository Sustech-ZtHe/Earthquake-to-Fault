function radius=radius_model(Mw)
%% this k-Mw relationship
L = [0.0005,0.002,0.005,0.02,0.06,0.15,0.5,2,6,25,100];
W = [0.0005,0.002,0.005,0.02,0.05,0.1,0.25,0.5,1,3.125,10];
k = [1.1,1.11,1.12,1.15,1.2,3,3.5,5,6.5,8,10];
Mag = [-3,-2,-1,0,1,2,3,4,5,6,7];
D = [0.001,0.0015,0.002,0.005,0.01,0.4,0.45,0.5,0.5,0.5,0.5].*10^-2;

%% D - fault displacement 
model = @(params, Mag)  params(1).*(1./(1+exp(params(2).*Mag+params(3))))+params(4);
initial_guess = [1, 1, 1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, Mag, D, [], [], options);
MDa = params_fit(1);
MDb = params_fit(2);
MDc = params_fit(3);
MDd = params_fit(4);
D= MDa.*(1./(1+exp(MDb.*Mw+MDc)))+MDd;

%% k - length to width ratio 
model = @(params, M) params(4).*log(1+exp(params(1).*Mag+params(2)))+params(3);
initial_guess = [1, 1, 1,1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, Mag, k, [], [], options);
Mka = params_fit(1);
Mkb = params_fit(2);
Mkc = params_fit(3);
Mkd = params_fit(4);

k=Mkd.*log(1+exp(Mka.*Mw+Mkb))+Mkc;
radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.^2.*3*10^10.*D),1/3)./k;


% %% Well 1994
% if 5<=Mw 
%     L=10.^((Mw-4.38)./1.49);
%     W=10.^((Mw-4.06)./2.25);
%     % k=L/W;
%     k=10.^((0.76*Mw-3.8)/3.3525);
% else
%     k=1;
% end
% radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.^2.*3*10^10.*D),1/3)./k;
end

