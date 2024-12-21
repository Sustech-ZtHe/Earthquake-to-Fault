function radius=radius_model(Mw)
% Mw=-3:0.01:7
L = [0.0005,0.002,0.005,0.02,0.06,0.15,0.5,2,6,25,100];
W = [0.0005,0.002,0.005,0.02,0.05,0.1,0.25,0.5,1,3.125,10];
% k = [1,1.01,1.02,1.05,1.5,3,3.5,4,5,5.5,6];
k = [1.1,1.11,1.12,1.15,1.2,3,3.5,5,6.5,8,10];
Mag = [-3,-2,-1,0,1,2,3,4,5,6,7];
D = [0.001,0.0015,0.002,0.005,0.01,0.8,0.85,0.9,1,1,1].*10^-2;
% D = [0.01,0.02,0.03,0.05,0.01,0.4,0.5,0.6,0.7,0.8,1].*10^-4;

% L = [0.02,0.06,0.15,0.5,2,6,25,100];
% W = [0.02,0.05,0.1,0.25,0.5,1,3.125,10];
% k = [1,1.1,1.5,2,4,6,8,10];
% Mag = [0,1,2,3,4,5,6,7];
% D = [0.002,0.006,0.015,0.05,0.1,0.2,0.4,1].*10^-2;

%% L - fault length
model = @(params, L)  params(1).*log10(L)+params(2);
initial_guess = [1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, L, Mag, [], [], options);
MLa = params_fit(1);
MLb = params_fit(2);
% M=MLa .*log10(L)+MLb;

%% W - fault width
model = @(params, W) params(1).*log10(W)+params(2);
initial_guess = [1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, W, Mag, [], [], options);
MWa = params_fit(1);
MWb = params_fit(2);
% M=MWa .*log10(W)+MWb;

%% D - fault displacement 
model = @(params, Mag)  params(1).*(1./(1+exp(params(2).*Mag+params(3))))+params(4);
% model = @(params, Mag)  params(1).*exp(params(2).*Mag+params(3))+params(4);
% model = @(params, Mag)  log(params(1)+exp(params(2).*Mag+params(3)));
initial_guess = [1, 1, 1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, Mag, D, [], [], options);
MDa = params_fit(1);
MDb = params_fit(2);
MDc = params_fit(3);
MDd = params_fit(4);
% D= MDa.*(1./(1+exp(MDb.*Mw+MDc)))+MDd;

% LM=(10.^((M-MLb)./MLa));
% WM=(10.^((M-MWb)./MWa));
%% k - length to width ratio 
model = @(params, M) params(4).*log(1+exp(params(1).*Mag+params(2)))+params(3);
% model = @(params, M) params(1).*log10((10.^((M-MLb)./MLa)) ./ (10.^((M-MWb)./MWa)))+params(2);
initial_guess = [1, 1, 1,1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, Mag, k, [], [], options);
Mka = params_fit(1);
Mkb = params_fit(2);
Mkc = params_fit(3);
Mkd = params_fit(4);
% k=Mka.*log10((10.^((Mw-MLb)./MLa)) ./ (10.^((Mw-MWb)./MWa)))+Mkb;

k=Mkd.*log(1+exp(Mka.*Mw+Mkb))+Mkc;
D= MDa.*(1./(1+exp(MDb.*Mw+MDc)))+MDd;
% D= MDa.*exp(MDb.*Mw+MDc)+MDd;
% D=log(MDa+exp(MDb.*Mw+MDc));
radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.^2.*3*10^10.*D),1/3)./k;


% figure
% LowerMag=2;
% UpMag=7;
% R=[];
% for i=1:length(Mw)
% if LowerMag<=Mw(i) & Mw(i)<UpMag
%     k=5;
%     radius=0.001*power(10^(1.5*Mw(i)+9.1)/(k^2*3*10^10/1000),1/3)/k;
% elseif Mw(i)<LowerMag
%     k=1;
%     radius=0.001*power(10^(1.5*Mw(i)+9.1)/(k^2*3*10^10/1000),1/3)/k;
% elseif Mw(i)>=UpMag
%     k=10;
%     radius=0.001*power(10^(1.5*Mw(i)+9.1)/(k^2*3*10^10/1000),1/3)/k;
% end
% R=[R;radius];
% end
% plot(Mw,R)
% 
% hold on
% % k=Mka.*log10((10.^((Mw-MLb)./MLa)) ./ (10.^((Mw-MWb)./MWa)))+Mkb;
% k=Mkd.*log(1+exp(Mka.*Mw+Mkb))+Mkc;
% % k=p(1)*Mw.^2+p(2)*Mw+p(3);
% % k=p(1)*Mw.^7+p(2)*Mw.^6+p(3)*Mw.^5+p(4)*Mw.^4+p(5)*Mw.^3+p(6)*Mw.^2+p(7)*Mw+p(8);
% % figure;plot(Mw,k)
% D=MDa.*(1./(1+exp(MDb.*Mw+MDc)))+MDd;
% % D= MDa.*exp(MDb.*Mw+MDc)+MDd;
% % D=log(MDa+exp(MDb.*Mw+MDc));
% radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.^2.*3*10^10.*D),1/3)./k;
% % radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.*3*10^10.*D),1/2)./2;
% plot(Mw,radius);

end

