function radius=radius_model(Mw)
L = [0.0005,0.002,0.005,0.02,0.06,0.15,0.5,2,6,25,100];
W = [0.0005,0.002,0.005,0.02,0.05,0.1,0.25,0.5,1,3.125,10];
k = [1.1,1.11,1.12,1.15,1.2,3,3.5,5,6.5,8,10];
Mag = [-3,-2,-1,0,1,2,3,4,5,6,7];
D = [0.001,0.0015,0.002,0.005,0.01,0.8,0.85,0.9,1,1,1].*10^-2;

%% L - fault length
model = @(params, L)  params(1).*log10(L)+params(2);
initial_guess = [1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, L, Mag, [], [], options);
MLa = params_fit(1);
MLb = params_fit(2);

%% W - fault width
model = @(params, W) params(1).*log10(W)+params(2);
initial_guess = [1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, W, Mag, [], [], options);
MWa = params_fit(1);
MWb = params_fit(2);

%% D - fault displacement 
model = @(params, Mag)  params(1).*(1./(1+exp(params(2).*Mag+params(3))))+params(4);
initial_guess = [1, 1, 1, 1];
options = optimset('Display', 'off');
params_fit = lsqcurvefit(model, initial_guess, Mag, D, [], [], options);
MDa = params_fit(1);
MDb = params_fit(2);
MDc = params_fit(3);
MDd = params_fit(4);

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
D= MDa.*(1./(1+exp(MDb.*Mw+MDc)))+MDd;
radius=0.001.*power(10.^(1.5.*Mw+9.1)./(k.^2.*3*10^10.*D),1/3)./k;

end

