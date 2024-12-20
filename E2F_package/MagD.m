function D=MagD(Mw)
Mag = [-3,-2,-1,0,1,2,3,4,5,6,7];
D = [0.001,0.0015,0.002,0.005,0.01,0.8,0.85,0.9,1,1,1].*10^-2;

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

end

