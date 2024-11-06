 function radius=radius_model(Mw,LowerMag,UpMag,MODEL)
if MODEL==1
    if LowerMag<=Mw & Mw<UpMag
        k=5;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
    elseif Mw<LowerMag
        k=1;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
    elseif Mw>=UpMag
        k=10;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
    end
elseif MODEL==2
    radius=power(power(10^(1.5*Mw+16.1)/(1.23*10^22),2/3)/(1.7^1.5),2/5);
end
end
