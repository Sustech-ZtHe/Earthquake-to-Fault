 function radius=radius_model(Mw,LowerMag,MidMag,UpMag,MODEL)
if MODEL==1
%     if LowerMag<=Mw & Mw<UpMag
% %         radius=0.001*power(10^(1.5*Mw+9)/(3*10^6*pi*2),1/3)/10;      
%         radius=sqrt(power(10^(1.5*Mw+16.1)/(1.23*10^22),2/3)/pi)/10;
%     elseif Mw<LowerMag
% %         radius=0.001*power(10^(1.5*Mw+9)/(3*10^7*pi*2),1/3);       
%         radius=sqrt(power(10^(1.5*Mw+16.1)/(1.23*10^22),2/3)/pi);
%     elseif Mw>=UpMag
% %         radius=0.001*power(10^(1.5*Mw+9)/(3*10^5*pi*2),1/3)/100;
%         radius=sqrt(power(10^(1.5*Mw+16.1)/(1.23*10^22),2/3)/pi)/100;
%     end


    if LowerMag<=Mw & Mw<UpMag
%         radius=0.001*power(10^(1.5*Mw+9)/(3*10^6*pi*2),1/3)/10;      

%         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/5)/5;
        k=5;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
%         radius=10^-5*power(10^(1.5*Mw+16.1)/(k^2*3*10^11/1000),1/3)/k;
    elseif Mw<LowerMag
%         radius=0.001*power(10^(1.5*Mw+9)/(3*10^7*pi*2),1/3);   

%         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/1)/1;
        k=1;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
%         radius=10^-5*power(10^(1.5*Mw+16.1)/(k^2*3*10^11/1000),1/3)/k;
    elseif Mw>=UpMag
%         radius=0.001*power(10^(1.5*Mw+9)/(3*10^5*pi*2),1/3)/100;

%         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/10)/10;
        k=10;
        radius=0.001*power(10^(1.5*Mw+9.1)/(k^2*3*10^10/1000),1/3)/k;
%         radius=10^-5*power(10^(1.5*Mw+16.1)/(k^2*3*10^11/1000),1/3)/k;
    end

% Mw=3;k=1    
% 10^-3*power( 10^(1.5*Mw+9.1)/(k^2*3*10^10), 1/3)
% 10^-5*power( 10^(1.5*Mw+16.1)/(k^2*3*10^11), 1/3)
% 10^(1.5*5.3+9.1)
% k=5
% Mw=6.2
% power(10^(1.5*(Mw+10.7))/(k/10*3*10^16),1/3)
% 0.001*power(10^(1.5*(Mw+6.07))/(k*k/100*3*10^10),1/3)
% sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/k)/k

%     sqrt(power(10^(1.5*(7+10.7))/(1.23*10^22),2/3)/pi/10)/10
%     sqrt(10^(7+10.7)/(1.148*10^(44/3)*pi*10))/10
% sqrt(power(10^(1.5*3+16.1)/(1.23*10^22),2/3)/pi/10)/sqrt(10)

%     uNM=6*10^10;k=1*10^-4;
%     if Mw<LowerMag %<2
%         %     1:1
%         c=1;
% %         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi)/2;
%         radius=0.001*power(10^(1.5*(Mw+6.07))/(uNM*(c^2)*k*pi),1/3)/4;
%     elseif LowerMag<=Mw & Mw<MidMag  %2< <5
%         %     1:5 
%         c=2;
% %         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/5)/2;
%         radius=0.001*power(10^(1.5*(Mw+6.07))/(uNM*(c^2)*k*pi),1/3)/4;
%     elseif MidMag<=Mw & Mw<UpMag  %5< <7
%         %     1：10 
%         c=5;
% %         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/10)/2;
%         radius=0.001*power(10^(1.5*(Mw+6.07))/(uNM*(c^2)*k*pi),1/3)/4;
%     elseif Mw>=UpMag %>7
%         %     1：20  
%         c=10;
% %         radius=sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi/20)/2;
%         radius=0.001*power(10^(1.5*(Mw+6.07))/(uNM*(c^2)*k*pi),1/3)/5;
%     end
% c=10;uNM=6*10^10;k=1*10^-4;
% %     0.001*power(10^(1.5*(3+6.07))/(uNM*(2*c^2)*k*pi),1/3)
%     0.001*power(10^(1.5*(3+6.07))/(uNM*(c^3)*k*pi/4),1/3)/10
%     sqrt(power(10^(1.5*3+16.1)/(1.23*10^22),2/3)/pi)/10
%     0.001*power(10^(1.5*(3+6.07))/(uNM*k*(1^2)*pi/4),1/3)
% sqrt(power(10^(1.5*6.5+16.1)/(1.23*10^22),2/3)/pi)^0.5
%        uNM=6*10^10; c=1;k=10^-4; Mw=2;
%     0.001*power(10^(1.5*(Mw+6.07))/(uNM*c*k*pi*c),1/3)
%     sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/pi)
%     c=5;Mw=3
%     0.001*power(10^(1.5*(Mw+6.07))/(uNM*c*k*pi*c),1/3)
%     sqrt(power(10^(1.5*(Mw+10.7))/(1.23*10^22),2/3)/5/pi)
%     0.001*power(10^(1.5*(5.1+6.07))/(uNM*10*k*0.1*pi),1/3)/2
%     sqrt(power(10^(1.5*(3+10.7))/(1.23*10^22),2/3)/pi)/10
% 10^((7-6.93)/0.82)

%     10^((1.5*5+9.1-log10(u*C2))/3)
%     sqrt(power(10^(1.5*3+16.1)/(1.23*10^22),2/3)/pi/2)
% 
%     L=10^((1.5*Mw+9.1-log10(u*C2))/3);
%     if L>5.5
%         radius=1.7*L^(2/3)/2;
%     else 
%         radius=L/2;
%     end

elseif MODEL==2
    radius=power(power(10^(1.5*Mw+16.1)/(1.23*10^22),2/3)/(1.7^1.5),2/5);
end
end
