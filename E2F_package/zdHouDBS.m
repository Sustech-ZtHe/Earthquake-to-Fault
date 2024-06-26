function [FE_km,FE]=zdHouDBS(C,l,f_sort,line_1,line_nums_2,la0,lo0,minPoi,Mc,hypotout,colortemplate,LowerMag,MidMag,UpMag,MODEL)
n=1;M1=0;M2=0;FE_km=[];
% figure;
for i=1:line_1(end)
    im=1;jj=1;
    if ismember(i,unique(f_sort(:,14)))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort(line_nums_2(k,1):line_nums_2(k,2)-1,:);
        nw=1;

%         scatter3(ftest4(:,1),ftest4(:,2),ftest4(:,3),10, ftest4(:,11:13), 'filled');hold on 
%         text(mean(ftest4(:,1)), mean(ftest4(:,2)), mean(ftest4(:,3)),num2str(i), ...
%         'FontSize',20,'FontName', 'Arial','Color',[ftest4(1,11) ftest4(1,12) ftest4(1,13)], ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
%         'FontWeight', 'bold','FontAngle', 'italic');
%         xlim([0 2]);ylim([0 3]);view(0,90)

        while 1
            MaxMag=max(ftest4(:,4));
%             CauMag=(log10(length(ftest4(:,1)))-2.46)/1.13;
%             if size(ftest4,1)~=1
%                 cp=median(ftest4(:,1:3));
%             else
%                 cp=ftest4(:,1:3);
%             end
%             if cp(3)<3.5
%                 if (-117.25<cp(1)&&cp(1)<-117.245 && 54.33<cp(2)&&cp(2)<54.35) ...
%                     || (-117.245<cp(1)&&cp(1)<-117.24 && 54.335<cp(2)&&cp(2)<54.34) ...
%                     || (-117.23<cp(1)&&cp(1)<-117.225 && 54.345<cp(2)&&cp(2)<54.355)
%                     b=0.92;a=2.23;
%                 else 
%                     b=1.69;a=1.72;
%                 end
%             else
%                 b=2.47;a=-0.89;
%             end
%             CauMag=(log10(length(ftest4(:,1)))-a)/b;
%             CauMag=(log10(length(ftest4(ftest4(:,4)>-1.6,1)))-1.13*1.6)/1.13;
%             CauMag=(log10(length(ftest4(ftest4(:,4)>-1.6,1)))+1.5)/1.13-1.6;
            if im==1
                [a,b,Y,CauMag]=Mmain(ftest4,Mc,nw,i,hypotout);
                im=2;
            elseif im==2
%                 CauMag=(log10(size(ftest4,1))-a)/b+mean(Y(:,1));
                CauMag=(log10(size(ftest4,1))-a)/b+Mc;
            end
%                 CM=[CM;[MaxMag,CauMag,size(ftest4,1),nw]];
            Mag=max(MaxMag,CauMag);
            if  MaxMag>CauMag
                M1=M1+1;
            else
                M2=M2+1;
            end
            c=C(1,l);
            radius=c*radius_model(Mag,LowerMag,MidMag,UpMag,MODEL);
%             radius=c*0.001*power(10^(1.5*(Mag+6))/(3*10^5*2*pi),1/3);
%                 0.001*power(10^(1.5*1+16.1)*10^(-7)/(3*10^7),1/3);
%                 0.001*power(10^(1.5*(3+6.033))/(3*10^7),1/3);
%             radius=c*0.001*(1/log(1.58+MaxMag))*power(10^(1.5*MaxMag+16.1)*10^(-7)/(3*10^5*pi),1/3);
%             radius=c*power(10^(1.5*Mag+16.05)*10^(-10)/(3*10^11*10^(-5)*10^(10)*0.001*0.1),1/3);
%             radius=c*0.001*power(10^(1.5*Mag+16.05)*10^(-7)/(3*10^5*pi),1/3);
            
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);

            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            nw=1;
            if all(idx_5 == -1) || isempty(idx_5)
                break
            end
            color_5=rand(length(unique(idx_5)),3);
%             for j=1:length(unique(idx_5))
%                 F3=ftest4(find(idx_5==j),:);
%                 F3(:,11)=color_5(j,1);F3(:,12)=color_5(j,2);F3(:,13)=color_5(j,3);F3(:,14)=n;
%                 FE_km=[FE_km;F3];
%                 n=n+1;
%             end

            F3=ftest4(find(idx_5==1),:);
%             F3(:,11)=color_5(1,1);F3(:,12)=color_5(1,2);F3(:,13)=color_5(1,3);
            F3(:,11)=colortemplate(i,3*(jj-1)+1);F3(:,12)=colortemplate(i,3*(jj-1)+2);F3(:,13)=colortemplate(i,3*(jj-1)+3);
            F3(:,14)=n;
            FE_km=[FE_km;F3];
            n=n+1;
            jj=jj+1;

            ftest4(idx_5==1,:)=[];
            if isempty(ftest4)
                break
            end
        end
    end
end
FE=KMtoDEG(FE_km,la0,lo0);
end


