function SUBfault_con(L,C,f_sort,line_1,line_nums_2,hypotout,FaultpaRameters,Mc,LowerMag,MidMag,UpMag,MODEL,minPoi)
figure('Position', [100, 100, 880, 1560]); 
colortemplate=rand(200,300);Kused=[];
M1=0;M2=0;
n=1;FE_km=[];CM=[];
nw=1;RabY=[];
ii=1;
for i=1:line_1(end)
    im=1;jj=1;
    if ismember(i,unique(f_sort(:,14)))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort(line_nums_2(k,1):line_nums_2(k,2)-1,:);
        while 1
            MaxMag=max(ftest4(:,4));
            if im==1
                [a,b,Y,CauMag]=Mmain(ftest4,Mc,nw,i,hypotout);
                im=2;
            elseif im==2
                CauMag=(log10(size(ftest4,1))-a)/b+Mc;
            end
            Ns=(log10(size(ftest4,1))-a)/b+Mc;
            Mag=max(MaxMag,CauMag);

            if  MaxMag>CauMag
                M1=M1+1;
            else
                M2=M2+1;
            end
            c=C(1,L);
            radius=c*radius_model(Mag,LowerMag,MidMag,UpMag,MODEL);
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);
            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            nw=nw+1;
            if all(idx_5 == -1) || isempty(idx_5)
                RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),Ns,Mag];
                break
            end
            color_5=rand(length(unique(idx_5)),3);              
            F3=ftest4(find(idx_5==1),:);
%                 F3(:,11)=color_5(1,1);F3(:,12)=color_5(1,2);F3(:,13)=color_5(1,3);
            F3(:,11)=colortemplate(i,3*(jj-1)+1);F3(:,12)=colortemplate(i,3*(jj-1)+2);F3(:,13)=colortemplate(i,3*(jj-1)+3);
            F3(:,14)=n;
            FE_km=[FE_km;F3];
            n=n+1;
            ii=ii+1;
            jj=jj+1;

            ftest4(idx_5==1,:)=[];
            if isempty(ftest4)
                break
            end
        end
    end
end 
l=find(FaultpaRameters(:,29)==L);
for i=1:size(l)
    x=[FaultpaRameters(l(i),7)',FaultpaRameters(l(i),7)',FaultpaRameters(l(i),10)',FaultpaRameters(l(i),10)'];
    y=[FaultpaRameters(l(i),8)',FaultpaRameters(l(i),8)',FaultpaRameters(l(i),11)',FaultpaRameters(l(i),11)'];
    z=[FaultpaRameters(l(i),15)',FaultpaRameters(l(i),18)',FaultpaRameters(l(i),18)',FaultpaRameters(l(i),15)'];
    fill3(x, y, z, [255 204 204]./255);hold on
end
scatter3(FE_km(:,1),FE_km(:,2),FE_km(:,3),2.*power(exp(FE_km(:,4)),1/5),FE_km(:,11:13),'filled')
set(gca,'ZDir','reverse')
view(0,90)
axis off