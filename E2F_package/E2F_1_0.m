clear all;clc
addpath('/home/me/Documents/E2F/H3D/hough-3d-lines-master/E2F_package/')
%% EARTHQUAKE CATALOGUE
% 'ToC2ME'
% 'qinghai' 
% 'foxcreek' 
% 'peaceriver' 
% 'ratonbasin' 
% 'reddeer' 
% 'weiyuan' 
% 'turkey_afad' 
% 'turkey_usgs'
% 'turkey_ml'
% 'qinghai'
% 'ridgecrest'
% 'wenchuan'
% 'SouthernC_12dev'
% 'SouthernC_9.5dev'
% 'Southernc2010'
% 'Pawnee'
% 'Grant'
% 'Okamc'


file_path='/home/me/Documents/E2F/H3D/hough-3d-lines-master/E2F_package/SaveData/wy_catalog_s1.dat'

MODEL=1;LowerMag=2;UpMag=7;
RC=[];[R,C,s,minPoi,Mc,na,hypo]=Detect_Region(file_path);
lon=hypo(:,9);
lat=hypo(:,8);
depth=hypo(:,10);
mag=hypo(:,11);

set(0,'defaultfigurecolor','w')
figure('Position', [100, 100, 1880, 1560]);
subplot(4,6,[1 2 7 8 13 14 19 20])
scatter3(lon,lat,depth,power(exp(mag),1/5),'filled','MarkerEdgeColor',[153/255 0 0],'MarkerFaceColor',[255/255 51/255 51/255]);view(0,90)
hold on;box on
xlim([min(lon) max(lon)]);
ylim([min(lat) max(lat)]);
zlim([min(depth) max(depth)]);
set(gca,'ZDir','reverse','FontSize',20);

subplot(4,6,[3 4 9 10]);
histogram(mag);
set(gca,'fontsize',15)

hypotout=hypo;
fdata = hypo;
la0=min(lat);lo0=min(lon);
evla = fdata(:,8);evlo = fdata(:,9);evdp = fdata(:,10);
evlo_km = deg2km(distance([la0*ones(length(evla),1),evlo],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
evla_km = deg2km(distance([evla,lo0*ones(length(evla),1)],[la0*ones(length(evla),1),lo0*ones(length(evla),1)]));
hypo_out_points=[evlo_km,evla_km,evdp,mag,hypo(:,2),hypo(:,3),hypo(:,4),hypo(:,5),hypo(:,6),hypo(:,7)];

%% 3D Hough Transform 

[datapath,dirname]=fileparts(mfilename('fullpath'))
PATH=[datapath,'/SaveData']

path_parts=strsplit(datapath,'/');
h3d_idex=find(strcmp(path_parts,'hough-3d-lines-master'),1)
h3d_path=strjoin(path_parts(1:h3d_idex),'/')
filename_path=[num2str(h3d_path),'/data/fault_try.dat'];
fopen(filename_path,'w');
fopen([num2str(h3d_path),'/FaultSegment/Fault_predict.txt'],'w')

dlmwrite(num2str(filename_path), hypo_out_points,'precision', 8);
cd(num2str(h3d_path));
command = './hough3dlines ./data/fault_try.dat > ./FaultSegment/Fault_predict.txt';
system(command)  
f_sort=load([num2str(h3d_path),'/sp.output.txt'])

line_1=unique(f_sort(:,end));
for i=1:length(unique(f_sort(:,end)))
   line_nums_1(i,1)=find(ismember(f_sort(:,end),line_1(i,1),'rows'),1);
   line_nums_1(i,2)=find(ismember(f_sort(:,end),line_1(i,1),'rows'),1,'last')+1;
end
disp(['lines: ',num2str(line_nums_1(:,1)')])
ab_output=readtable([num2str(h3d_path),'/ab.output.txt']);
A=table2array(ab_output(:,1:3));
B=table2array(ab_output(:,4:6));

% DBSCAN  
HY=[];L=[[15 21];[16 22];[17 23];[18 24]];IN=[];HY_km=[];NEvents=size(hypo_out_points,1);
for j=1:4
    n=0;
    subplot(4,6,L(j,:))
    Points=hypo_out_points;
    if j==4
        MaxMag=max(Points(:,4));
        while 1
            Mag=max(Points(:,4));
            if Mag>2 && Mag<7
                k=5;
                radius=median(C(1,:))*0.001*power(10^(1.5*Mag+9.1)/(k^2*3*10^10/1000),1/3)/k;
            elseif Mag<2
                k=1;
                radius=median(C(1,:))*0.001*power(10^(1.5*Mag+9.1)/(k^2*3*10^10/1000),1/3)/k;
            elseif Mag>7
                k=10;
                radius=median(C(1,:))*0.001*power(10^(1.5*Mag+9.1)/(k^2*3*10^10/1000),1/3)/k;
            end
            dbs=dbscan(Points(:,1:3),radius,minPoi);
            if isempty(Mag) || all(dbs==-1)
                break;
            end
            for i=1:length(unique(dbs,'rows'))
                in=find(dbs==i);n=n+1;
                rs=rand(1,3);
                Points(in,11)=rs(1);Points(in,12)=rs(2);Points(in,13)=rs(3);Points(in,14)=n;
                HY_km=[HY_km;Points(in,:)];
                IN=[IN;in];       
            end
            Points(IN,:)=[];
            if isempty(Points)
                break
            end
        end
        HY=KMtoDEG(HY_km,la0,lo0);
        Rmag=round(radius_model(MaxMag,LowerMag,UpMag,MODEL)*median(C(1,:)), 3);
        Scatter3_E2F(hypotout,HY,NEvents,[' \in = ',num2str(Rmag),'km'],0,1);
    else
        dbs=dbscan(Points(:,1:3),R(1,j),minPoi);
        for i=1:length(unique(dbs,'rows'))
            rs=rand(1,3);
            in=find(dbs==i);n=n+1;
            Points(in,11)=rs(1);Points(in,12)=rs(2);Points(in,13)=rs(3);Points(in,14)=n;
            HY_km=[HY_km;Points(in,:)];
        end
        HY=KMtoDEG(HY_km,la0,lo0);
        Scatter3_E2F(hypotout,HY,NEvents,['\in =',num2str(R(1,j)),' km'],0,1);
    end
    HY_km=[];
end
 
sita=acos(sqrt(B(:,1).^2+B(:,2).^2))*180/pi;
fai=real(asin(B(:,2)./cos(sita./180*pi))*180/pi);
fai_m=median(fai);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
MAD=median(abs(fai-fai_m));
% 绘制sita直方图
sita=acos(sqrt(B(:,1).^2+B(:,2).^2))*180/pi;
subplot(4,6,[5 6]);histogram(real(sita),9,'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6],'LineWidth',2);
xlabel('sita(°)');ylabel('Fault lines');ax=gca;ax.FontSize=10;ax.Box="on";
% 绘制fai直方图
[counts, edges,number]=histcounts(fai,-90:10:90);
yfai=max(counts)+10;
faiup=fai_m+na*MAD;faidown=fai_m-na*MAD;
if faiup>90
    faiup=90;
end
if faidown<-90
    faidown=-90;
end
uplim=find(faiup<edges|faiup==edges,1)-1;
downlim=find(edges<faidown|faidown==edges==1,1,'last');
h=subplot(4,6,[11 12]);pos = get(h,'Position');
for i=1:length(counts(1,:))
    if i<downlim || uplim<i 
        bar(edges(i)+(edges(end)-edges(end-1))/2,counts(i),edges(end)-edges(end-1),'FaceColor',[1 0.2 0],'EdgeColor',[0.6 0 0],'LineWidth',2);hold on
    else 
        bar(edges(i)+(edges(end)-edges(end-1))/2,counts(i),edges(end)-edges(end-1),'FaceColor',[0 0.8 1],'EdgeColor',[0 0.2 0.6],'LineWidth',2);hold on
    end
end
acc_num=size(fai(faidown<=fai & fai<=faiup),1);
unacc_num=size(fai(fai<faidown | faiup<fai),1);
xticks([edges(1,1):edges(end)-edges(end-1):edges(1,end)]);ylim([0 yfai]);
xlabel('Azimuth(°)');ylabel('Fault lines');ax=gca;ax.FontSize=10;ax.Box="on";
annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 0.8*pos(3), 0.89*pos(4)], 'LineStyle','none',...
            'String', ['A=',num2str(acc_num)], ...
            'FontSize', 13, 'Color', [0 0.2 0.6]);
annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, 0.8*pos(3), 0.75*pos(4)], 'LineStyle','none',...
            'String', ['U=',num2str(unacc_num)], ...
            'FontSize', 13, 'Color', [0.6 0 0]);
% 
figure('Position', [100, 100, 1880, 1560]);
faia= fai(-90+10*downlim<=fai&fai<=-90+10*uplim);
faiu= fai(fai<-90+10*downlim|-90+10*uplim<fai);
% Hough for Accept and Unaccept
TIME=[];MEDIAN0=[];MEDIAN1=[];f_sort0=[];f_sort1=[];
for I=1:2
    if I==2,subplot(4,3,[1 4 7]);else subplot(4,3,[2 5 8]); end
    for i=1:line_1(end)
        SITA=acos(sqrt(B(i,1).^2+B(i,2).^2))*180/pi;
        FAI=real(asin(B(i,2)./cos(SITA./180*pi))*180/pi);
        LENGTH=line_nums_1(i,2)-line_nums_1(i,1);
        index=line_nums_1(i,1):line_nums_1(i,2)-1;
        f_sort(index,11:13)=rand(1,3).*ones(LENGTH,3);
        if (FAI<faidown || faiup<FAI) && I==1    
            MEDIAN0=[MEDIAN0;[median(f_sort(index,1:3)),f_sort(index(1),14),size(f_sort(index,1:3),1)]];
            f_sort0=[f_sort0;f_sort(index,:)];
        elseif faidown<=FAI && FAI<=faiup && I==2
            MEDIAN1=[MEDIAN1;[median(f_sort(index,1:3)),f_sort(index(1,1),14),size(f_sort(index,1:3),1)]];
            f_sort1=[f_sort1;f_sort(index,:)];
       end
    end
    if I==2
        f_sort1=KMtoDEG(f_sort1,la0,lo0);
        Scatter3_E2F(hypotout,f_sort1,NEvents,['Accepted'],1)
    else
        f_sort0=KMtoDEG(f_sort0,la0,lo0);
        Scatter3_E2F(hypotout,f_sort0,NEvents,['Unaccepted'],1);hold on
    end
end

subplot(4,3,[3 6 9])
uf=unique(f_sort0(:,14));f_sort_deg=KMtoDEG(f_sort,la0,lo0);
Ys=size(MEDIAN0,1)+size(MEDIAN1,1)+10;Ye=Ys/2;
disp('Hough Finished')
%% Hough Intergate
for i=1:size(MEDIAN0,1)
    guiyi1=max(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3)));
    guiyi2=max(MEDIAN0(i,5)./MEDIAN1(:,5));
    [d,idxmerge]=min(pdist2(MEDIAN0(i,1:3),MEDIAN1(:,1:3))./guiyi1);
    M1_cho=MEDIAN1(idxmerge,4);
    M1_color=find(f_sort(:,14)==M1_cho);
    rs=unique(f_sort(M1_color,11:13),'rows');
    M0_color=find(f_sort(:,14)==MEDIAN0(i,4));
    f_sort(M0_color,11)=rs(1,1);
    f_sort(M0_color,12)=rs(1,2);
    f_sort(M0_color,13)=rs(1,3);
    f_sort(M0_color,14)=M1_cho;
    MEDIAN1(find(MEDIAN1(:,4)==M1_cho),:)=[];
end

f_sort= sortrows(f_sort, size(f_sort, 2));
line_nums_2=[];
line_2=unique(f_sort(:,end));
for i=1:length(line_2)
   line_nums_2(i,1)=find(ismember(f_sort(:,end),line_2(i,1),'rows'),1);
   line_nums_2(i,2)=find(ismember(f_sort(:,end),line_2(i,1),'rows'),1,'last')+1;
   line_nums_2(i,3)=line_2(i,1);
end
disp('Hough&DBSCAN Finished')
F1=KMtoDEG(f_sort,la0,lo0);
Scatter3_E2F(hypotout,F1,NEvents,[' '],1);
subplot(4,3,10)
POLAR_gram(faia,0)
subplot(4,3,11)
POLAR_gram(faiu,1)
subplot(4,3,12)
POLAR_gram(faia,2,faiu)

hypotout2=[hypotout(:,9),hypotout(:,8),hypotout(:,10),hypotout(:,11),...
    hypotout(:,2:7),ones(size(hypotout(:,1),1),1).*[255/255 51/255 51/255 1]];
figure('Position', [100, 100, 1880, 1560]);
subplot(1,4,1)
Scatter3_E2F(hypotout,hypotout2,NEvents,[''],1)
subplot(1,4,2)
Scatter3_E2F(hypotout,f_sort1,NEvents,['Accepted'],1)
subplot(1,4,3)
Scatter3_E2F(hypotout,f_sort0,NEvents,['Unaccepted'],1);
subplot(1,4,4)
Scatter3_E2F(hypotout,F1,NEvents,[' '],1);

%% Segment clustering %%
figure('Position', [100, 100, 1880, 1560]);
colortemplate=rand(200,300);Kused=[];
for l=1:8
    subplot(4,4,[C(2,l),C(3,l)])
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
                c=C(1,l);
                radius=c*radius_model(Mag,LowerMag,UpMag,MODEL);

                [~,s4]=sort(ftest4(:,4),'descend');
                ftest4=ftest4(s4,:);

                idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
                nw=nw+1;
                if all(idx_5 == -1) || isempty(idx_5)
                    RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,0,];
                    break
                end
                if  MaxMag>CauMag
                    M1=M1+1;
                else
                    M2=M2+1;
                end
                RabY=[RabY;i,a,b,mean(Y(:,1)),size(ftest4,1),MaxMag,CauMag,Ns,Mag,1];
                Kused=[Kused;Mag,1,l];
                color_5=rand(length(unique(idx_5)),3);

                F3=ftest4(find(idx_5==1),:);
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
    FE=KMtoDEG(FE_km,la0,lo0);
    Scatter3_E2F(hypotout,FE,NEvents,['C= ',num2str(C(1,l))],0,l);
    disp([' M1= ',num2str(M1),' M2= ',num2str(M2)])
    RC=[RC;[M1,M2]];
end
RC_per=sum(RC,2)./line_1(end);
disp('Diff C Hough&DBSCAN Finished');

%% Elliptic Fitting %%
figure('Position', [100, 100, 1880, 1560]);number_EPC=[];number_Eve=[];ECC=[];
ratio=[1;1;1];Azimuth=[];Elevation=[];FaultpaRameters=[];Magrecord=[];
for l=1:8  
    subplot(4,4,[C(2,l),C(3,l)])
    [FE_km,FE]=zdHouDBS(C,l,f_sort,line_1,line_nums_2,la0,lo0,minPoi,Mc,hypotout,colortemplate,LowerMag,UpMag,MODEL);
    EP=[];EP_km=[];X=[];Y=[];Z=[];n=0;ell=[];
    line_e=unique(FE(:,end));
    for i=1:length(line_e)
       line_nums_e(i,1)=find(ismember(FE(:,end),line_e(i,1),'rows'),1);
       line_nums_e(i,2)=find(ismember(FE(:,end),line_e(i,1),'rows'),1,'last')+1;
       line_nums_e(i,3)=line_e(i,1);
    end
    for i=1:length(line_e)
        sp_km=FE_km(line_nums_e(i,1):line_nums_e(i,2)-1,:);

        w=power(10.^(11.8+1.5*sp_km(:,4))/max(10.^(11.8+1.5*sp_km(:,4))),1/3);
        L=w*ratio(1,1);W=w*ratio(2,1);H=w*ratio(3,1);weight=[L,W,H];
        m1 = mean(sp_km(:,1:3));
        C1 = cov(sp_km(:,1:3));
        loca=median(hypo_out_points(:,1:3));
        [a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,el,linea,Vol,elliptic_points]=plotcov_3d(C1, m1,sp_km,weight,ratio,0,s,loca);
        n=n+1;
        ell=[ell;[elliptic_points,ones(size(elliptic_points,1),1).*n]];
        Azimuth=[Azimuth;[az,l]];
        Elevation=[Elevation;[el,l]];
        mass=KMtoDEG(sp_km,la0,lo0);
        Magrecord=[Magrecord;[sp_km(:,4),ones(size(sp_km,1),1).*l,ones(size(sp_km,1),1).*n]];
        FaultpaRameters=[FaultpaRameters;[median(mass(:,1:3)),a,b,c,p1_u,p1_d,p2_u,p2_d,p3_u,p3_d,az,el,linea,size(sp_km,1)/Vol,l,size(sp_km,1)]];
        EP_km=[EP_km;elliptic_points,sp_km(1,11:13).*ones(size(elliptic_points,1),3)];
        ECC=[ECC;[linea,l,size(sp_km,1)/Vol,length(sp_km(:,1))]];
    end
    Scatter3_E2F(hypotout,FE,NEvents,['C= ',num2str(C(1,l))],0,l);hold on
    ellsurf(ell,la0,lo0);

    EP=KMtoDEG(EP_km,la0,lo0);  
    number_EPC=[number_EPC;length(unique(EP(:,4:6),'rows'))];
    number_Eve=[number_Eve;length(FE(:,1))];
end
FaultpaRameters(FaultpaRameters(:,28)>1,28)=1;
ECC(ECC(:,3)>1,3)=1;
ratio_Eve=number_Eve./NEvents;
disp('Fitting Finished')
 
%% Subsurface fault
L=4;
SUBfault_con(L,C,f_sort,line_1,line_nums_2,hypotout,FaultpaRameters,Mc,LowerMag,UpMag,MODEL,minPoi)

%% RSM
figure('Position', [100, 100, 1880, 1560]);Energyratio=[];
mi=2;
for j=1:8
    subplot(4,4,[C(2,j),C(3,j)]);
    SubRLength=[];SubRWidth=[];Moreal=[];k=[];
    num0=find(FaultpaRameters(:,29)==j);
    FR=FaultpaRameters(num0,:);
    Ku=Kused(num0,:);
    for i=1:size(FR,1)
        dist1=pdist2([FR(i,7:9)],[FR(i,10:12)]);
        dist2=pdist2([FR(i,13:15)],[FR(i,16:18)]);
        SubRLength=[SubRLength;dist1];
        SubRWidth=[SubRWidth;dist2];
        if LowerMag<=Ku(i,1) & Ku(i,1)<UpMag
            k=[k;5];
        elseif Ku(i,1)<LowerMag
            k=[k;1];
        elseif Ku(i,1)>=UpMag
            k=[k;10];
        end
    end
    Motheo=3*10^10.*(10^3.*SubRLength).*(10^3.*SubRWidth).*SubRLength;
    num1=find(Magrecord(:,2)==j);
    MR=Magrecord(num1,:);
    for i1=1:MR(end,3)       
        num2=find(MR(:,3)==i1);
        Moreal=[Moreal;sum(10.^(1.5.*MR(num2,1)+9.1))];
    end
    loglog(Moreal,Moreal,'LineWidth',2,'LineStyle','-','Color',[0 0 0]);hold on
    d1=min(Moreal);
    d2=max(Moreal);
    dx=(d2-d1)/20;
    loglog(d1:dx:d2,10^mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle','--','Color',[0 0 0]);
    loglog(d1:dx:d2,10^-mi.*(d1:dx:d2),'LineWidth',1.5,'LineStyle',':','Color',[0 0 0]);
    ax=gca;
    set(ax,'FontSize',15,'XMinorTick','on','YMinorTick','off','TickLength',[0.02 0.1],'LineWidth',1);box on
    loglog(Moreal ,Motheo,'+','MarkerSize',8,'LineWidth',1.5);
    legend('RSM=1',['RSM=',num2str(1+mi/10)],['RSM=',num2str(1-mi/10)],'location','southeast','box','off');
    Energyratio=[Energyratio;log10(Motheo)./log10(Moreal),j.*ones(size(FR,1),1)];
    pos = get(ax, 'Position');
    annotation('textbox', [pos(1)+pos(3)*0.08, pos(2)+pos(4)*0.12, pos(3), 0.8*pos(4)], 'LineStyle','none',...
        'String', ['C= ',num2str(C(1,j))], ...
        'FontAngle','normal','FontSize', 20, 'Color', 'black');
 xlim([10^7 10^15])
 ylim([10^2 10^20])
 xticks([10^7 10^11 10^15])
 yticks([10^2 10^11 10^20])
end

%% Azimuth
figure('Position', [100, 100, 1880, 1560]);
for l=1:8  
    subplot(4,4,[C(2,l),C(3,l)])
    Azimuth=FaultpaRameters(FaultpaRameters(:,29)==l,25);
    AzimuthStat([Azimuth;Azimuth-180]);
    title(['C=',num2str(C(1,l))],'FontWeight','normal');
    rlim([0 20]);
end
%% Elevation
figure('Position', [100, 100, 1880, 1560]);
for l=1:8  
    subplot(4,4,[C(2,l),C(3,l)])
    DipA=FaultpaRameters(FaultpaRameters(:,29)==l,26);
    ElevationStat(DipA);
    title(['C=',num2str(C(1,l))],'FontWeight','normal');
    rlim([0 30]); 
end
%% Evalulate
EvaLF(ECC,Energyratio,NEvents,C,(10-mi)/10,2)
l61=[];l62=[];
for i=1:8
    numi=find(ECC(:,2)==i);
    l61=[l61;length(find(ECC(numi,1)>0.3))];
    numj=find(Energyratio(:,2)==i);
    l62=[l62;length(find(Energyratio(numj,1)>0.6))];
end
