function [MaxMw,FE_km,FE]=zdHouDBS(C,l,f_sort,line_nums_2,la0,lo0,minPoi,Mc,hypotout,colortemplate)
n=1;FE_km=[];c=C(1,l);MaxMw=[];
for i=1:f_sort(end)
    im=1;
    if ismember(i,line_nums_2(:,3))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort(line_nums_2(k,1):line_nums_2(k,2),:);
        while 1
            MaxMag=max(ftest4(:,4));
            if im==1
                [a,b,Y,CauMag]=Mmain(ftest4,Mc,hypotout);
                im=2;
            elseif im==2
                CauMag=(log10(size(ftest4,1))-a)/b+Mc;
            end
            Mag=max(MaxMag,CauMag);
            radius=c*radius_model(Mag);            
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);
            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            if all(idx_5 == -1) || isempty(idx_5)
                break
            end
            MaxMw=[MaxMw;Mag];
            F3=ftest4(find(idx_5==1),:);
            F3(:,11)=colortemplate(i,3*(n-1)+1);F3(:,12)=colortemplate(i,3*(n-1)+2);F3(:,13)=colortemplate(i,3*(n-1)+3);
            F3(:,14)=n;
            FE_km=[FE_km;F3];
            n=n+1;
            ftest4(idx_5==1,:)=[];
            if isempty(ftest4)
                break
            end
        end
    end
end
if isempty(FE_km)
    disp('C is too small, please increase it !')
    FE=[];
else
    FE=KMtoDEG(FE_km,la0,lo0);
end
end


