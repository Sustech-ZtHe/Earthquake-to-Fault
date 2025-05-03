function [MaxMw,FE_km,FE,F3rsm]=zdHouDBS(C,l,f_sort,line_nums_2,la0,lo0,minPoi,hypotout,colortemplate)
n=1;FE_km=[];c=C(1,l);MaxMw=[];F3rsm=[];
for i=1:f_sort(end)
    im=1;
    if ismember(i,line_nums_2(:,3))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort(line_nums_2(k,1):line_nums_2(k,2),:);
        while 1
            Mobs=max(ftest4(:,4));
            if length(ftest4)>100
                [a,b,Mthre]=Mmain(ftest4,hypotout);
                Mw=max(Mobs,Mthre);
            else
                Mw=Mobs;
            end
            radius=c*radius_model(Mw);            
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);
            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            if all(idx_5 == -1) || isempty(idx_5)
                break
            end
            MaxMw=[MaxMw;Mw];
            F3=ftest4(find(idx_5==1),:);
            F3(:,11)=colortemplate(i,3*(n-1)+1);F3(:,12)=colortemplate(i,3*(n-1)+2);F3(:,13)=colortemplate(i,3*(n-1)+3);
            F3(:,14)=n;
            FE_km=[FE_km;F3];
            F3rsm=[F3rsm;F3,n.*ones(size(F3,1),1),C.*ones(size(F3,1),1)];
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


