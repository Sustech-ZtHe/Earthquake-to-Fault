function [FE_km,FE]=zdHouDBS(C,l,f_sort,line_nums_2,la0,lo0,minPoi,Mc,hypotout,colortemplate,LowerMag,UpMag,MODEL)
n=1;M1=0;M2=0;FE_km=[];

for i=1:f_sort(end)
    im=1;jj=1;
    if ismember(i,unique(f_sort(:,14)))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort(line_nums_2(k,1):line_nums_2(k,2)-1,:);
        nw=1;
        while 1
            MaxMag=max(ftest4(:,4));
            if im==1
                [a,b,Y,CauMag]=Mmain(ftest4,Mc,nw,i,hypotout);
                im=2;
            elseif im==2
                CauMag=(log10(size(ftest4,1))-a)/b+Mc;
            end
            Mag=max(MaxMag,CauMag);
            if  MaxMag>CauMag
                M1=M1+1;
            else
                M2=M2+1;
            end
            c=C(1,l);
            radius=c*radius_model(Mag,LowerMag,UpMag,MODEL);            
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);

            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            nw=1;
            if all(idx_5 == -1) || isempty(idx_5)
                break
            end
            color_5=rand(length(unique(idx_5)),3);
            F3=ftest4(find(idx_5==1),:);
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


