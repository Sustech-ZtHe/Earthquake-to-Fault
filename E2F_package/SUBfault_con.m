function Rectangle_points=SUBfault_con(Optimal_C,f_sort_E2F,line_nums_2,hypotout,FaultpaRameters,minPoi,colortemplate)
n=1;FE_km=[];c=Optimal_C;
for i=1:f_sort_E2F(end)
    im=1;
    if ismember(i,unique(f_sort_E2F(:,14)))
        k=find(line_nums_2(:,3)==i);
        ftest4=f_sort_E2F(line_nums_2(k,1):line_nums_2(k,2),:);
        while 1
            MaxMag=max(ftest4(:,4));
            if length(ftest4)>100
                [a,b,CauMag]=Mmain(ftest4,hypotout);
                Mag=max(MaxMag,CauMag);
            else
                Mag=MaxMag;
            end
            radius=c*radius_model(Mag);
            [~,s4]=sort(ftest4(:,4),'descend');
            ftest4=ftest4(s4,:);
            idx_5=dbscan(ftest4(:,1:3),radius,minPoi);
            if all(idx_5 == -1) || isempty(idx_5)
                break
            end            
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
else
    % FaultpaRameters(:,29)
    % Optimal_C
    Rectangle_points=[];
    l=find(FaultpaRameters(:,28)==Optimal_C);
    % l=find(FaultpaRameters(:,29)==1);
    % for i=1:size(l)
    %     x=[FaultpaRameters(l(i),7)',FaultpaRameters(l(i),7)',FaultpaRameters(l(i),10)',FaultpaRameters(l(i),10)'];
    %     y=[FaultpaRameters(l(i),8)',FaultpaRameters(l(i),8)',FaultpaRameters(l(i),11)',FaultpaRameters(l(i),11)'];
    %     z=[FaultpaRameters(l(i),15)',FaultpaRameters(l(i),18)',FaultpaRameters(l(i),18)',FaultpaRameters(l(i),15)'];
    %     fill3(x, y, z, [96 96 96]./255);hold on
    % end
    % scatter3(FE_km(:,1),FE_km(:,2),FE_km(:,3),10,FE_km(:,11:13),'filled')
    % set(gca,'ZDir','reverse','Color',[0.9 0.9 0.9],'FontSize',15);box on;view(-30,60)
    % xlabel('Easting(km)');ylabel('Northing(km)');zlabel('Depth(km)')
    for i=1:size(l)
        A1=[FaultpaRameters(l(i),7),FaultpaRameters(l(i),8),FaultpaRameters(l(i),9)];
        A2=[FaultpaRameters(l(i),10),FaultpaRameters(l(i),11),FaultpaRameters(l(i),12)];
        B1=[FaultpaRameters(l(i),13),FaultpaRameters(l(i),14),FaultpaRameters(l(i),15)];
        B2=[FaultpaRameters(l(i),16),FaultpaRameters(l(i),17),FaultpaRameters(l(i),18)];
        L = A2 - A1; 
        S = B2 - B1; 
        L_unit = L / norm(L); 
        S_unit = S / norm(S); 
        S_half = 0.5 * norm(S);
        V1 = A1 - S_half * S_unit; 
        V2 = A1 + S_half * S_unit; 
        V3 = A2 + S_half * S_unit; 
        V4 = A2 - S_half * S_unit; 
        rectangle_points = [V1; V2; V3; V4; V1]; 
        % rectangle_points_deg=KMtoDEG(rectangle_points,la0,lo0);
        Rectangle_points=[Rectangle_points,[rectangle_points,i.*ones(length(rectangle_points),1)]];
        % plot3(rectangle_points_deg(:,1), rectangle_points_deg(:,2), rectangle_points_deg(:,3), 'k-', 'LineWidth', 2);
        plot3(rectangle_points(:,1), rectangle_points(:,2), rectangle_points(:,3), 'k-', 'LineWidth', 2);hold on;
    end
    grid on
    % num1=find(FE_km(:,14)==4);
    % FE_km(num1,11:13)=[187 121 186]./255.*ones(length(num1),3);
    % num2=find(FE_km(:,14)==3);
    % FE_km(num2,11:13)=[86 118 168]./255.*ones(length(num2),3);
    % FE_km(:,14)
    % length(num1)
    % length(num2)

    scatter3(FE_km(:,1),FE_km(:,2),FE_km(:,3),10,FE_km(:,11:13),'filled');axis equal;
    set(gca,'ZDir','reverse','Color',[0.9 0.9 0.9],'FontSize',15,'Box','on','LineWidth',1.5);box on;view(-30,60);
    xlabel('Easting(km)');ylabel('Northing(km)');zlabel('Depth(km)'); 
    annotation('textbox', [0.4 0.01 0.7 0.4], 'String', ['Fault candidates=',num2str(FE_km(end,14)),newline,'Clustered events=',num2str(size(FE_km,1)),newline,'Remaining events=',num2str(size(hypotout,1)-size(FE_km,1))], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

    % FE=KMtoDEG(FE_km,la0,lo0);
    % scatter3(FE(:,1),FE(:,2),FE(:,3),10,FE(:,11:13),'filled');axis equal;
    % set(gca,'ZDir','reverse','Color',[0.9 0.9 0.9],'FontSize',15,'Box','on','LineWidth',1.5);box on;view(-30,60);
    % xlabel('Longtitude ');ylabel('Latitude ');zlabel('Depth(km)'); axis equal; view(0,90)
end