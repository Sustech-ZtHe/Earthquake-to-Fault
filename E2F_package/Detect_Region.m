function [R,C,scale,minPoi,Mc,na,hypo,hypo_out]=Detect_Region(name)
    if strcmp(name, 'Toc2me')
       hypo_out = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/ToC2ME_test.txt');
    elseif strcmp(name, 'chuannan')
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/chuannan/chuanan.txt');
    elseif strcmp(name, 'foxcreek')  
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/Foxcreek/Foxcreek.txt');
    elseif strcmp(name, 'peaceriver')  
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/peaceriver/peaceriver.txt');
    elseif strcmp(name, 'ratonbasin') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/ratonbasin/Ratonbasin.txt');
    elseif strcmp(name, 'reddeer') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/reddeer/Catalog_Box.txt');
    elseif strcmp(name, 'weiyuan')  
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/weiyuan/weiyuan.txt');
    elseif strcmp(name, 'turkey_afad')   
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/turkey/TURKEY_test.txt');
    elseif strcmp(name, 'turkey_usgs') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/turkey/tk_cata_usgs.txt');
    elseif strcmp(name, 'turkey_ml') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/turkey/TK_ml.txt');
    elseif strcmp(name, 'qinghai') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/qinghai/qinghai.txt');
    elseif strcmp(name, 'ridgecrest') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/ridgecrest/ridgecrest.txt');
    elseif strcmp(name, 'wenchuan') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/wenchuan/wenchuan.txt');
    elseif strcmp(name, 'SouthernC_12dev') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/SouthernCalifornia/SouthernC_12dev.txt');
    elseif strcmp(name, 'SouthernC_9.5dev') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/SouthernCalifornia/SouthernC_9.5dev.txt');
    elseif strcmp(name, 'Southernc2010') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/SouthernCalifornia_loc/Southernc2010.txt');
    elseif strcmp(name, 'Pawnee') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/Pawnee/Pawnee.txt');
    elseif strcmp(name, 'Grant') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/Grantcounty/Grantcounty.txt');
    elseif strcmp(name, 'Okamc') 
       hypo = readtable('/home/me/Documents/H3D/hough-3d-lines-master/OriData/Magcatalog/Okamc/Okamc.txt');
    end
    if strcmp(name, 'Toc2me') 
        hypo=hypo_out;
        hypo.Var1=hypo_out.Var1;
        hypo.Var2=hypo_out.Var5;
        hypo.Var3=hypo_out.Var6;
        hypo.Var4=hypo_out.Var7;
        hypo.Var5=hypo_out.Var8;
        hypo.Var6=hypo_out.Var9;
        hypo.Var7=hypo_out.Var10;
        hypo.Var8=hypo_out.Var2;
        hypo.Var9=hypo_out.Var3;
        hypo.Var10=hypo_out.Var4;
        hypo.Var11=hypo_out.Var11;
    elseif strcmp(name, 'foxcreek') || strcmp(name, 'peaceriver')... 
                                    || strcmp(name, 'ratonbasin')...
                                    || strcmp(name, 'reddeer')... 
                                    || strcmp(name, 'weiyuan')...  
                                    || strcmp(name, 'turkey_afad')...  
                                    || strcmp(name, 'turkey_usgs')... 
                                    || strcmp(name, 'qinghai')...
                                    || strcmp(name, 'chuannan')...
                                    || strcmp(name, 'ridgecrest')...
                                    || strcmp(name, 'wenchuan')...
                                    || strcmp(name, 'turkey_ml')...
                                    || strcmp(name, 'SouthernC_12dev')...
                                    || strcmp(name, 'SouthernC_9.5dev')...
                                    || strcmp(name, 'Southernc2010')...
                                    || strcmp(name, 'Pawnee')...
                                    || strcmp(name, 'Okamc')...
                                    || strcmp(name, 'Grant')
        hypo_out=hypo;
        hypo_out.Var2=hypo.Var8;
        hypo_out.Var3=hypo.Var9;
        hypo_out.Var4=hypo.Var10;
        hypo_out.Var5=hypo.Var2;
        hypo_out.Var6=hypo.Var3;
        hypo_out.Var7=hypo.Var4;
        hypo_out.Var8=hypo.Var5;
        hypo_out.Var9=hypo.Var6;
        hypo_out.Var10=hypo.Var7;
    end

    [counts, edges]=histcounts(table2array(hypo_out(:,11)));
    [maxCount, maxIndex] = max(counts);
    Mc=edges(maxIndex + 1);
    na=2;minPoi=10;
    if max(table2array(hypo_out(:,3)))-min(table2array(hypo_out(:,3)))<=0.1
        R=[0.01,0.1,0.2];C=[1,2,3,4,5,6,7,8;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=0.1;
    elseif 0.1<max(table2array(hypo_out(:,3)))-min(table2array(hypo_out(:,3)))&&...
        max(table2array(hypo_out(:,3)))-min(table2array(hypo_out(:,3)))<=1
        R=[0.1,0.5,1];C=[5,8,10,12,15,18,20,25;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=1;
    else
%         R=[3,6,10];C=[13,14,15,16,17,18,19,20;
%                                        1,2,3,4,9,10,11,12;
%                                        5,6,7,8,13,14,15,16];scale=10;
        R=[3,6,10];C=[15,20,25,30,35,40,45,50;
                                       1,2,3,4,9,10,11,12;
                                       5,6,7,8,13,14,15,16];scale=50;
    end
end

% well1=table2array(readtable('/home/me/Documents/HOU3D-CODE/well-1.txt'));
% well2=table2array(readtable('/home/me/Documents/HOU3D-CODE/well-2.txt'));
% well3=table2array(readtable('/home/me/Documents/HOU3D-CODE/well-3.txt'));
% well4=table2array(readtable('/home/me/Documents/HOU3D-CODE/well-4.txt'));

% % figure;scatter3(table2array(hypo_out(:,3)),table2array(hypo_out(:,2)),table2array(hypo_out(:,4)),'filled','MarkerEdgeColor',[0.29 0.43 0.54],'MarkerFaceColor',[0.42 0.65 0.8]);
% % ax=gca;set(ax,'Zdir','reverse');grid on;grid minor;box on;view(-50,0);zlim([2.5 4.5])
% % histogram(hypo_out.Var11);

% figure;scatter(table2array(hypo_out(:,3)),table2array(hypo_out(:,2)),table2array(hypo_out(:,11))) 
% cc=table2array(hypo_out)
% n=find(cc(:,11)==max(cc(:,11)));
% cc(n,:)
% hold on
% scatter(cc(n,3),cc(n,2),cc(n,11),"red") 

% figure;scatter3(table2array(hypo_out(:,3)),table2array(hypo_out(:,2)),table2array(hypo_out(:,4))) 