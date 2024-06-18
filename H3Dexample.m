clc;clear all
addpath 'D:\20230318\hough-3d-lines-master\data';
data_ori=readtable('D:\20230318\hough-3d-lines-master\data\synthetic_c.dat');
data_test=table2array(data_ori);
data_test=data_test(data_test(:,3)<0&data_test(:,1)<3,:);

figure('Position', [600, 400, 500, 500]);
% scatter3(data_test(:,1),data_test(:,2),data_test(:,3),50,'MarkerFaceColor',[.2 .6 .8],'MarkerfaceAlpha',0.4,'MarkerEdgeColor','none')
scatter3(data_test(:,1),data_test(:,2),data_test(:,3),50,'black','filled','MarkerFaceAlpha',0.4)
ax=gca;set(ax,'FontSize',15);grid on;grid minor;box on;view(-45,50)
xticks(-2:1:4);yticks(-5:2:5);zticks(-4:1:0);
data_test=[data_test,ones(126,8)];
% save('data_test.mat',"data_test");
data_testf1=fopen('data_test_txt.txt','w');
fprintf(data_testf1,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',data_test.');
fclose(data_testf1);