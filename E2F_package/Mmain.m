function [a,b,Y,M]=Mmain(ftest4,Mc,nw,I,hypotout)
    if size(ftest4,1)<10 || size(find(ftest4(:,4)>Mc),1)<10
%         M=max(ftest4(:,4))-0.1;
        M=Mc;
        a=nan;b=nan;Y=nan;
%         scatter(ftest4(:,1), ftest4(:,2), 300,'k.');hold on
%         title(['n= ',num2str(size(ftest4,1)),' option 1,I=',num2str(I),'  nw=',num2str(nw),'  M=',num2str(M)])

    else
%         dt=(max(ftest4(:,4))-Mc)/size(ftest4,1);
        dt=(max(ftest4(:,4))-min(ftest4(:,4)))/size(ftest4,1);
        % ddt=t/size(ftest4,1);
        Y=[];Y2=[];
        if size(ftest4,1)>500000
            n=1:size(ftest4,1);
            st=size(ftest4,1)/5;
            for i=1:5
                n1=(i-1)*st+1:i*st;
                Y2=[Y2;sum(ftest4(:,4)'>(Mc+dt.*n1)',2)];
            end
            Y=[(Mc+dt.*n)',Y2];
        else
            Y(:,1)=Mc+dt.*(1:size(ftest4,1));
            Y(:,2)=sum(ftest4(:,4)>Mc+dt.*(1:size(ftest4,1)))';
        end

%         Z=[];
%         Z(:,1)=Mc+dt.*(1:size(ftest4,1));
%         Z(:,2)=sum(ftest4(:,4)<Mc+dt.*(1:size(ftest4,1)))';

% figure;histogram(ftest4(:,4))
%         for i=1:size(ftest4,1)
%             Y(i,1)=Mc+dt*(i-1);
%             Y(i,2)=size(ftest4(ftest4(:,4)>Mc+dt*(i-1),:),1);
% %             Y(i,1)=min(ftest4(:,4))+dt*(i-1);
% %             Y(i,2)=size(ftest4(ftest4(:,4)>min(ftest4(:,4))+dt*(i-1),:),1);
%         % Y(i,3)=0+ddt*(i-1);
%         end

        % 定义拟合方程
        % fun = @(params, x) 10.^(params(1)+params(3).* (params(2) +Y(:,1)))./((0.01+Y(:,3)).^1.08);
        fun = @(params, x) 10.^(params(1)+params(3).* (params(2) -Y(:,1)));
        % 初始参数猜测值
        params0 = [-1, Mc+1, 1];
%         la=[-inf,Mc,0];
%         lb=[0,max(hypotout(:,11))+1.2,2];
        la=[-inf,Mc,0];
        lb=[2,max(hypotout(:,11))+1.2,2];
        % 进行拟合
        params = lsqcurvefit(fun, params0, Y(:,1), Y(:,2),la,lb);
%         params = lsqcurvefit(fun, params0, Y(:,1), Y(:,2));
        % 提取拟合参数
        a = params(1);
        M = params(2);
        b = params(3);

%         figure
% %         plot((log10(1:Z(end,2))-a)/b+Mc,1:Z(end,2))
% 
%         scatter(Y(:,1), Y(:,2), 300,'k.');hold on
%         plot(Y(:,1),10.^(a+b.* (M -Y(:,1))),'Color',[1 0 0],'LineWidth',2)
%          scatter(Z(:,1), Z(:,2), 300,'k.');hold on
%          plot(Z(:,1),10.^(a+b.* Z(:,1)-Mc),'Color',[1 0 0],'LineWidth',2)

%         figure
%         semilogy((log10(1:Z(end,2))-a)/b+Mc,1:Z(end,2),'Color',[153/255 0 0],'LineWidth',4);hold on   
%         fill([(log10(1)-a)/b+Mc (log10(Z(end,2))-a)/b+Mc (log10(Z(end,2))-a)/b+Mc (log10(1)-a)/b+Mc], ...
%                                         [1 1 Z(end,2) 1], 'r', 'FaceAlpha', 0.3,'LineStyle','none')
%         xlim([(log10(1)-a)/b+Mc (log10(Z(end,2))-a)/b+Mc])
%         ax=gca;
%         ax.FontSize=18;ax.XMinorTick="on";ax.YMinorTick="on";ax.ZMinorTick="on";
%         ax.TickLength=[0.02 0.1];grid minor;box on
%         xlabel('Theoretical magnitude');ylabel('Numbers of event');

%         title(['n= ',num2str(size(ftest4,1)),' option 2,I=',num2str(I),'  nw=',num2str(nw),'  M=',num2str(M),'  a=',num2str(a),'  b=',num2str(b)])
%         title('M_{main}=',num2str(M))
%         legend('Event','Fit')
%          figure
      
%         pos = get(ax, 'Position');
%         annotation('textbox', [pos(1)+pos(3)*0.55, pos(2)+pos(4)*0.1, 0.1*pos(3), 0.55*pos(4)], 'LineStyle','none',...
%             'String', ['M_{main}=',num2str(M),'  ', char(9679)+string('Events ')], ...
%             'FontWeight','bold','FontSize', 12, 'Color', [0 0 0]); 
    end
end

% i=1
% (log10(5000)-RabY(i,1))/RabY(i,2)+RabY(i,3)
% figure
% for i=1:5
%     x=1:40000;
%     plot(x,(log10(x)-RabY(i,2))/RabY(i,3)+RabY(i,4),'linewidth',2);hold on
%     drawnow
%     pause(1)
% end

% figure;histogram(Y(:,1))
% dt=(max(ftest4(:,4))-Mc)/size(ftest4,1);
% Y=[];
% for i=1:size(ftest4,1)
%     Y(i,1)=Mc+dt*(i-1);
%     Y(i,2)=size(ftest4(ftest4(:,4)>Mc+dt*(i-1),:),1);
% end
% % 定义拟合方程
% fun = @(params, x) exp(params(1).* (params(2) - Y(:,1)));
% % 初始参数猜测值
% params0 = [1, 1];
% % 进行拟合
% params = lsqcurvefit(fun, params0, Y(:,1), Y(:,2));
% % 提取拟合参数
% a = params(1)
% M = params(2)


