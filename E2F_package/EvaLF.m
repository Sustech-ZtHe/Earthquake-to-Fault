function EvaLF(ECC,Energyratio,NEvents,C,up,style)
if style==1
    figure('Position', [100, 100, 1880, 1560]);sumevents=[];
    t=tiledlayout(7,1)
    for j=1:8
       numse=find(ECC(:,2)==j);
       subECC=ECC(numse);
       sumevents=[sumevents;sum(ECC(numse,4))];
       subER=Energyratio(numse);
       subER0=subER(subER<=1);
       subER1=[size(subER(subER>1),1),0,0,0,0,0,0,0,0,0];
       [countse,edgese,numberse]=histcounts(subECC,[0:0.1:1]);
       [counter0,edgeer0,numberer0]=histcounts(subER0,[0:0.1:1]);
       countse=flip(countse);
       counter0=flip(counter0);
       edge=edgese(j):0.01:edgese(j+1);
       subplot(7,1,[1 2 3])
       for i=1:8
           bar(edge(1:end-1)+0.005,countse,1,'FaceColor',[1 0.6 0.2],'EdgeColor',[0.4 0.2 0],'LineWidth',2);hold on;box off
           xlim([0 0.8]);
           set(gca,'XTickLabel','','FontSize',15,'YMinorTick','on','LineWidth',2,'TickLength',[0.005 0.05], ...
               'YColor',[0.4 0.2 0],'XColor',[0.4 0.2 0]);
       end
       subplot(7,1,[5 6 7])
       for i=1:8
    %        figure
           bar(edge(1:end-1)+0.005,counter0,1,'FaceColor',[0 1 0.4],'EdgeColor',[0 0.4 0.2],'LineWidth',2);hold on;
           bar(edge(1)+0.005,subER1(1),0.01,'FaceColor',[1 0 0],'EdgeColor',[204 0 0]./255,'LineWidth',2);
           xlim([0 0.8]);
           set(gca, 'YDir', 'reverse','XTickLabel','','FontSize',15,'YMinorTick','on','LineWidth',2,'TickLength',[0.005 0.05], ...
               'YColor',[0 0.4 0.2],'XColor',[0 0.4 0.2],'xtick','');box off
       end   
    end
    subplot(7,1,4)
    plot(0:0.01:0.8,0,'.','LineWidth',2,'Color',[0 0 0])
    set(gca,'YTickLabel','','FontSize',12,'XMinorTick','on','LineWidth',2, ...
        'XAxisLocation','origin','YAxisLocation','origin','TickDir','both','TickLength',[0.005 0.05],'YTick','');box off
    xticks([0:0.05:0.8]);
    AX=gca;
    pos = get(AX, 'Position');
    AX.XTickLabel={'1','0.5', ...
                          '1','0.5', ...
                          '1','0.5', ...
                          '1','0.5', ...
                          '1','0.5', ...
                          '1','0.5', ...
                          '1','0.5', ...
                         '1','0.5',}
    for i=1:8
        text(0.035+(i-1)*0.1, 0.6, ['C=',num2str(C(1,i))], 'FontSize', 15,'Color', 'black');
    end
    ylabel('Fault lines','FontSize',20)
    subplot(7,1,[1 2 3])
    ax1=gca;
    pos = get(ax1, 'Position');
    axes('position',pos,'color','none','YColor',[0 102 204]./255,'yaxislocation', ...
        'right','YLim',[0 100],'FontSize',15,'XTickLabel','','XTick','','YMinorTick','on','LineWidth',2);hold on
    plot(0.05:0.1:0.75,sumevents./NEvents.*100,'-o','LineWidth',2,'Color',[0 102 204]./255,'MarkerFaceColor',[0 1 1])
    ylabel('Percentage (%)','FontSize',20)
    xlim([0 0.8])
elseif style==2
    figure('Position', [100, 100, 1880, 1560]);sumevents=[];sorx=[];
    for j=1:8
       numse=find(ECC(:,2)==j);
       subECC=ECC(numse);
       sumevents=[sumevents;sum(ECC(numse,4))];
       subER=Energyratio(numse);
       [countse,edgese,numberse]=histcounts(subECC,[0,0.6,1]);
       [counter0,edgeer0,numberer0]=histcounts(subER,[0,up,1,1-up+1]);
       countse=flip(countse);
       counter0=flip(counter0);
       edge=[0,0.05,0.1]+(j-1)*0.2;

       bar(edge(1)-0.025,countse(1),0.05,'FaceColor',[203 238 249]./255,'EdgeColor',[66 146 197]./255,'LineWidth',5);hold on;
       bar(edge(2)-0.025,counter0(1)+counter0(2),0.05,'FaceColor',[161 217 156]./255,'EdgeColor',[5 165 158]./255,'LineWidth',5)
       set(gca,'XTickLabel',{['C=',num2str(C(1,1))],['C=',num2str(C(1,2))],['C=',num2str(C(1,3))],['C=',num2str(C(1,4))], ...
                            ['C=',num2str(C(1,5))],['C=',num2str(C(1,6))],['C=',num2str(C(1,7))],['C=',num2str(C(1,8))],''},'FontSize',30,'YMinorTick','on','LineWidth',2,'TickLength',[0.005 0.05], ...
           'YColor',[0 0 0],'XColor',[0 0 0]);

       yup=max(counter0(1)+counter0(2),countse(2));

       xlim([-0.075 1.475]);
       ylim([0 yup*1.5]);
       ylabel('Number of faults','FontSize',30);
       box off
    end
    h1=legend('Aspect ratio','RSM','box','off')
    h2=get(h1,'Position')
    h1.Position=[h2(1)-0.05, h2(2)-0.05, 0.9*h2(3), 0.01*h2(4)]
    ax1=gca;
    pos = get(ax1, 'Position');
    axes('position',pos,'color','none','YColor',[143 38 126]./255,'yaxislocation', ...
        'right','YLim',[0 100],'FontSize',30,'XTickLabel','', ...
            'XTick','','YMinorTick','on','LineWidth',2);hold on
    plot(0.07:0.18:1.33,sumevents./NEvents.*100,'-^','LineWidth',3,'Color',[143 38 126]./255,'MarkerFaceColor',[177 49 51]./255,'MarkerSize',20)
    h3=legend('Events utilization','box','off')
    h3.Position=[h2(1)-0.05, h2(2)-0.12, 0.7*h2(3), 0.01*h2(4)]
    ylabel('Percentage (%)','FontSize',30)
    end
end
