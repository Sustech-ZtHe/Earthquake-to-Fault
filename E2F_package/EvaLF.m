function EvaLF(ECC,Energyratio,NEvents,C,up,style)
if style==1
    sumevents=[];
    t=tiledlayout(7,1);
    for j=1:8
       numse=find(ECC(:,2)==j);
       subflattening=ECC(numse);
       sumevents=[sumevents;sum(ECC(numse,4))];
       subRSM=Energyratio(numse);
       subER0=subRSM(subRSM<=1);
       subER1=[size(subRSM(subRSM>1),1),0,0,0,0,0,0,0,0,0];
       [countse,edgese,numberse]=histcounts(subflattening,[0:0.1:1]);
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
    sumevents=[];sorx=[];XTicklab={};XTICK=[];YUPmax=[];EDGE=[];
    for j=1:length(C)
       numse=find(ECC(:,2)==j);
       subflattening=ECC(numse,1);
       sumevents=[sumevents;sum(ECC(numse,4))];
       subRSM=Energyratio(numse,1);
       [countse,edgese,numberse]=histcounts(subflattening,[0,0.8,1]);
       [counter0,edgeer0,numberer0]=histcounts(subRSM,[0,up,1,1-up+1]);
       countse=flip(countse);
       counter0=flip(counter0);
       yup=max(counter0(1)+counter0(2),countse(2));
       YUPmax=[YUPmax;yup*1.5];
    end
    for j=1:length(C)
       numse=find(ECC(:,2)==j);
       subflattening=ECC(numse,1);
       sumevents=[sumevents;sum(ECC(numse,4))];
       subRSM=Energyratio(numse,1);
       [countse,edgese,numberse]=histcounts(subflattening,[0,0.8,1]);
       [counter0,edgeer0,numberer0]=histcounts(subRSM,[0,up,1,1-up+1]);
       countse=flip(countse);
       counter0=flip(counter0);
       edge=0+(j-1)*0.5;
       XTICK=[XTICK;edge];
       EDGE=[EDGE;edge+0.2];
       bar(edge,countse(1),0.1,'FaceColor',[203 238 249]./255,'EdgeColor',[66 146 197]./255,'LineWidth',3);hold on;
       bar(edge+0.1,counter0(1)+counter0(2),0.1,'FaceColor',[161 217 156]./255,'EdgeColor',[5 165 158]./255,'LineWidth',3);
       bar(edge+0.2,sum(ECC(numse,4))/NEvents.*max(YUPmax),0.1,'FaceColor',[177 49 51]./255,'EdgeColor',[143 38 126]./255,'LineWidth',3)
       XTicklab={XTicklab{:},['C=',num2str(C(j))]};
    end
    % yTicks = get(gca, 'YTick');
    % yMaxTick = max(yTicks)
    ylim([0 max(YUPmax)]);
    % bar(EDGE,sumevents./NEvents.*max(YUP),0.2,'FaceColor',[177 49 51]./255,'EdgeColor',[143 38 126]./255,'LineWidth',3)
    h1=legend('Linearity (0.8~1)','RSM (0.6~1.4)','Utilization rate of events','box','off');
    box off
    set(gca,'XTick',XTICK,'XTickLabel',XTicklab,'FontSize',20,'YMinorTick','on','LineWidth',2,'TickLength',[0.005 0.05],'YColor',[0 0 0],'XColor',[0 0 0]);
    ylabel('High RSM/Linearity Fault Counts','FontSize',20);
    h2=get(h1,'Position');
    h1.Position=[h2(1)-0.5, h2(2)+0.02, h2(3), 0.01*h2(4)];
    ax1=gca;
    pos = get(ax1, 'Position');
    axes('position',pos,'color','none','YColor',[143 38 126]./255,'yaxislocation', ...
        'right','YLim',[0 100],'FontSize',20,'XTickLabel','', ...
            'XTick','','YMinorTick','on','LineWidth',2);hold on
    % bar(EDGE,sumevents./NEvents.*100,0.2,'FaceColor',[177 49 51]./255,'EdgeColor',[143 38 126]./255,'LineWidth',3)
    % plot(XTICK./100,sumevents./NEvents.*100,'-^','LineWidth',3,'Color',[143 38 126]./255,'MarkerFaceColor',[177 49 51]./255,'MarkerSize',20)
    % h3=legend('Events utilization','box','off');
    % h3.Position=[h2(1)-0.05, h2(2)-0.135, h2(3), 0.01*h2(4)];
    ylabel('Percentage (%)','FontSize',20);
    end
end
