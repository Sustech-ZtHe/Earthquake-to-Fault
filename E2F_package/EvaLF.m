function EvaLF(LINEAR,Energyratio,NEvents,C,up)
aspect3d=[];U=[];R=[];
sumevents=[];sorx=[];XTicklab={};XTICK=[];YUPmax=[];EDGE=[];
for j=1:length(C)
   numse=find(LINEAR(:,1)==j);
   subflattening=LINEAR(numse,3);
   sumevents=[sumevents;sum(LINEAR(numse,2))];
   subRSM=Energyratio(numse,1);
   [countse,edgese,numberse]=histcounts(subflattening,[0,0.8,1]);
   [counter0,edgeer0,numberer0]=histcounts(subRSM,[0,up,1,1-up+1]);
   countse=flip(countse);
   counter0=flip(counter0);
   yup=max(counter0(1)+counter0(2),countse(2));
   YUPmax=[YUPmax;yup*1.5];
end
for j=1:length(C)
   numse=find(LINEAR(:,1)==j);
   subflattening=LINEAR(numse,3);
   sumevents=[sumevents;sum(LINEAR(numse,2))];
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
   bar(edge+0.2,sum(LINEAR(numse,2))/NEvents.*max(YUPmax),0.1,'FaceColor',[143 38 126]./255,'EdgeColor',[177 49 51]./255,'LineWidth',3)
   XTicklab={XTicklab{:},['C=',num2str(C(j))]};

    aspect3d=[aspect3d;countse(1)];
    R=[R;counter0(1)+counter0(2)];
    U=[U;sum(LINEAR(numse,2))/NEvents];
end
ylim([0 max(YUPmax)]);
h1=legend('3-D Aspect Ratio(0.8~1)','RSM (0.6~1.4)','Utilization rate of events','box','off');
box off
set(gca,'XTick',XTICK,'XTickLabel',XTicklab,'FontSize',20,'YMinorTick','on','LineWidth',2,'TickLength',[0.005 0.05],'YColor',[0 0 0],'XColor',[0 0 0]);
ylabel('High RSM/3-D Aspect Ratio Counts','FontSize',20);
h2=get(h1,'Position');
h1.Position=[h2(1)-0.5, h2(2)+0.02, h2(3), 0.01*h2(4)];
ax1=gca;
pos = get(ax1, 'Position');
axes('position',pos,'color','none','YColor',[143 38 126]./255,'yaxislocation', ...
    'right','YLim',[0 100],'FontSize',20,'XTickLabel','', ...
        'XTick','','YMinorTick','on','LineWidth',2);hold on
ylabel('Utilization rate(%)','FontSize',20);
end
