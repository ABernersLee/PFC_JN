function decode_choice_by_thetaphase_all(dirs,savefolder,toplot,igroup,numshuff)
% can I decode choice by the PFC FR during certain parts of the theta cycle
% yes
indinds = [0 .06];
indinds2 = [-.06 0];
% indinds = [0 .06]; %Nonlocal
% indinds2 = [-.06 0]; %Local, 0 is the end of the theta cycle
% numshuff = 1000; % was 2000
if igroup==3
    toplot = false;
end
% tit = 'Time_Since_Reward_750ms_';
tit = 'Future_Choice_eachdaythetasep_'; % from original version of this, different way of excluding cells (mean<.02 vs sum<5) and so different sessions too

cd(dirs.homedir)
accd = [];accdv2 = [];accdv = []; A = []; AS = []; A2 = []; AS2 = []; L = []; LS = []; L2 = []; LS2 = []; LapsNL = []; LapsL = [];
d2 = dir('*.mat');
% toplot = false;
if ~isfolder([savefolder '\Classify\'])
    mkdir([savefolder '\Classify\'])
end
pid = [];
for id = 1:size(d2,1)
%     if id==4; continue; end %matched sessions with time since reward %id==1 || || id==7
    load(d2(id).name,'other_cells_touse')
    if sum(other_cells_touse(:,igroup))==0 && igroup == 3
        continue
    end
    disp(num2str(id))
    [accdiff1,accdiff2,accdiff1v2,accdiff2v2,accall,accSall,accallL,accSallL,p1,lapaccL,lapaccA] = decode_choice_by_thetaphase(d2(id).name,numshuff,indinds,indinds2,toplot,savefolder);
%     [accdiff1,accdiff2,accall,accSall] = decode_time_since_last_reward_by_thetaphase(d2(id).name,numshuff,indinds,indinds2);    
    accd = [accd;[accdiff1' accdiff2']];
    accdv = [accdv;[accdiff1v2 accdiff2v2]];
    accdv2 = [accdv2;[nanmean(accdiff1v2) nanmean(accdiff2v2)]];
    A = [A;accall];
    AS = [AS;accSall];
    A2 = [A2;nanmean(accall)];
    AS2 = [AS2;nanmean(accSall,1)];
    L = [L;accallL];
    LS = [LS;accSallL];
    L2 = [L2;nanmean(accallL)];
    LS2 = [LS2;nanmean(accSallL,1)];
    LapsNL = [LapsNL;lapaccA];
    LapsL = [LapsL;lapaccL];
    pid = [pid;p1];
    disp(num2str([p1 (sum(nanmean(AS)>=nanmean(A))+1)/(numshuff+1) diff(accdv2(end,:))]))
end


p2 = (sum(nanmean(AS)>=nanmean(A))+1)/(numshuff+1);
figure; histogram(nanmean(AS),.43:.005:.57,'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(A) mean(A)],yl,'r-','LineWidth',3)
title(['All sweep foward, p = ' num2str(p2) ' all'])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' tit 'NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff) '_all'])



p2 = (sum(nanmean(AS)>=nanmean(A))+1)/(numshuff+1);
figure; %histogram(nanmean(AS),.43:.005:.57,'FaceColor','k'); 
hold on; 
[h,c] = hist(nanmean(AS),.43:.005:.57); 
excl = [1:find(h~=0,1,'first')-1 find(h~=0,1,'last'):length(h)];
h(excl) = []; c(excl) = [];
plot(c,h,'.-k','LineWidth',2,'MarkerSize',20)
yl = get(gca,'ylim');
plot([mean(A) mean(A)],yl,'r-','LineWidth',3)
set(gca,'xlim',[.45 .58])
xt = get(gca,'xtick');
set(gca,'xticklabel',xt*100)
title(['All sweep foward, p = ' num2str(p2) ' all'])
set(gcf,'renderer','Painters')
set(gca,'FontSize',18)
ylabel('Shuffle Counts')
xlabel('Accuracy (%)')
helper_saveandclosefig([savefolder '\Classify\' tit 'NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff) '_all_nobars'])

p2 = (sum(nanmean(AS2)>=nanmean(A2))+1)/(numshuff+1);
figure; histogram(nanmean(AS2),.43:.005:.57,'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(A2) nanmean(A2)],yl,'r-','LineWidth',3)
title(['All sweep foward, p = ' num2str(p2) ' means'])
set(gcf,'renderer','Painters')
set(gcf,'Position',[2323         244         816         645])
helper_saveandclosefig([savefolder '\Classify\' tit 'NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff) '_means'])

p2 = (sum(nanmean(LS)>=nanmean(L))+1)/(numshuff+1);
figure; histogram(nanmean(LS),.43:.005:.57,'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(L) mean(L)],yl,'r-','LineWidth',3)
title(['All sweep foward, p = ' num2str(p2) ' all'])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' tit 'Local_' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'ms_numshuff_' num2str(numshuff) '_all'])


p2 = (sum(nanmean(LS)>=nanmean(L))+1)/(numshuff+1);
figure;% histogram(nanmean(LS),.43:.005:.57,'FaceColor','k'); 
hold on; 
[h,c] = hist(nanmean(LS),.43:.005:.57); 
excl = [1:find(h~=0,1,'first')-1 find(h~=0,1,'last'):length(h)];
h(excl) = []; c(excl) = [];
plot(c,h,'.-k','LineWidth',2,'MarkerSize',20)
yl = get(gca,'ylim');
plot([mean(L) mean(L)],yl,'r-','LineWidth',3)
set(gca,'xlim',[.45 .58])
xt = get(gca,'xtick');
set(gca,'xticklabel',xt*100)
title(['All sweep foward, p = ' num2str(p2) ' all'])
set(gcf,'renderer','Painters')
set(gca,'FontSize',18)
ylabel('Shuffle Counts')
xlabel('Accuracy (%)')
helper_saveandclosefig([savefolder '\Classify\' tit 'Local_' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'ms_numshuff_' num2str(numshuff) '_all_nobars'])

p2 = (sum(nanmean(LS2)>=nanmean(L2))+1)/(numshuff+1);
figure; histogram(nanmean(LS2),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(L2) nanmean(L2)],yl,'r-','LineWidth',3)
title(['All sweep foward, p = ' num2str(p2) ' means'])
set(gcf,'renderer','Painters')
set(gcf,'Position',[2323         244         816         645])
helper_saveandclosefig([savefolder '\Classify\' tit 'Local_' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'ms_numshuff_' num2str(numshuff) '_means'])




figure; hold on
subplot(1,2,1); hold on
histogram(LapsL,25,'FaceColor','k'); 
histogram(LapsNL,25,'FaceColor','r'); 
legend('Local','NonLocal')
xlabel('Improvement over Shuffle')
ylabel('Counts')
set(gca,'FontSize',18)

subplot(1,2,2); hold on
% plot(1,accd(:,2),'ok')
plot([.8 4.2],[0 0],'k--')
errorbar(1,nanmean(LapsL),nanstd(LapsL)/sqrt(length(LapsL)),'k','LineWidth',3)
% plot(2,LapsNL,'or')
errorbar(2,nanmean(LapsNL),nanstd(LapsNL)/sqrt(length(LapsNL)),'r','LineWidth',3)
aa = errorbar(3,nanmedian(LapsL),nanstd(LapsL)/sqrt(length(LapsL)),'k','LineWidth',3);
bb=errorbar(4,nanmedian(LapsNL),nanstd(LapsNL)/sqrt(length(LapsNL)),'r','LineWidth',3);
xlim([.8 4.2])
yl = get(gca,'ylim');
% ylim([-.01 yl(2)])

legend([aa bb],'Local','NonLocal','Location','northwest')
ylabel('Improvement over Shuffle')
set(gca,'xtick',[1.5 3.5],'xticklabel',{'Mean';'Median'})
set(gca,'FontSize',18)

set(gcf,'Position',[ 1979         150        1467         706])
p2 = ranksum(accdv(:,1),accdv(:,2),'tail','right');
p3 = signrank(LapsNL,LapsL,'tail','right');
% p3 = (sum(accd(:,2)>=nanmean(accd(:,1)))+1)/(sum(~isnan(accd(:,2)))+1);
suptitle(['Differential PFC decoding of time since reward depending on HP Phase, all p = ' num2str(round(p2,2,'significant')) ',lap means p = ' num2str(round(p3,2,'significant'))])
% suptitle(['Differential PFC decoding of choice depending on HP Phase, p = ' num2str(round(p3,2,'significant'))])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' tit 'ImprovementOverShuffleByPhase_NonLocal' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_Local' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'msnumshuff_' num2str(numshuff)])


% histogram(accd(:,2),50,'FaceColor','k'); 
% histogram(accd(:,1),50,'FaceColor','r'); 
% [h1,c2] = histcounts(accd(:,2));
% [h2,c] = histcounts(accd(:,1)); 

figure; hold on
h1 = histcounts(LapsL,-.6:.1:.6)';
[h2,c] = histcounts(LapsNL,-.6:.1:.6);
b = bar([h1 h2'],'LineWidth',1,'BarWidth',1);
b(2).FaceColor = [0 0 0];
b(1).FaceColor = [1 1 1];
set(gca,'xtick',1.5:5:length(c)-.5,'xticklabel',c(2:5:end)*100)
% set(gca,'xtick',4:9:length(c),'xticklabel',c(4:9:end)*100)
legend('Local','NonLocal')
xlabel('Improvement over Shuffle (%)')
ylabel('Counts')
set(gca,'FontSize',18)
set(gcf,'Position',[2323         244         816         645])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' tit 'ImprovementOverShuffleByPhase2_NonLocal' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_Local' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'msnumshuff_' num2str(numshuff)])


figure; hold on
h1 = hist(LapsL,-.6:.1:.6)';
[h2,c] = hist(LapsNL,-.6:.1:.6);
% b = bar([h1 h2'],'LineWidth',1,'BarWidth',1);
% b(2).FaceColor = [0 0 0];
% b(1).FaceColor = [1 1 1];
plot(c,h1,'.-','color',[.5 .5 .5],'LineWidth',2,'MarkerSize',20)
plot(c,h2,'.-k','LineWidth',2,'MarkerSize',20)
set(gca,'xtick',c(2:5:end),'xticklabel',c(2:5:end)*100)
% set(gca,'xtick',4:9:length(c),'xticklabel',c(4:9:end)*100)
legend('Local','NonLocal')
xlabel('Improvement over Shuffle (%)')
ylabel('Counts')
set(gca,'FontSize',18)
set(gcf,'Position',[2323         244         816         645])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' tit 'ImprovementOverShuffleByPhase2_NonLocal' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_Local' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'msnumshuff_' num2str(numshuff) '_nobars'])

