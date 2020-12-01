
% can I decode choice by the HP FR during certain parts of the theta cycle
% yes
% indinds = [0 .06];
% indinds2 = [-.07 -.01];
indinds = [-.06 0];
indinds2 = [-.12 -.06];
numshuff = 50;

tit = 'Mean.2_WithinCycle';

cd(dirs.homedir)
accd = []; A = []; AS = [];
d2 = dir('*.mat');

for id = 1:size(d2,1)
%     if id ==1 || id==6 || id == 7 || id==4 || id==10; continue; end %these are all with mean(.015)
        %if use .1, can use all sessions
        
%         make_HP_thetaeventtriggeredmat(d2(id).name)
%         disp(['Made for ' num2str(id)])

    disp(id)
    [accdiff1,accdiff2,accall,accSall] = decode_choice_by_thetaphase_hp(d2(id).name,numshuff,indinds,indinds2);
    accd = [accd;[accdiff1' accdiff2']];
    A = [A;accall];
    AS = [AS;accSall];
end

p2 = (sum(nanmean(AS)>=nanmean(A))+1)/(numshuff+1);
figure; histogram(nanmean(AS),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(A) mean(A)],yl,'r-','LineWidth',3)
title(['All sweep foward, p = ' num2str(p2)])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Classify\HP' tit '_NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff)])




figure; hold on
subplot(1,2,1); hold on
histogram(accd(:,2),'FaceColor','k'); 
histogram(accd(:,1),'FaceColor','r'); 
legend('Local','NonLocal')
xlabel('Improvement over Shuffle')
ylabel('Counts')
set(gca,'FontSize',18)

subplot(1,2,2); hold on
plot(1,accd(:,2),'ok')
errorbar(1.2,mean(accd(:,2)),std(accd(:,2))/sqrt(length(accd(:,2))),'k','LineWidth',3)
errorbar(1.3,median(accd(:,2)),std(accd(:,2))/sqrt(length(accd(:,2))),'k','LineWidth',3)
plot(2,accd(:,1),'or')
errorbar(2.2,mean(accd(:,1)),std(accd(:,1))/sqrt(length(accd(:,1))),'r','LineWidth',3)
errorbar(2.3,median(accd(:,1)),std(accd(:,1))/sqrt(length(accd(:,1))),'r','LineWidth',3)
xlim([.7 2.5])
ylabel('Improvement over Shuffle')
set(gca,'xtick',1:2,'xticklabel',{'Local';'NonLocal'})
set(gca,'FontSize',18)

set(gcf,'Position',[ 1979         150        1467         706])
p2 = ranksum(accd(:,1),accd(:,2));
p3 = signrank(accd(:,1),accd(:,2));
suptitle(['Differential HP decoding of choice depending on HP Phase, p = ' num2str(round(p3,2,'significant'))])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Classify\HP' tit  '_ImprovementOverShuffleByPhase_NonLocal' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_Local' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'msnumshuff_' num2str(numshuff)])

