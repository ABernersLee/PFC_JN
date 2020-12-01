function plot_Velocity_Vs_Value(label)
cd('E:\XY_matdata\AllDays')
d2 = dir('*.mat');


replaylab = {'AllArmEvents';'Wc3ArmCov3';'WC4ArmCov5mj7';'WC5ArmCov5mj4'};
ilab = 1; 
grouplabels = {['AllCells_' replaylab{ilab}];['BestTTDays_' replaylab{ilab}];['BestDays_' replaylab{ilab}]};
igroup = 1;
savefolder = ['E:\XY_matdata\Figures\ForPaper\' grouplabels{igroup}];
if ~isfolder(['E:\XY_matdata\Figures\ForPaper\' grouplabels{igroup} '\CompLinearModels\'])    
    mkdir(['E:\XY_matdata\Figures\ForPaper\' grouplabels{igroup} '\CompLinearModels\']) 
end
% bn = [.025;.05;.1;.2;.4;.8];
% bn = [.02;.04;.06;.1;.2;.4;.8];
bn = [.02;.04;.06;.1;.2;.4;.8;1];
% bn = [.2];
% bn = [.04;.1;.2;.8];
% bn = .2;
% bn = [.1;.2];
% bn = NaN;
% ps = NaN(size(d2,1),3,length(bn));
Rs = NaN(size(d2,1),5,length(bn));
CA = []; newR = []; newC = [];
% errorRall = []; errorSall = [];
for ib = 1:length(bn)
    for id = 1:size(d2,1)
        thisdir = d2(id).name;
%         [RsV2,RsT2,RsS2,RsP2,RsTR,coeff,armsig,R1,R2,coeff1,coeff2] = Velocity_Vs_Value(thisdir,bn(ib),label);
        [RsV2,RsT2,RsS2,RsP2,RsTR,RsA,RsProp,pT,pV,pTrev] = Velocity_Vs_Value(thisdir,bn(ib),label,savefolder);
%         ps(id,1,ib) = pT;
%         ps(id,2,ib) = pV;
%         ps(id,3,ib) = pTrev; 
        Rs(id,1,ib) = mean(RsT2);
        Rs(id,7,ib) = mean(RsP2);
        Rs(id,3,ib) = mean(RsV2);
        Rs(id,4,ib) = mean(RsS2);
        Rs(id,2,ib) = mean(RsTR);
        Rs(id,5,ib) = mean(RsA);
        Rs(id,6,ib) = mean(RsProp);
%         CA = cat(1,CA,[coeff armsig]);
%         newR = cat(1,newR,[R1 R2]);
%         newC = cat(1,newC,[coeff1 coeff2]);
    %     errorRall = cat(1,errorRall,errorT);
    %     errorSall = cat(1,errorSall,errorTS);
    end
    disp(['done with ib ' num2str(ib)])
end

if 0 

signrank(newR(:,1),newR(:,2),'tail','right')

figure; hold on; 
subplot(2,1,1), hold on
histogram(abs(CA(CA(:,2)>=.05 & CA(:,3)==1,1)),20,'Normalization','probability'); 
histogram(abs(CA(CA(:,2)<.05  & CA(:,3)==1,1)),20,'Normalization','probability');  
legend('Not','Replay Info Modulated')
p = ranksum(abs(CA(CA(:,2)>=.05,1)),abs(CA(CA(:,2)<.05,1)));
text(10,.18,['ranksum p = ' num2str(round(p,2,'significant'))])
ylabel('Probability (No Low-FR Cells)')
set(gca,'FontSize',18)
subplot(2,1,2), hold on

% legend('Not','Significant')
histogram(abs(CA(CA(:,2)>=.05,1)),20,'Normalization','probability'); 
histogram(abs(CA(CA(:,2)<.05,1)),20,'Normalization','probability');  
p = ranksum(abs(CA(CA(:,2)>=.05 & CA(:,3)==1,1)),abs(CA(CA(:,2)<.05 & CA(:,3)==1,1)));
xlabel('Abs(T-value of Correlation from LM)')
ylabel('Probability (All Cells)')
text(10,.18,['ranksum p = ' num2str(round(p,2,'significant'))])
set(gca,'FontSize',18)
set(gcf,'Position',  [2181          61         888         725])
% helper_savefig(['D:\XY_matdata\Figures\ForPaper\CompLinearModels\Tvalue_TimeFromRewardArmSig_New' label])
helper_saveandclosefig([savefolder '\CompLinearModels\Tvalue_TimeFromRewardArmSig_New' label])

for ihalf = 1:2
    figure; hold on; 
    subplot(2,1,1), hold on
    histogram(abs(newC(CA(:,2)>=.05 & CA(:,3)==1,ihalf)),20,'Normalization','probability'); 
    histogram(abs(newC(CA(:,2)<.05  & CA(:,3)==1,ihalf)),20,'Normalization','probability');  
    legend('Not','Replay Info Modulated')
    p = ranksum(abs(newC(CA(:,2)>=.05,ihalf)),abs(newC(CA(:,2)<.05,ihalf)));
    text(10,.18,['ranksum p = ' num2str(round(p,2,'significant'))])
    ylabel('Probability (No Low-FR Cells)')
    set(gca,'FontSize',18)
    subplot(2,1,2), hold on

    % legend('Not','Significant')
    histogram(abs(newC(CA(:,2)>=.05,ihalf)),20,'Normalization','probability'); 
    histogram(abs(newC(CA(:,2)<.05,ihalf)),20,'Normalization','probability');  
    p = ranksum(abs(newC(CA(:,2)>=.05 & CA(:,3)==1,ihalf)),abs(newC(CA(:,2)<.05 & CA(:,3)==1,ihalf)));
    xlabel('Abs(T-value of Correlation from LM)')
    ylabel('Probability (All Cells)')
    xl = get(gca,'xlim');
    text(xl(1)*1.2,.18,['ranksum p = ' num2str(round(p,2,'significant'))])
    set(gca,'FontSize',18)
    set(gcf,'Position',  [2181          61         888         725])
    helper_saveandclosefig([savefolder '\CompLinearModels\Tvalue_TimeFromRewardArmSig_New' label '_half' num2str(ihalf)])
end


figure; hold on; 
subplot(2,1,1), hold on

legend('Not','Replay Info Modulated')


histogram(CA(CA(:,2)>=.05,1),20,'Normalization','probability'); 
histogram(CA(CA(:,2)<.05,1),20,'Normalization','probability');  

subplot(2,1,2), hold on
histogram(CA(CA(:,2)>=.05 & CA(:,3)==1,1),20,'Normalization','probability'); 
histogram(CA(CA(:,2)<.05  & CA(:,3)==1,1),20,'Normalization','probability');  
% legend('Not','Significant')
end

% clear a
% figure; hold on
% for itype = 1:4
%     for ib = 1:length(bn)
%         m = squeeze(mean(Rs(:,itype,:),1));
%         sem = squeeze(std(Rs(:,itype,:)))./sqrt(size(Rs,1));
%         a(itype) = errorbar(1:size(Rs,3),m,sem,'LineWidth',2);
%     end
%     if itype == 4
%     legend(a,{'Time to Reward';'Velocity'})
%     end
% end
% set(gca,'xtick',1:length(bn),'xticklabel',bn)
% title('Start to Reward, All Laps, No Overfitting')



clear a
figure; hold on
cc = varycolor(size(Rs,2));
for itype = [1:size(Rs,2)]
    for ib = 1:length(bn)
        m = squeeze(mean(Rs(:,itype,:),1));
        sem = squeeze(std(Rs(:,itype,:)))./sqrt(size(Rs,1));
        a(itype) = errorbar(1:size(Rs,3),m,sem,'LineWidth',2,'Color',cc(itype,:));
    end
    if itype == size(Rs,2) 
    legend(a,{'Time to Reward';'Time Since Left Reward';'Velocity';'Speed';'Acceleration';'Proportion of Way to Reward';'Distance To Reward'},'Location','NorthWest')
%     legend(a,{'Time to Reward';'Proprotional Distance';'Velocity';'Speed';'Time Since Left Reward'})
%     legend({'Time to Reward';'Time Since Left Reward'})
    end
end
set(gca,'xtick',1:length(bn),'xticklabel',bn)
tit = 'End to Reward, All Laps, 10-fold crossval, all';
title(tit)
xlabel('Time Window/Bin (Seconds)')
ylabel('R-squared')
set(gca,'FontSize',18)
set(gcf,'Position',[   2329          77         906         696])
helper_saveandclosefig([savefolder '\CompLinearModels\CopLinearModels_newPFC ' tit])
