function armcombs_acrosspvals(dirs,cutoff,lab,zs)

if zs == 1
    lab = [lab '_zscore'];
elseif zs == 2
    lab = [lab '_mediansubtract'];
end
label = 'RP';
cd(dirs.homedir)
d2 = dir('*.mat');
pvals = [1:-.025:.05];
r1 = NaN(length(pvals),3,2); p1 = r1; numcells = NaN(length(pvals),1);
pp = p1;
pDiag = NaN(6,length(pvals));
% nS = 1000;
for ip = 1:length(pvals)
Th = []; Fr = []; Re = []; %FrS = [];% Nl = []; Re = []; Plain = [];% allarmsig = []; cellname = []; 
ReShuff = [];ThShuff = []; FrShuff = [];
for id = 1:size(d2,1) %1:size(d2,1)
    load(d2(id).name,'prospectiveFR_armcombs_shuff','prospectiveTheta_armcombs','prospectiveTheta_armcombs_shuff','prospectiveFR_armcombs',...
        [label  '_pfcfwdreplaymodu_armXpos'],[label  '_pfcfwdreplaymodu_armXpos_shuff'],'prospectiveFR_meanFR','thetatrigFR_armcombs') %,'other_cells')
    
    load(d2(id).name, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'])
    eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])
    
    eval(['pfcfwdreplaymodu_armXpos = ' label '_pfcfwdreplaymodu_armXpos;'])
    eval(['pfcfwdreplaymodu_armXpos_shuff = ' label '_pfcfwdreplaymodu_armXpos_shuff;'])

    if ip>1
        touse = armsig(:,1)<pvals(ip);
    else
        touse = armsig(:,3)==1;
    end


%     Fr = [Fr; prospectiveFR_armcombs(touse,:)];
    %new, try out 
    Fr = [Fr; thetatrigFR_armcombs(touse,:)];
    
    
    FrShuff = [FrShuff;prospectiveFR_armcombs_shuff(touse,:,:)];
    
    Re = [Re; pfcfwdreplaymodu_armXpos(touse,:)];
    ReShuff = [ReShuff; pfcfwdreplaymodu_armXpos_shuff(touse,:,:)];
    
    
    Th = [Th; prospectiveTheta_armcombs(touse,:)];    
    ThShuff = [ThShuff; prospectiveTheta_armcombs_shuff(touse,:,:)];
end

if zs == 1
FrZ = nanzscore(Fr,[],2);
ReZ = nanzscore(Re,[],2);
ThZ = nanzscore(Th,[],2);
ThShuffZ = nanzscore(ThShuff,[],2);
ReShuffZ = nanzscore(ReShuff,[],2);
FrShuffZ = nanzscore(FrShuff,[],2); 

elseif zs == 2
    ThZ = Th-nanmedian(Th,2);
ThShuffZ = ThShuff-nanmedian(ThShuff,2);
ReShuffZ = ReShuff-nanmedian(ReShuff,2);
FrShuffZ = FrShuff-nanmedian(FrShuff,2);
FrZ = Fr-nanmedian(Fr,2);
ReZ = Re-nanmedian(Re,2);
end


% FrShuffZ = repmat(FrShuffZ,[1 1 10]); % TAKE OUT WHEN FIX FR SHUFFLEs

corrdat = NaN(6,6,size(ThZ,1),3); 
corrdatS = NaN(6,6,size(ThZ,1),6,max([size(ThShuff,3) size(ReShuff,3)]));
for icell = 1:size(ThZ,1)
    corrdat(:,:,icell,1) = repmat(ReZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    corrdat(:,:,icell,2) = repmat(ReZ(icell,:),[6 1])-(repmat(ThZ(icell,:),[6 1])');
    corrdat(:,:,icell,3) = repmat(ThZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    corrdat(:,:,icell,4) = repmat(ThZ(icell,:),[6 1])-(repmat(ReZ(icell,:),[6 1])');
    corrdat(:,:,icell,5) = repmat(ReZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    corrdat(:,:,icell,6) = repmat(ThZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    for ishuff = 1:size(ReShuffZ,3)
        corrdatS(:,:,icell,1,ishuff) = repmat(ReShuffZ(icell,:,ishuff),[6 1])-(repmat(FrZ(icell,:),[6 1])');
        corrdatS(:,:,icell,2,ishuff) = repmat(ReShuffZ(icell,:,ishuff),[6 1])-(repmat(ThZ(icell,:),[6 1])');
    end
    for ishuff = 1:size(ThShuffZ,3)
        corrdatS(:,:,icell,3,ishuff) = repmat(ThShuffZ(icell,:,ishuff),[6 1])-(repmat(FrZ(icell,:),[6 1])');
        corrdatS(:,:,icell,4,ishuff) = repmat(ThShuffZ(icell,:,ishuff),[6 1])-(repmat(ReZ(icell,:),[6 1])');
    end
    for ishuff = 1:size(FrShuffZ,3)
        corrdatS(:,:,icell,5,ishuff) = repmat(FrShuffZ(icell,:,ishuff),[6 1])-(repmat(ReZ(icell,:),[6 1])');
        corrdatS(:,:,icell,6,ishuff) = repmat(FrShuffZ(icell,:,ishuff),[6 1])-(repmat(ThZ(icell,:),[6 1])');
    end
end


% x3 = [.5; 1.5;1.5; .5];
% x2 = x3*ones(1,6)+ones(4,1)*(0:5);
% y3 = [.5; .5 ; 1.5;1.5];
% y2 = y3*ones(1,6)+ones(4,1)*(0:5);

for idattype = 1:6
    j = squeeze(nanmean(abs(corrdat(:,:,:,idattype)),3));
    real = (nansum(j(eye(6)==1)))./(nansum(nansum(j)));
    fake = NaN(size(ThShuffZ,3),1);
    for in = 1:size(ThShuffZ,3)
        j2 = squeeze(nanmean(abs(corrdatS(:,:,:,idattype,in)),3));
%         j2 = squeeze(nanmean(abs(corrdat(randperm(size(corrdat,1)),:,:,idattype)),3));
        fake(in,1) = (nansum(j2(eye(6)==1)))./(nansum(nansum(j2)));
    end
    pDiag(idattype,ip) = (sum(fake<=real)+1)/(sum(~isnan(squeeze(nanmean(nanmean(nanmean(corrdatS(:,:,:,idattype,:),1),2),3))))+1);
end

% Replay shuffles
[r1(ip,1,1),p1(ip,1,1)] = corr(ReZ(:),ThZ(:),'rows','complete','type','Kendall');
[r1(ip,1,2),p1(ip,1,2)] = corr(ReZ(:),FrZ(:),'rows','complete','type','Kendall');

rshuff = NaN(size(ReShuffZ,3),2);
for ii = 1:size(ReShuffZ,3)
    ReShuffZ2 = ReShuffZ(:,:,ii);
    rshuff(ii,1) = corr(ReShuffZ2(:),ThZ(:),'rows','complete','type','Kendall');
    rshuff(ii,2) = corr(ReShuffZ2(:),FrZ(:),'rows','complete','type','Kendall');
end

pp(ip,1,1) = (sum(rshuff(:,1)>=r1(ip,1,1))+1)/(1+size(ReShuffZ,3));
pp(ip,1,2) = (sum(rshuff(:,2)>=r1(ip,1,2))+1)/(1+size(ReShuffZ,3));

% Theta shuffles
[r1(ip,2,1),p1(ip,2,1)] = corr(ThZ(:),ReZ(:),'rows','complete','type','Kendall');
[r1(ip,2,2),p1(ip,2,2)] = corr(ThZ(:),FrZ(:),'rows','complete','type','Kendall');

rshuff = NaN(size(ThShuffZ,3),2);
for ii = 1:size(ThShuffZ,3)
    ThShuffZ2 = ThShuffZ(:,:,ii);
    rshuff(ii,1) = corr(ThShuffZ2(:),ReZ(:),'rows','complete','type','Kendall');
    rshuff(ii,2) = corr(ThShuffZ2(:),FrZ(:),'rows','complete','type','Kendall');
end

pp(ip,2,1) = (sum(rshuff(:,1)>=r1(ip,2,1))+1)/(1+size(ThShuffZ,3));
pp(ip,2,2) = (sum(rshuff(:,2)>=r1(ip,2,2))+1)/(1+size(ThShuffZ,3));

% Fr shuffles
[r1(ip,3,1),p1(ip,3,1)] = corr(FrZ(:),ReZ(:),'rows','complete','type','Kendall');
[r1(ip,3,2),p1(ip,3,2)] = corr(FrZ(:),ThZ(:),'rows','complete','type','Kendall');

rshuff = NaN(size(FrShuffZ,3),2);
for ii = 1:size(FrShuffZ,3)
    FrShuffZ2 = FrShuffZ(:,:,ii);
    rshuff(ii,1) = corr(FrShuffZ2(:),ReZ(:),'rows','complete','type','Kendall');
    rshuff(ii,2) = corr(FrShuffZ2(:),ThZ(:),'rows','complete','type','Kendall');
end

pp(ip,3,1) = (sum(rshuff(:,1)>=r1(ip,3,1))+1)/(1+size(FrShuffZ,3));
pp(ip,3,2) = (sum(rshuff(:,2)>=r1(ip,3,2))+1)/(1+size(FrShuffZ,3));




numcells(ip) = size(ReZ,1);
end
%%
label = 'RP';
datlab = {'Replays vs Prospective FR','Replays vs Theta','Theta vs Prospective FR','Theta vs Replays','FR vs Replays','FR Vs Theta'};
figure; hold on
for itype = 1:6
    plot(1:length(pvals),pDiag(itype,:),'.-','LineWidth',3,'MarkerSize',20)
end
plot(1:length(pvals),.05*ones(size(pvals)),'k--')
set(gca,'xtick',1:2:length(pvals),'xticklabels',pvals(1:2:end))
legend(datlab)
title('Diagonal across pvals')
set(gcf,'Position',[246         212        1566         730])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Diagonal_acrosspvals_RP_cutoff' num2str(cutoff) '_' lab]) 
%%
axlab = {'Replays';'Theta';'FR'};
for it = 1:3
    oth = setdiff(1:3,it);
    for iy = 1:2
        
    figure; 
    hold on
    yyaxis right
    rr1 = r1(:,it,iy);
    plot(rr1,'k.-','LineWidth',2,'MarkerSize',15)
    aa = plot(find(p1(:,it,iy)<.1),rr1(p1(:,it,iy)<.1),'ko','MarkerSize',20,'LineWidth',3);
    bb = plot(find(pp(:,it,iy)<.1),rr1(pp(:,it,iy)<.1),'k*','MarkerSize',20,'LineWidth',3);
    cc = plot(find(p1(:,it,iy)<.05),rr1(p1(:,it,iy)<.05),'ro','MarkerSize',20,'LineWidth',3);
    dd = plot(find(pp(:,it,iy)<.05),rr1(pp(:,it,iy)<.05),'r*','MarkerSize',20,'LineWidth',3);
    ylabel([axlab{it} ' vs ' axlab{oth(iy)}])

    yyaxis left
    plot(numcells(2:end),'b.','MarkerSize',15)
    ylabel('Number of PFC Cells')
    xlabel('Significance Cutoff for PFC cells')
    set(gca,'xtick',1:2:length(pvals),'xticklabels',pvals(1:2:end))
    set(gca,'FontSize',18)
    % xlim([.7 length(pvals)+.3])
    legend([aa cc bb dd],{'Correlation Trending','Correlation Significant',[axlab{it} ' Shuffle of Correlation Trending'],[axlab{it} ' Shuffle of Correlation Significant']},'Location','north')
    set(gcf,'Position',[246         212        1566         730])
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\' axlab{it} '_vs_' axlab{oth(iy)} '_acrosspvals_RP_cutoff' num2str(cutoff) '_' lab]) 
    end
end







%%
% figure; 
% hold on
% yyaxis right
% plot(r2,'k.-','LineWidth',2,'MarkerSize',15)
% aa = plot(find(p2<.1),r2(p2<.1),'ko','MarkerSize',20,'LineWidth',3);
% bb = plot(find(pp2<.1),r2(pp2<.1),'k*','MarkerSize',20,'LineWidth',3);
% cc = plot(find(p2<.05),r2(p2<.05),'ro','MarkerSize',20,'LineWidth',3);
% dd = plot(find(pp2<.05),r2(pp2<.05),'r*','MarkerSize',20,'LineWidth',3);
% ylabel('Theta - Replay Relationship (r)')
% 
% yyaxis left
% plot(numcells,'b.','MarkerSize',15)
% ylabel('Number of PFC Cells')
% xlabel('Significance Cutoff for PFC cells')
% set(gca,'xtick',1:2:length(pvals),'xticklabels',pvals(1:2:end))
% set(gca,'FontSize',18)
% xlim([.7 length(pvals)+.3])
% legend([aa cc bb dd],{'Correlation Trending','Correlation Significant','Shuffle of Correlation Trending','Shuffle of Correlation Significant'},'Location','north')
% set(gcf,'Position',[246         212        1566         730])
% helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\ForwardReplaysVsThetaFR_acrosspvals_RP_mod.2_nobaseline_allevents_CandEventSigToo_cutoff' num2str(cutoff) '_' lab])  