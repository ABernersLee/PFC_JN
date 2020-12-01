function plot_armcombs_relationship(dirs,label,sigmodu,addlabel,zs,igroup,savefolder)
if ~isfolder([savefolder '\ArmCombinations\'])        
    mkdir([savefolder '\ArmCombinations\'])        
end

if zs == 1
    addlabel = [addlabel '_zscore'];
elseif zs == 2
    addlabel = [addlabel '_norm'];
end
%%
Nl = []; Th = []; Fr = []; Re = []; Plain = []; La = [];% allarmsig = []; 
cellname = []; Th_same = []; Th_diff = [];
ThShuff = []; ReShuff = [];
cd(dirs.homedir)
d2 = dir('*.mat');
SigP = [];
for id = 1:size(d2,1)
    
%     make_prospective_theta_triggered(d2(id).name,0); %.34);
    load(d2(id).name,'prospectiveFR_armcombslast','prospectiveTheta_armcombs2','prospectiveTheta_armcombs','prospectiveTheta_armcombs_shuffx','prospectiveFR_armcombs',...
        [label  '_pfcfwdreplaymodu_armXpos'],[label  '_pfcfwdreplaymodu_armXpos_shuff'],'prospectiveFR_meanFR',...
        'prospectiveTheta_armcombs_same','prospectiveTheta_armcombs_diff','thetatrigFR_armcombs') %,'other_cells')
    
    
    load(d2(id).name, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'],'other_cells_touse','other_cells')
    eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])
    
    eval(['pfcfwdreplaymodu_armXpos = ' label '_pfcfwdreplaymodu_armXpos;'])
    eval(['pfcfwdreplaymodu_armXpos_shuff = ' label '_pfcfwdreplaymodu_armXpos_shuff;'])

    if sigmodu>0
        if sigmodu==3        
            touse = armsig(:,1)<.05;
            sigmodulabel = ['_ArmReplaySigModCellsOnly' addlabel];
        elseif sigmodu==2        
            touse = armsig(:,1)<.1;
            sigmodulabel = ['_ArmReplayTrendingModCellsOnly' addlabel];
        elseif sigmodu==1
            sigmodulabel = ['_CandEventSigModCellsOnly' addlabel];
            touse = armsig(:,3)==1;
        end
    else
        sigmodulabel = addlabel;
        touse = true(size(armsig,1),1);
    end
%     if id ==3
%         touse(16,1) = false;
%     end
    touse = touse & other_cells_touse(:,igroup);
    cellname = cat(1,cellname,[other_cells(touse) ones(sum(touse),1)*id]);
    SigP = [SigP; armsig(:,1:2)];
    
    Th = [Th; prospectiveTheta_armcombs2(touse,:)];
    Th_same = [Th_same; prospectiveTheta_armcombs_same(touse,:)];
    Th_diff = [Th_diff; prospectiveTheta_armcombs_diff(touse,:)];
%         Fr = [Fr; prospectiveFR_armcombs(touse,:)];
    %new, try out 
    Fr = [Fr; thetatrigFR_armcombs(touse,:)];
    Re = [Re; pfcfwdreplaymodu_armXpos(touse,:)];
    ReShuff = [ReShuff; pfcfwdreplaymodu_armXpos_shuff(touse,:,:)];
    La = [La; prospectiveFR_armcombslast(touse,:)];
    ThShuff = [ThShuff; prospectiveTheta_armcombs_shuffx(touse,:,:)];
    Nl = [Nl;nanmean(cat(3,prospectiveFR_armcombslast(touse,:),prospectiveFR_armcombs(touse,:)),3)];
    Plain = [Plain;prospectiveFR_meanFR(touse,:)];                
%     else
%         Th = [Th; prospectiveTheta_armcombs];
% %         Fr = [Fr; prospectiveFR_armcombs];        
%         %new, try out 
%         Fr = [Fr; thetatrigFR_armcombs];
%         Re = [Re; pfcfwdreplaymodu_armXpos];    
%         La = [La; prospectiveFR_armcombslast];
%         Nl = [Nl;nanmean(cat(3,prospectiveFR_armcombslast,prospectiveFR_armcombs),3)];
%         Plain = [Plain;prospectiveFR_meanFR];
%         ThShuff = [ThShuff; prospectiveTheta_armcombs_shuff];
%         ReShuff = [ReShuff; pfcfwdreplaymodu_armXpos_shuff];

%     end
%     allarmsig = [allarmsig;armsig];
%     cellname = [cellname; other_cells repmat(str2num(d2(id).name(3:end-4)),[length(other_cells) 1])];
    disp(num2str(id))
    if zs == 1
    ThZ = nanzscore(Th,[],2);
    ReZ = nanzscore(Re,[],2);
    elseif zs == 2
%     ThZ = Th-nanmedian(Th,2);
%     ReZ = Re-nanmedian(Re,2);
    ThZ = (Th-min(Th,[],2))./(range(Th,2));
    ReZ = (Re-min(Re,[],2))./(range(Re,2));
    end
    if ~isempty(ThZ)
        [r1,p1] = corr(ThZ(:),ReZ(:),'rows','complete','type','Spearman')
    end
end

Re2 = [nanmean([Re(:,2) Re(:,3)],2) nanmean([Re(:,1) Re(:,6)],2) nanmean([Re(:,4) Re(:,5)],2)];
Fr2 = [nanmean([Fr(:,2) Fr(:,3)],2) nanmean([Fr(:,1) Fr(:,6)],2) nanmean([Fr(:,4) Fr(:,5)],2)];
Th2 = [nanmean([Th(:,2) Th(:,3)],2) nanmean([Th(:,1) Th(:,6)],2) nanmean([Th(:,4) Th(:,5)],2)];

if zs == 1
ThZ = nanzscore(Th,[],2);
ThZ_same = nanzscore(Th_same,[],2);
ThZ_diff = nanzscore(Th_diff,[],2);
ReZ = nanzscore(Re,[],2);
ThShuffZ = nanzscore(ThShuff,[],2);
ReShuffZ = nanzscore(ReShuff,[],2);
NlZ = nanzscore(Nl,[],2);
FrZ = nanzscore(Fr,[],2);

Re2Z = nanzscore(Re2,[],2);
Fr2Z = nanzscore(Fr2,[],2);
Th2Z = nanzscore(Th2,[],2);
PlainZ = nanzscore(Plain,[],2);
LaZ = nanzscore(La,[],2);
elseif zs == 2
% ThZ = Th-nanmedian(Th,2);
% ThShuffZ = ThShuff-nanmedian(ThShuff,2);
% ReShuffZ = ReShuff-nanmedian(ReShuff,2);
% NlZ = Nl-nanmedian(Nl,2);
% FrZ = Fr-nanmedian(Fr,2);
% ReZ = Re-nanmedian(Re,2);
% Re2Z = Re2-nanmedian(Re2,2);
% Fr2Z = Fr2-nanmedian(Fr2,2);
% Th2Z = Th2-nanmedian(Th2,2);
% PlainZ = Plain-nanmedian(Plain,2);
% LaZ = La-nanmedian(La,2);

ThZ = (Th-min(Th,[],2))./(range(Th,2));
ThZ_same = (Th_same-min(Th_same,[],2))./(range(Th_same,2));
ThZ_diff = (Th_diff-min(Th_diff,[],2))./(range(Th_diff,2));
ReZ = (Re-min(Re,[],2))./(range(Re,2));
Th2Z = (Th2-min(Th2,[],2))./(range(Th2,2));
Re2Z = (Re2-min(Re2,[],2))./(range(Re2,2));
ThShuffZ = (ThShuff-min(ThShuff,[],2))./(range(ThShuff,2));
ReShuffZ = (ReShuff-min(ReShuff,[],2))./(range(ReShuff,2));
NlZ = (Nl-min(Nl,[],2))./(range(Nl,2));
FrZ = (Fr-min(Fr,[],2))./(range(Fr,2));
Fr2Z = (Fr2-min(Fr2,[],2))./(range(Fr2,2));
PlainZ = (Plain-min(Plain,[],2))./(range(Plain,2));
LaZ = (La-min(La,[],2))./(range(La,2));
end
%%
rvalues = NaN(size(Re,1),1);
rvaluesS = NaN(size(Re,1),size(ReShuff,3));
% ind = sum(isnan(Re),2)==0;
for icell = 1:size(Re,1)  
%     if ~ind(icell)
%         continue
%     end
    rvalues(icell,1) = corr(Re(icell,:)',Th(icell,:)','rows','complete','type','Spearman'); 
    rvaluesS(icell,:) = corr(squeeze(ReShuffZ(icell,:,:)),Th(icell,:)','rows','complete','type','Spearman'); 
end
p = signrank(rvalues,0,'tail','right')
p2 = (sum(nanmean(rvaluesS)>=nanmean(rvalues))+1)./(size(rvaluesS,2)+1)
 %%
 [h,c] = histcounts(rvalues,[-1:.25:1]);     
 figure; hold on; 
 b = bar(c(1:end-1)+.125,h,'LineWidth',3,'FaceColor','w');
%  set(gca,'xtick',.5:length(h)+.5,'xticklabel',c)
 b(1).EdgeColor = 'k';
 xlabel('Spearman''s rho')
 ylabel('PFC cells')
 yl = get(gca,'ylim');
 yl(2) = 8;
 plot([0 0],yl,'k--')
 plot([nanmedian(rvalues) nanmedian(rvalues)],yl,'r-','LineWidth',3)
 title(['onesided sign-rank p = ' num2str(p) ', permutation test ' num2str(p2)])
 helper_savefig([savefolder '\ArmCombinations\Figure3o_ReplayVSTheta_Bars' sigmodulabel '_Spearman'])
 helper_saveandclosefig([savefolder '\Figure3\Figure3o_ReplayVSTheta_Bars' sigmodulabel '_Spearman'])
 %%
 figure; hold on
 c = cdfplot(rvalues);
 c.LineWidth = 3;
 c.Color = 'r';
 d = cdfplot(rvaluesS(:));
 d.LineWidth = 3;
 d.Color = 'k';
 xlabel('Spearman''s rho')
 helper_savefig([savefolder '\ArmCombinations\Figure3o_ReplayVSTheta_cdf' sigmodulabel '_Spearman'])
 helper_saveandclosefig([savefolder '\Figure3\Figure3o_ReplayVSTheta_cdf' sigmodulabel '_Spearman'])

 %%
 
[r1,p1] = corr(ThZ(:),ReZ(:),'rows','complete','type','Spearman')
[r1_same,p1_same] = corr(ThZ_same(:),ReZ(:),'rows','complete','type','Spearman')
[r1_diff,p1_diff] = corr(ThZ_diff(:),ReZ(:),'rows','complete','type','Spearman')
 %%
 
 
rvalues = NaN(size(Re,1),1);
rvaluesS = NaN(size(Re,1),size(ReShuff,3));
% ind = sum(isnan(Re),2)==0;
for icell = 1:size(Re,1)  
%     if ~ind(icell)
%         continue
%     end
    rvalues(icell,1) = corr(Re(icell,:)',Fr(icell,:)','rows','complete','type','Pearson'); 
    rvaluesS(icell,:) = corr(squeeze(ReShuff(icell,:,:)),Fr(icell,:)','rows','complete','type','Spearman'); 
end
p = signrank(rvalues,0,'tail','right')
p2 = (sum(nanmean(rvaluesS)>=nanmean(rvalues))+1)./(size(rvaluesS,2)+1)
 %%
 [h,c] = histcounts(rvalues,[-1:.25:1]);     
 figure; hold on; 
 b = bar(c(1:end-1)+.125,h,'LineWidth',3,'FaceColor','w');
%  set(gca,'xtick',.5:length(h)+.5,'xticklabel',c)
 b(1).EdgeColor = 'k';
 xlabel('Spearman''s rho')
 ylabel('PFC cells')
 yl = get(gca,'ylim');
 yl(2) = 8;
 plot([0 0],yl,'k--')
 plot([nanmedian(rvalues) nanmedian(rvalues)],yl,'r-','LineWidth',3)
 title(['onesided sign-rank p = ' num2str(p) ', permutation test ' num2str(p2)])
 helper_savefig([savefolder '\ArmCombinations\Figure3o_ReplayVSTheta_Bars' sigmodulabel '_Spearman_100ms'])
 helper_saveandclosefig([savefolder '\Figure3\Figure3o_ReplayVSTheta_Bars' sigmodulabel '_Spearman_100ms'])
 %%
 figure; hold on
 c = cdfplot(rvalues);
 c.LineWidth = 3;
 c.Color = 'r';
 d = cdfplot(rvaluesS(:));
 d.LineWidth = 3;
 d.Color = 'k';
 xlabel('Spearman''s rho')
helper_savefig([savefolder '\ArmCombinations\Figure3o_ReplayVSTheta_cdf' sigmodulabel '_Spearman_100ms'])
 helper_saveandclosefig([savefolder '\Figure3\Figure3o_ReplayVSTheta_cdf' sigmodulabel '_Spearman_100ms'])


%%
% rshuff = NaN(size(ReShuffZ,3),1);
% for ii = 1:size(ReShuffZ,3)
%     ReShuffZ2 = ReShuffZ(:,:,ii);
%     rshuff(ii) = corr(ThZ(:),ReShuffZ2(:),'rows','complete','type','Spearman'); % 'Pearson'); % 'Spearman');
% end
% rshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1)
% rshuff = NaN(size(ThShuffZ,3),1);
% for ii = 1:size(ThShuffZ,3)
%     ReShuffZ2 = ThShuffZ(:,:,ii);
%     rshuff(ii) = corr(ReZ(:),ReShuffZ2(:),'rows','complete','type','Spearman'); % 'Pearson'); % 'Spearman');
% end
% tshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1)
% % rshuff = NaN(size(ThShuffZ,3),size(ReShuffZ,3));
% % for ii = 1:size(ThShuffZ,3)
% %     ThShuffZ2 = ThShuffZ(:,:,ii);
% %     for iii = 1:size(ReShuffZ,3)
% %         ReShuffZ2 = ReShuffZ(:,:,ii);
% %         rshuff(ii,iii) = corr(ReShuffZ2(:),ThShuffZ2(:),'rows','complete','type','Spearman');
% %     end
% % end
% % rtshuffp = (sum(rshuff(:)>=r1)+1)/(length(rshuff(:))+1)

[r1,p1] = corr(ThZ(:),ReZ(:),'rows','complete','type','Spearman')
[r1_same,p1_same] = corr(ThZ_same(:),ReZ(:),'rows','complete','type','Spearman')
[r1_diff,p1_diff] = corr(ThZ_diff(:),ReZ(:),'rows','complete','type','Spearman')
% rshuff = NaN(size(ReShuffZ,3),1);
% for ii = 1:size(ReShuffZ,3)
%     ReShuffZ2 = ReShuffZ(:,:,ii);
%     rshuff(ii) = corr(ThZ(:),ReShuffZ2(:),'rows','complete'); % 'Spearman'); % 'Spearman');
% end
% rshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1)
% rshuff = NaN(size(ThShuffZ,3),1);
% for ii = 1:size(ThShuffZ,3)
%     ReShuffZ2 = ThShuffZ(:,:,ii);
%     rshuff(ii) = corr(ReZ(:),ReShuffZ2(:),'rows','complete'); % 'Pearson'); % 'Spearman');
% end
% tshuffp = (sum(rshuff>=r1)+1)/(size(ThShuffZ,3)+1)
% % rshuff = NaN(size(ThShuffZ,3),size(ReShuffZ,3));
% % for ii = 1:size(ThShuffZ,3)
% %     ThShuffZ2 = ThShuffZ(:,:,ii);
% %     for iii = 1:size(ReShuffZ,3)
% %         ReShuffZ2 = ReShuffZ(:,:,ii);
% %         rshuff(ii,iii) = corr(ReShuffZ2(:),ThShuffZ2(:),'rows','complete'); % 'Pearson'); % 'Spearman');
% %     end
% % end
% % rtshuffp = (sum(rshuff(:)>=r1)+1)/(length(rshuff(:))+1)

% [r1,p1] = corr(Th2Z(:),Re2Z(:),'rows','complete','type','Spearman')
% [r1,p1] = corr(Th2Z(:),Re2Z(:),'rows','complete','type','Pearson')
% [r1,p1] = corr(FrZ(:),ReZ(:),'rows','complete','type','Spearman')
ThZS = [ThZ(:,3) ThZ(:,2) ThZ(:,6) ThZ(:,1) ThZ(:,5) ThZ(:,4)];
% nS = 5000;
%%

% [mx,mm] = nanmax(Th2,[],2);
% mn = nanmin(Th2,[],2);
% [~,si] = sort(mm);
% mx2 = nanmax(Re2,[],2);
% mn2 = nanmin(Re2,[],2);
% 
% figure; hold on
% subplot(1,2,1); hold on
% imagesc((Th2(si,:)-mn(si))./(mx(si)-mn(si)))
% colormap gray
% subplot(1,2,2); hold on
% imagesc((Re2(si,:)-mn2(si))./(mx2(si)-mn2(si)))
% colormap gray
% %%
% 
% [mx,mm] = max(Re2,[],2);
% mn = min(Re2,[],2);
% [~,si] = sort(mm);
% mx2 = max(Th2,[],2);
% mn2 = min(Th2,[],2);
% 
% figure; hold on
% subplot(1,2,1); hold on
% imagesc((Re2(si,:)-mn(si))./(mx(si)-mn(si)))
% subplot(1,2,2); hold on
% imagesc((Th2(si,:)-mn2(si))./(mx2(si)-mn2(si)))

%%
Th2Zord = NaN(size(Re2,1),3);
for icell = 1:size(Re2,1)
    [~,cellorder] = sort(Re2(icell,:));
    Th2Zord(icell,:) = Th2(icell,cellorder);
end

figure; hold on; bar(1:3,nanmedian(Th2Zord),'FaceColor','w','LineWidth',2)

% figure; hold on; 
for iarm = 1:3
    errorbar(iarm,nanmedian(Th2Zord(:,iarm)),nanstd(Th2Zord(:,iarm))./sqrt(sum(~isnan(Th2Zord(:,iarm)))),'k','LineWidth',2)
end
xlabel('Arm (in rank order of replay modulation)')
set(gca,'xtick',1:3,'xticklabels',{'Lowest';'Middle';'Highest'})
ylabel('Theta modulation')
set(gca,'FontSize',18','FontName','Ariel')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVsReplay_newfig' sigmodulabel])
%%
% figure; hold on
% boxplot(Th2Zord)

% ThetaVsReplay_newfig
%
Th2Zord = NaN(size(Re2Z,1),3);
for icell = 1:size(Re2Z,1)
    [~,cellorder] = sort(Re2Z(icell,:));
    Th2Zord(icell,:) = Th2Z(icell,cellorder);
end

figure; hold on; bar(1:3,nanmean(Th2Zord),'FaceColor','w','LineWidth',2)

% figure; hold on; 
for iarm = 1:3
    errorbar(iarm,nanmean(Th2Zord(:,iarm)),nanstd(Th2Zord(:,iarm))./sqrt(sum(~isnan(Th2Zord(:,iarm)))),'k','LineWidth',2)
end
% plot(Th2Zord','k')
xlabel('Arm (in rank order of replay modulation)')
set(gca,'xtick',1:3,'xticklabels',{'Lowest';'Middle';'Highest'})
ylabel('Z-scored theta modulation')
set(gca,'FontSize',18','FontName','Ariel')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVsReplay_newfigZ' sigmodulabel])
% figure; hold on
% boxplot(Th2Zord)
%%
nn = num2str(size(Th,1));
[r,p] = corr(Th(:),Fr(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ThZ(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Th,Fr,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by HP Theta Representation')
ylabel('Modulation by Prospective Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
rshuff = NaN(size(ThShuffZ,3),1);
for ii = 1:size(ThShuffZ,3)
    ThShuffZ2 = ThShuffZ(:,:,ii);
    rshuff(ii) = corr(ThShuffZ2(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
end
rshuffp = (sum(rshuff>=r1)+1)/(size(ThShuffZ,3)+1);
t1 = text(xl(2)*.5,yl(2)*.7,['Shuffle permutation test p = ' num2str(round(rshuffp,2,'significant'))]);
if rshuffp<.05
    t1.Color = 'r';
end
title(['Theta vs Prospective FR n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVSProspectiveFR' sigmodulabel])



%%
[r,p] = corr(Th(:),Re(:),'rows','complete','type','Spearman'); % 'Spearman');
% [r1,p1] = corr(ThZ(:),ReZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Th,Re,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(Th(~isnan(Th) & ~isnan(Re)),Re(~isnan(Th) & ~isnan(Re)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)
xlabel('Modulation by HP Theta Representation')
ylabel('Modulation by HP Replays')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
% t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
% if p1<.05
%     t1.Color = 'r';
% end
% rshuff = NaN(size(ThShuffZ,3),1);
% for ii = 1:size(ThShuffZ,3)
%     ThShuffZ2 = ThShuffZ(:,:,ii);
%     rshuff(ii) = corr(ThShuffZ2(:),ReZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% end
% rshuffp = (sum(rshuff>=r1)+1)/(size(ThShuffZ,3)+1);
% t1 = text(xl(2)*.5,yl(2)*.7,['Theta shuffle perm test p = ' num2str(round(rshuffp,2,'significant'))]);
% if rshuffp<.05
%     t1.Color = 'r';
% end
% rshuff = NaN(size(ReShuffZ,3),1);
% for ii = 1:size(ReShuffZ,3)
%     ReShuffZ2 = ReShuffZ(:,:,ii);
%     rshuff(ii) = corr(ThZ(:),ReShuffZ2(:),'rows','complete','type','Spearman'); % 'Spearman');
% end
% rshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1);
% t1 = text(xl(2)*.5,yl(2)*.6,['Replay shuffle perm test p = ' num2str(round(rshuffp,2,'significant'))]);
% if rshuffp<.05
%     t1.Color = 'r';
% end
title(['Theta vs Replays ' label ' n=' nn])
set(gcf,'Position',[680   431   928   547])
set(gcf,'renderer','Painters')
%%
helper_savefig([savefolder '\ArmCombinations\ThetaVSReplays_' label sigmodulabel])
helper_saveandclosefig([savefolder '\Figure3\ThetaVSReplays_' label sigmodulabel])
%%
% [r,p] = corr(Th(:),Re(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ThZ(:),ReZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(ThZ,ReZ,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(ThZ(~isnan(ThZ) & ~isnan(ReZ)),ReZ(~isnan(ThZ) & ~isnan(ReZ)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)
xlabel('Modulation by HP Theta Representation')
ylabel('Modulation by HP Replays')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
% t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
% if p<.05
%     t1.Color = 'r';
% end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
% rshuff = NaN(size(ThShuffZ,3),1);
% for ii = 1:size(ThShuffZ,3)
%     ThShuffZ2 = ThShuffZ(:,:,ii);
%     rshuff(ii) = corr(ThShuffZ2(:),ReZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% end
% rshuffp = (sum(rshuff>=r1)+1)/(size(ThShuffZ,3)+1);
% t1 = text(xl(2)*.5,yl(2)*.7,['Theta shuffle perm test p = ' num2str(round(rshuffp,2,'significant'))]);
% if rshuffp<.05
%     t1.Color = 'r';
% end
rshuff = NaN(size(ReShuffZ,3),1);
for ii = 1:size(ReShuffZ,3)
    ReShuffZ2 = ReShuffZ(:,:,ii);
    rshuff(ii) = corr(ThZ(:),ReShuffZ2(:),'rows','complete','type','Spearman'); % 'Pearson'); % 'Spearman');
end
rshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1);
t1 = text(xl(2)*.5,yl(2)*.6,['Replay shuffle perm test p = ' num2str(round(rshuffp,2,'significant'))]);
if rshuffp<.05
    t1.Color = 'r';
end
title(['Theta vs Replays ' label ' n=' nn])
set(gcf,'Position',[680   431   928   547])
set(gcf,'renderer','Painters')
helper_savefig([savefolder '\ArmCombinations\ThetaVSReplaysZ_' label sigmodulabel])
helper_saveandclosefig([savefolder '\Figure3\ThetaVSReplaysZ_' label sigmodulabel])
%%
[r,p] = corr(Re(:),Fr(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ReZ(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Re,Fr,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Prospective Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
rshuff = NaN(size(ReShuffZ,3),1);
for ii = 1:size(ReShuffZ,3)
    ReShuffZ2 = ReShuffZ(:,:,ii);
    rshuff(ii) = corr(ReShuffZ2(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
end
rshuffp = (sum(rshuff>=r1)+1)/(size(ReShuffZ,3)+1);
t1 = text(xl(2)*.5,yl(2)*.7,['Replay shuffle perm test p = ' num2str(round(rshuffp,2,'significant'))]);
if rshuffp<.05
    t1.Color = 'r';
end
title([label ' Forward Replays vs Prospective FR n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ForwardReplaysVsProspectiveFR_' label sigmodulabel])
%% retrospective
[r,p] = corr(Re(:),La(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ReZ(:),LaZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Re,La,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Retrospective Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Forward Replays vs Retrospective FR'])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ForwardReplaysVsRetrospectiveFR_' label sigmodulabel])

[r,p] = corr(Th(:),La(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ThZ(:),LaZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Th,La,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by Theta')
ylabel('Modulation by Retrospective Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Theta vs Retrospective FR  n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVsRetrospectiveFR_' label sigmodulabel])
%% nonlocal

[r,p] = corr(Re(:),Nl(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ReZ(:),NlZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Re,Nl,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Non-local Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Forward Replays vs Non-local FR n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ForwardReplaysVs_NonLocalFR_' label sigmodulabel])

[r,p] = corr(Th(:),Nl(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(ThZ(:),NlZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Th,Nl,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by Theta')
ylabel('Modulation by Non-local Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Theta vs Non-local FR n=' nn])
set(gcf,'Position',[680   431   928   547])

set(gcf,'Position',[680   431   928   547])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVs_NonLocalFR_' label sigmodulabel])
%%
% [r1,p1] = corr(ThZ(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% [r1,p1] = corr(ThZS(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% [r,p] = partialcorr(ThZ(:),FrZ(:),ThZS(:),'rows','complete');
% [r,p] = partialcorr(ThZS(:),FrZ(:),ThZ(:),'rows','complete');
% 
% [r1,p1] = corr(ReZ(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% numshuff = 5000;
% rs = NaN(numshuff,1);
% for ii = 1:numshuff 
%     R2 = ReZ(randperm(size(ReZ,1)),:); 
%     rs(ii) = corr(R2(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman'); 
% end
% (sum(rs>=r1)+1)/(numshuff+1)
% 
% [r1,p1] = corr(ReZ(:),ThZ(:),'rows','complete','type','Spearman'); % 'Spearman');
% numshuff = 5000;
% rs = NaN(numshuff,1);
% for ii = 1:numshuff 
%     R2 = ReZ(randperm(size(ThZ,1)),:); 
%     rs(ii) = corr(R2(:),FrZ(:),'rows','complete','type','Spearman'); % 'Spearman'); 
% end
% (sum(rs>=r1)+1)/(numshuff+1)
%% Plain FR

[r,p] = corr(Re2(:),Plain(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
% plot(Re2,Plain,'ok','MarkerSize',10,'LineWidth',3)
plot(Re2,Plain,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(Re2(~isnan(Re2) & ~isnan(Plain)),Plain(~isnan(Re2) & ~isnan(Plain)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)

xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Plain Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
% t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
% if p1<.05
%     t1.Color = 'r';
% end
title([label ' Replays vs Plain FR n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ReplaysVsPlainFR_' label sigmodulabel])


[r1,p1] = corr(Re2Z(:),PlainZ(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
% plot(Re2Z,PlainZ,'ok','MarkerSize',10,'LineWidth',3)
plot(Re2Z,PlainZ,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(Re2Z(~isnan(Re2Z) & ~isnan(PlainZ)),PlainZ(~isnan(Re2Z) & ~isnan(PlainZ)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)
xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Plain Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
% t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
% if p<.05
%     t1.Color = 'r';
% end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Replays vs Plain FR n=' nn])
set(gcf,'Position',[680   431   928   547])

set(gcf,'Position',[680   431   928   547])
set(gcf,'renderer','Painters')
helper_savefig([savefolder '\ArmCombinations\ReplaysVsPlainFRZ_' label sigmodulabel])
helper_saveandclosefig([savefolder '\Figure3\ReplaysVsPlainFRZ_' label sigmodulabel])

%%


[r,p] = corr(Re2(:),Th2(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
% plot(Re2,Th2,'ok','MarkerSize',10,'LineWidth',3)
plot(Re2,Th2,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(Re2(~isnan(Re2) & ~isnan(Th2)),Th2(~isnan(Re2) & ~isnan(Th2)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)

xlabel('Modulation by Replays')
ylabel('Modulation by Theta Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
% t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
% if p1<.05
%     t1.Color = 'r';
% end
title([label ' Replays vs Theta (collapsed) n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ReplaysVsTheta_collapsed_' label sigmodulabel])


[r1,p1] = corr(Re2Z(:),Th2Z(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
% plot(Re2Z,Th2Z,'ok','MarkerSize',10,'LineWidth',3)
plot(Re2Z,Th2Z,'ok','MarkerSize',10,'LineWidth',3)
lm = polyfit(Re2Z(~isnan(Re2Z) & ~isnan(Th2Z)),Th2Z(~isnan(Re2Z) & ~isnan(Th2Z)),1);
xl = get(gca,'xlim');
y2 = polyval(lm,xl);
plot(xl,y2,'r','LineWidth',3)
xlabel('Modulation by Replays')
ylabel('Modulation by Theta')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
% t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
% if p<.05
%     t1.Color = 'r';
% end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Replays vs Theta, collapsed, n=' nn])
set(gcf,'Position',[680   431   928   547])

set(gcf,'Position',[680   431   928   547])
set(gcf,'renderer','Painters')
helper_savefig([savefolder '\ArmCombinations\ReplaysVsTheta_collapsed_' label sigmodulabel])
helper_saveandclosefig([savefolder '\Figure3\ReplaysVsTheta_collapsed_' label sigmodulabel])
%%

[r,p] = corr(Re2(:),Fr2(:),'rows','complete','type','Spearman'); % 'Spearman');
[r1,p1] = corr(Re2Z(:),Fr2Z(:),'rows','complete','type','Spearman'); % 'Spearman');
figure; hold on
plot(Re2,Fr2,'ok','MarkerSize',10,'LineWidth',3)
xlabel('Modulation by HP Foward Replays')
ylabel('Modulation by Future Arm')
set(gca,'FontSize',18)
xl = get(gca,'xlim'); yl = get(gca,'ylim');
t1 = text(xl(2)*.5,yl(2)*.9,['Across Cell corr: r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    t1.Color = 'r';
end
t1 = text(xl(2)*.5,yl(2)*.8,['Within Cell corr: r = ' num2str(round(r1,2,'significant')) ' p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end
title([label ' Replays vs Future FR n=' nn])
set(gcf,'Position',[680   431   928   547])
helper_saveandclosefig([savefolder '\ArmCombinations\ReplaysVsFutureFR_' label sigmodulabel])

% %%
% [~,a] = max(Re,[],2);
% [~,b] = max(Fr,[],2);
% j = confusionmat(a,b);
% real = (sum(j(eye(6)==1)))./(sum(sum(j)));
% 
% fake = NaN(nS,1);
% for in = 1:nS
%     j = confusionmat(a(randperm(length(a))),b);
%     fake(in,1) = (sum(j(eye(6)==1)))./(sum(sum(j)));
% end
% p = (sum(fake>=real)+1)/(nS+1);
% figure; histogram(fake,'FaceColor','k')
% yl = get(gca,'ylim'); hold on;
% plot([real real],yl,'r','LineWidth',3)
% title([label ' Forward Replays vs Prospective FR, p = ' num2str(p)])
% 
% helper_saveandclosefig([savefolder '\ArmCombinations\ForwardReplaysVsProspectiveFR_' label '_max' sigmodulabel])
% 
% [~,a] = max(Re,[],2);
% [~,b] = max(Th,[],2);
% j = confusionmat(a,b);
% real = (sum(j(eye(6)==1)))./(sum(sum(j)));
% 
% fake = NaN(nS,1);
% for in = 1:nS
%     j = confusionmat(a(randperm(length(a))),b);
%     fake(in,1) = (sum(j(eye(6)==1)))./(sum(sum(j)));
% end
% p = (sum(fake>=real)+1)/(nS+1);
% figure; histogram(fake,'FaceColor','k')
% yl = get(gca,'ylim'); hold on;
% plot([real real],yl,'r','LineWidth',3)
% title([label ' Forward Replays vs Theta, p = ' num2str(p)])
% 
% helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVSForwardReplays_' label '_max' sigmodulabel])
% 
% [~,a] = max(Fr,[],2);
% [~,b] = max(Th,[],2);
% j = confusionmat(a,b);
% real = (sum(j(eye(6)==1)))./(sum(sum(j)));
% 
% fake = NaN(nS,1);
% for in = 1:nS
%     j = confusionmat(a(randperm(length(a))),b);
%     fake(in,1) = (sum(j(eye(6)==1)))./(sum(sum(j)));
% end
% p = (sum(fake>=real)+1)/(nS+1);
% figure; histogram(fake,'FaceColor','k')
% yl = get(gca,'ylim'); hold on;
% plot([real real],yl,'r','LineWidth',3)
% title(['Theta vs Prospective FR, p = ' num2str(p)])
% 
% helper_saveandclosefig([savefolder '\ArmCombinations\ThetaVSProspectiveFR_max' sigmodulabel])

%%
if 1
corrdat = NaN(6,6,size(ThZ,1),3);
corrdatS = NaN(6,6,size(ThZ,1),size(ReShuffZ,3));
nS = size(ReShuffZ,3);
for icell = 1:size(ThZ,1)
    corrdat(:,:,icell,1) = repmat(ReZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    corrdat(:,:,icell,2) = repmat(ReZ(icell,:),[6 1])-(repmat(ThZ(icell,:),[6 1])');
    corrdat(:,:,icell,3) = repmat(ThZ(icell,:),[6 1])-(repmat(FrZ(icell,:),[6 1])');
    corrdatS(:,:,icell,:) = repmat(ReShuffZ(icell,:,:),[6 1])-(repmat(repmat(ThZ(icell,:),[6 1])',[1 1 size(ReShuffZ,3)]));
end

p = NaN(3,1);
x3 = [.5; 1.5;1.5; .5];
x2 = x3*ones(1,6)+ones(4,1)*(0:5);
y3 = [.5; .5 ; 1.5;1.5];
y2 = y3*ones(1,6)+ones(4,1)*(0:5);

xlab = {'Forward Replays (arm on, arm representing)','Forward Replays (arm on, arm representing)','Theta (arm on, arm representing)'};
ylab = {'Prospective FR (arm on, arm representing)','Theta (arm on, arm representing)','Prospective FR (arm on, arm representing)'};
datlab = {[label ' Forward Replays vs Prospective FR, p = '],[label ' Forward Replays vs Theta, p = '],'Theta vs Prospective FR, p = '};
for idattype = 1:3
    j = squeeze(nanmean(abs(corrdat(:,:,:,idattype)),3));
    real = (nansum(j(eye(6)==1)))./(nansum(nansum(j)));
    fake = NaN(nS,1);
    for in = 1:nS

        if idattype==2
            j2 = squeeze(nanmean(abs(corrdatS(:,:,:,in)),3));
        else
            j2 = squeeze(nanmean(abs(corrdat(randperm(size(corrdat,1)),:,:,idattype)),3));
        end
        fake(in,1) = (nansum(j2(eye(6)==1)))./(nansum(nansum(j2)));
    end
    p(idattype) = (sum(fake<=real)+1)/(nS+1);
    
    figure; hold on
    subplot(1,2,1); hold on
    imagesc(j); colormap gray; cb = colorbar; axis xy
    cb.Label.String = 'Difference in Z-scored value';
    patch(x2,y2,'r','FaceColor','none','EdgeColor','r','LineWidth',3)
    set(gca,'xtick',[1:6],'xticklabel',{'1,2','1,3','2,1','2,3','3,1','3,2'})
    set(gca,'ytick',[1:6],'yticklabel',{'1,2','1,3','2,1','2,3','3,1','3,2'})
    xlabel(xlab{idattype})
    ylabel(ylab{idattype})
    xlim([.5 6.5]);ylim([.5 6.5])
    set(gca,'FontSize',18)
    
    subplot(1,2,2)
    a1 = histogram(fake,'FaceColor','k');
    yl = get(gca,'ylim'); 
    hold on;
    a2 = plot([real real],yl,'r','LineWidth',3);
    xlabel('Sum of the Diagonal')
    ylabel('Count')
    legend([a1 a2],{'Shuffle','Data'})
    set(gca,'FontSize',18)
    
    suptitle([datlab{idattype} num2str(round(p(idattype),3,'significant')) ' n=' nn])   
    set(gcf,'Position',[1958         138        1523         669])
    helper_saveandclosefig([savefolder '\ArmCombinations\Allcombinations_Type_' num2str(idattype) '_' label sigmodulabel])
end
end
    
