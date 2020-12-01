function Figure3_make_4sets_barplots(dirs,igroup,savefolder,toplot)
%%
cd(dirs.homedir)
d2 = dir('*.mat');
repall = []; sigcells = []; othall = [];
othlabel = {'Behavioral FR';'Local';'Non-Local'};
% newcol = [75 0 130;255 130 0;34 139 34]/255;
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'ArmBarPlot_LocalNonLocal_MeanSEM','other_cells','RP_moduarm','RP_pSSDarm','other_cells_touse')
    replay = RP_moduarm';
    
    
    fr = NaN(size(other_cells,1),3);
    
    %from linear fields
    load(thisdir,'InFR','OutFR','armposindex')    
    FR = InFR(other_cells,:)+OutFR(other_cells,:);
    for iarm = 1:3
        fr(:,iarm) = nanmean(FR(:,armposindex(:,iarm)),2);
    end    
    
    %from open fields
%     load(thisdir,'OpenFR','pos','armpos')    
%     binsize = 20;    
%     pos2 = pos(:,2:3);
%     pos2(:,1) = ceil(pos2(:,1)/binsize);
%     pos2(:,2) = ceil(pos2(:,2)/binsize);
%     bp = pos2-min(pos2)+1;    
%     pind = max(bp);
%     ba = armpos;
%     x = sub2ind(pind,bp(:,1),bp(:,2)); %real linear indicies    
%     ind1 = setdiff(unique(x(ba==1)),unique([x(ba==2);x(ba==3)]));
%     ind2 = setdiff(unique(x(ba==2)),unique([x(ba==1);x(ba==3)]));
%     ind3 = setdiff(unique(x(ba==3)),unique([x(ba==1);x(ba==2)]));
%     for icell = 1:length(other_cells)
%         dat = OpenFR(:,:,other_cells(icell));        
%         datnan = isnan(dat);
%         dat(datnan) = 0;
%         dat = filter2(Two_D_Filter,dat);
%         dat(datnan) = NaN;
%         fr(icell,:) = [nanmean(dat(ind1)) nanmean(dat(ind2)) nanmean(dat(ind3))];
%     end
    
    local = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,1,1);
    nonlocal = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,2,1);
    clear FR InFR OutFR RP_moduarm ArmBarPlot_LocalNonLocal_MeanSEM
    % use replay, fr and local and lonlocal
%     replay = replay./sum(replay,2);
    if toplot
        for icell = 1:length(other_cells)
            if RP_pSSDarm(icell)>=.05
                continue
            end
            figure; hold on
            subplot(2,2,1); hold on
            bar(replay(icell,:),'FaceColor','w')
            subplot(2,2,2); hold on
            bar(fr(icell,:),'FaceColor','w')
            subplot(2,2,3); hold on
            bar(local(icell,:),'FaceColor','w')
            subplot(2,2,4); hold on
            bar(nonlocal(icell,:),'FaceColor','w')
            set(gcf,'renderer','Painters')
            suptitle([thisdir(1) ' ' thisdir(3:end-4) ', Cell' num2str(other_cells(icell)) ' Raw'])
            helper_saveandclosefig([savefolder '\Figure3\related\BarPlots_ByArm_' thisdir(1:end-4) '_Cell' num2str(other_cells(icell)) '_raw_new'])            

            figure; hold on
            subplot(2,2,1); hold on
            bar(replay(icell,:)./sum(replay(icell,:)),'FaceColor','w')
            subplot(2,2,2); hold on
            bar(fr(icell,:)./sum(fr(icell,:)),'FaceColor','w')
            subplot(2,2,3); hold on
            bar(local(icell,:)./sum(local(icell,:)),'FaceColor','w')
            subplot(2,2,4); hold on
            bar(nonlocal(icell,:)./sum(nonlocal(icell,:)),'FaceColor','w')
            set(gcf,'renderer','Painters')
            suptitle([thisdir(1) ' ' thisdir(3:end-4) ', Cell' num2str(other_cells(icell)) ' Norm'])
            helper_saveandclosefig([savefolder '\Figure3\related\BarPlots_ByArm_' thisdir(1:end-4) '_Cell' num2str(other_cells(icell)) '_norm_new'])            
        end
    end
    repall = cat(1,repall,replay(other_cells_touse(:,igroup),:));
    othall = cat(1,othall,cat(3,fr(other_cells_touse(:,igroup),:),local(other_cells_touse(:,igroup),:),nonlocal(other_cells_touse(:,igroup),:)));
    sigcells = cat(1,sigcells,RP_pSSDarm(other_cells_touse(:,igroup))<.05);
end
%%
siglab = {'All Cells';'Sig Cells'};
for isig = 2
r = repall;
o = othall;
if isig == 1; rr = r; oo = o; 
elseif isig == 2
rr = r(sigcells==1,:); oo = o(sigcells==1,:,:);
end
% rr(sum(rr<0,2)==3,:) = [];
% rr(sum(rr<0,2)==3,:) = -rr(sum(rr<0,2)==3,:);
rvalues = NaN(size(rr,1),3);
col = {'r';'b';'y'};


 for ioth = 1:3    
    for icell = 1:size(rr,1); rvalues(icell,ioth) = corr(rr(icell,:)',squeeze(oo(icell,:,ioth))','type','Spearman'); end
        
    figure; hold on; 
%     histogram(rvalues(:,ioth),4,'FaceColor',col{ioth})
    h = histc(rvalues(:,ioth),-1:.5:1);     
    b = bar([-1 -.5 .5 1],h(1:end-1,:),'LineWidth',3,'FaceColor','w');
    b(1).EdgeColor = col{ioth};
%     plot([mean(rvalues(:,ioth)) mean(rvalues(:,ioth))],[0 12],'LineWidth',2,'Color','k')
%     plot([median(rvalues(:,ioth)) median(rvalues(:,ioth))],[0 12],'--','LineWidth',2,'Color','k')
%     ylim([0 15])
    xlabel('Spearman''s rho')
    ylabel('PFC cells')
    p = signrank(rvalues(:,ioth));
    title(['Modulation by ' othlabel{ioth} ' p = ' num2str(round(p,2,'significant'))])
    set(gcf,'renderer','Painters')
    
    helper_saveandclosefig([savefolder '\Figure3\BarPlots_Bars_ReplayVS_' ...
               othlabel{ioth} ' ' siglab{isig} '_new'])            
 end


 for ioth = 1:2
  figure; hold on; 
  

    bar(1,mean(rvalues(:,ioth)),'LineWidth',3,'FaceColor','w');
    errorbar(1,mean(rvalues(:,ioth)),std(rvalues(:,ioth))./sqrt(size(rvalues,1)),'k','LineWidth',3)
    bar(2,mean(rvalues(:,3)),'LineWidth',3,'FaceColor','w');
    errorbar(2,mean(rvalues(:,3)),std(rvalues(:,3))./sqrt(size(rvalues,1)),'k','LineWidth',3)  
    
%     xlabel('Spearman''s rho')
    set(gca,'xtick',1:2,'xticklabel',{othlabel{ioth};othlabel{3}})
    ylabel('Spearman''s rho for each PFC cell')
    yl = get(gca,'ylim');
%     ylim([yl(1) .6])
    
    p = signrank(rvalues(:,ioth),rvalues(:,3),'tail','left');
    title([' p = ' num2str(round(p,2,'significant'))])
    set(gcf,'renderer','Painters')
    
    helper_saveandclosefig([savefolder '\Figure3\BarPlots_Means_ReplayVS_' ...
               othlabel{ioth}  ' ' siglab{isig} '_new'])            
 end

 %%
 
%  figure; hold on
 h = NaN(5,3);
 for ioth = 1:3
     h(:,ioth) = histc(rvalues(:,ioth),-1:.5:1);     
%      bar([-1 -.5 .5 1],h(1:end-1,ioth),'LineWidth',3,'FaceColor','w')
 end
 figure; hold on; 
 b = bar(h(1:end-1,:),'LineWidth',3,'FaceColor','w');
 b(1).EdgeColor = 'r';
 b(2).EdgeColor = 'b';
 b(3).EdgeColor = 'y';
 legend(othlabel,'Location','northwest')
 helper_saveandclosefig([savefolder '\Figure3\BarPlots_Bars_ReplayVS_' siglab{isig}])
 
 figure; hold on
 for ioth = 1:3
     c = cdfplot(rvalues(:,ioth));
     c.LineWidth = 3;
 end
 legend(othlabel,'Location','northwest')
 helper_saveandclosefig([savefolder '\Figure3\BarPlots_cdf_ReplayVS_' siglab{isig} '_new'])   
end
 %%
zlab = {'Across Cells';'Median Subtracted';'Zscored';'Norm'};
siglab = {'All Cells';'Sig Cells'};
for iz = 1:4
    if iz == 1
        r = repall;
        o = othall;
    elseif iz == 2
        r = repall-median(repall,2);
        o = othall-median(othall,2);
    elseif iz == 3
        r = zscore(repall,[],2);
        o = zscore(othall,[],2);
    elseif iz == 4
        r = (repall-min(repall,[],2))./range(repall,2);
        o = (othall-min(othall,[],2))./range(othall,2);
    end
    
    for isig = 2
       if isig == 1
           rr = r; oo = o;
       elseif isig == 2
           rr = r(sigcells==1,:); oo = o(sigcells==1,:,:);
       end
       dat1 = rr(:); 
       for ioth = 1:3
           dat2 = oo(:,:,ioth); dat2 = dat2(:);
           figure; hold on
           plot(dat1,dat2,'ok','MarkerSize',10,'LineWidth',3)
           [rho,p] = corr(dat1,dat2,'rows','complete','type','Spearman'); %,'type','Kendall');
           xlabel('Modulation by Replay')
           ylabel(['Modulation by ' othlabel{ioth}])
           set(gca,'FontSize',18)
           xl = get(gca,'xlim'); yl = get(gca,'ylim');
           t1 = text(xl(2)*.5,yl(2)*.9,['r = ' num2str(round(rho,2,'significant')) ' p = ' num2str(round(p,2,'significant'))],'FontSize',18);
           if p<.05
                t1.Color = 'r';
           end
           lm = polyfit(dat1,dat2,1);           
           sh = range(xl)/10;
           xll = [xl(1)+sh xl(2)-sh];
           y2 = polyval(lm,xll);           
           plot(xll,y2,'r','LineWidth',3)
           title([zlab{iz} ' ' siglab{isig} ' - Replay vs ' othlabel{ioth}])
           set(gcf,'Position',[680   431   928   547])
           set(gcf,'renderer','Painters')           
           helper_saveandclosefig([savefolder '\Figure3\BarPlots_Scatter_ReplayVS' ...
               othlabel{ioth} ' ' siglab{isig} ' ' zlab{iz} '_new'])            
       end
    end    
end

