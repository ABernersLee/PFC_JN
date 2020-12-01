function NonLocalThetaDifference_SummaryFigs(dirs,cutoff,savelab,igroup,savefolder)
cd(dirs.homedir)
d2 = dir('*.mat');
pcells1 = []; pcells2 = []; celldat = []; shuffdat = []; armsig = []; label = 'RP';


for id = 1:size(d2,1)
    thisdir = d2(id).name;
%     make_prospective_theta_triggered(thisdir,cutoff)
    load(thisdir,'prospectiveTheta_pcells_ranksum','prospectiveTheta_pcells_perms','prospectiveTheta_celldat','prospectiveTheta_shuffdat')
    load(thisdir,'other_cells_touse')
    prospectiveTheta_pcells_ranksum = prospectiveTheta_pcells_ranksum(other_cells_touse(:,igroup),:);
    prospectiveTheta_pcells_perms = prospectiveTheta_pcells_perms(other_cells_touse(:,igroup),:);   
    prospectiveTheta_celldat = prospectiveTheta_celldat(:,:,other_cells_touse(:,igroup));
    prospectiveTheta_shuffdat = prospectiveTheta_shuffdat(:,:,other_cells_touse(:,igroup),:);
    
%     numshuff = size(prospectiveTheta_shuffdat,4);
    pcells1 = cat(1,pcells1,prospectiveTheta_pcells_ranksum);
    pcells2 = cat(1,pcells2,prospectiveTheta_pcells_perms);
    
    
    numshuff = size(prospectiveTheta_shuffdat,4);
        
    if numshuff>1
        celldat1 = NaN(size(prospectiveTheta_pcells_perms,1),2);
        shuffdat1 = NaN(size(prospectiveTheta_pcells_perms,1),2,numshuff);
        if ~isempty(prospectiveTheta_celldat)
        for icell = 1:size(prospectiveTheta_pcells_perms,1)
            celldat1(icell,:) = [nanmean(prospectiveTheta_celldat(prospectiveTheta_celldat(:,2,icell)==0,1,icell)) ...
                nanmean(prospectiveTheta_celldat(prospectiveTheta_celldat(:,2,icell)==1,1,icell))];
            for ii = 1:numshuff
                shuffdat1(icell,:,ii) = [nanmean(prospectiveTheta_shuffdat(prospectiveTheta_shuffdat(:,2,icell,ii)==0,1,icell,ii)) ...
                nanmean(prospectiveTheta_shuffdat(prospectiveTheta_shuffdat(:,2,icell,ii)==1,1,icell,ii))];
            end
        end
        end
        celldat = cat(1,celldat,celldat1);        
        shuffdat = cat(1,shuffdat,shuffdat1);
    end
    load(d2(id).name,[label '_SSDarm'],[label '_moduarm'],[label '_moduarm2'],[label '_pSSDarm'],[label '_SSDarm2'],[label '_Cand_sig_modu_include'])
    eval(['aa = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include ' label '_SSDarm2];'])
    armsig = cat(1,armsig,aa(other_cells_touse(:,igroup),:));
end

%%


realdat = celldat(:,2)-celldat(:,1);
shuff = squeeze(shuffdat(:,2,:)-shuffdat(:,1,:));
figure; hold on
violins = violinplot1([realdat mean(shuff,2)]);
violins(1).ViolinColor = 'k'; 
violins(1).EdgeColor = 'k'; 
violins(1).MedianColor = 'k';
violins(2).ViolinColor = 'k'; 
violins(2).EdgeColor = 'k'; 
violins(2).MedianColor = 'k';
% plot(1,realdat,'+','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',20)
% errorbar(1,mean(realdat),std(realdat)./sqrt(sum(~isnan(realdat))),'k','LineWidth',3)
p1 = signrank(celldat(:,1),celldat(:,2));
% shuff2 = shuff(:);
p2 = (sum(mean(shuff)>=mean(realdat))+1)/(size(shuff,2)+1);
% plot(2,mean(shuff,2),'+','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',20)
% errorbar(2,mean(mean(shuff,2)),std(mean(shuff,2))./sqrt(sum(~isnan(mean(shuff,2)))),'k','LineWidth',3)

if p2<.001
    yl = get(gca,'ylim');
    plot([1 2],[yl(2)*1.05 yl(2)*1.05],'k','LineWidth',3)
    plot([1.47 1.5 1.53],[yl(2)*1.15 yl(2)*1.15 yl(2)*1.15],'r*','MarkerSize',10)
elseif p2<.05
    yl = get(gca,'ylim');
    plot([1 2],[yl(2)*1.05 yl(2)*1.05],'k','LineWidth',3)
    plot([1.5],[ yl(2)*1.15],'r*','MarkerSize',10)
end
if p1<.01
    yl = get(gca,'ylim');
    plot([.98 1.02],[yl(2)*1.2 yl(2)*1.2],'r*','MarkerSize',10)
    
elseif p1<.05    
    yl = get(gca,'ylim');
    plot(1,yl(2)*1.2,'r*','MarkerSize',10)
end
ylim([yl(1) yl(2)*2])
xlim([.7 2.3])
set(gca,'xtick',1:2,'xticklabel',{'Data';'Shuffle'})
ylabel('mPFC firing rate difference between theta sweeps of two non-local arms')
% title(['ranksum p = ' num2str(p1) ',permutation test p = ' num2str(p2)])


axes('Position',[.7 .72 .2 .2])
box on
hold on
pe = pie([sum(~isnan(pcells2)) sum(pcells2<.05)]);
pe(1).FaceColor = 'w';
pe(3).FaceColor = 'k';
p= 1-binocdf(sum(pcells2<.05)-1,sum(~isnan(pcells2)),.05);
% title(['p = ' num2str(p)])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])

set(gcf,'Position',[2193         151         954         732])

%%

if ~isfolder([savefolder '\Theta\'])
    mkdir([savefolder '\Theta\'])
end
if ~isfolder([savefolder '\Figure3\'])
    mkdir([savefolder '\Figure3\'])
end
set(gcf,'renderer','Painters')    
helper_savefig([savefolder '\Theta\NonLocalThetaDifference_pfc_SummaryFigs2_' num2str(cutoff) '_' savelab])
helper_saveandclosefig([savefolder '\Figure3\NonLocalThetaDifference_pfc_SummaryFigs2_' num2str(cutoff) '_' savelab])
%%
figure; hold on
% h1 = histc(realdat,-6:.5:6); 
[h2,c] = hist(nanmean(shuff),50); 
% b = bar([h1./sum(h1) h2./sum(h2)],'LineWidth',2,'BarWidth',1);
% b(2).FaceColor = [0 0 0];
% b(1).FaceColor = [1 1 1];
% bar(-6:.1:6,h1./sum(h1),'FaceColor','k','FaceAlpha',.5,'EdgeAlpha',1)
bar(c,h2./sum(h2),'FaceColor','w','FaceAlpha',.5,'LineWidth',2)
yl = get(gca,'ylim');
plot([nanmean(realdat) nanmean(realdat)],yl,'k','LineWidth',3)
text(.35,.06,['p = ' num2str(round(p2,2,'significant'))])
xlabel('Firing Rate Difference Between Theta Sweeps')
ylabel('Count')
set(gca,'FontName','Arial','FontSize',18)

axes('Position',[.16 .65 .2 .2])
box on
hold on
pe = pie([sum(~isnan(pcells2)) sum(pcells2<.05)]);
pe(1).FaceColor = 'w';
pe(3).FaceColor = 'k';
p= 1-binocdf(sum(pcells2<.05)-1,sum(~isnan(pcells2)),.05);
title(['Cells Significantly Modulated' newline ' by Direction of Theta Sweep'])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])
set(gcf,'Position',[  680   388   816   590])

set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_savefig([savefolder '\Theta\NonLocalThetaDifference_pfc_SummaryFigs3_' num2str(cutoff) '_' savelab])
helper_saveandclosefig([savefolder '\Figure3\NonLocalThetaDifference_pfc_SummaryFigs3_' num2str(cutoff) '_' savelab])
%%
p2 = (sum(mean(shuff(armsig(:,1)<.05,:))>=mean(realdat(armsig(:,1)<.05)))+1)/(size(shuff(armsig(:,1)<.05,:),2)+1);

figure; hold on
[h2,c] = hist(nanmean(shuff(armsig(:,1)<.05,:)),50); 
bar(c,h2./sum(h2),'FaceColor','w','FaceAlpha',.5,'LineWidth',2)
yl = get(gca,'ylim');
plot([nanmean(realdat(armsig(:,1)<.05)) nanmean(realdat(armsig(:,1)<.05))],yl,'k','LineWidth',3)
text(.35,.06,['p = ' num2str(round(p2,2,'significant'))])
xlabel('Firing Rate Difference Between Theta Sweeps')
ylabel('Count')
set(gca,'FontName','Arial','FontSize',18)

axes('Position',[.16 .65 .2 .2])
box on
hold on
pe = pie([sum(~isnan(pcells2(armsig(:,1)<.05))) sum(pcells2(armsig(:,1)<.05)<.05)]);
pe(1).FaceColor = 'w';
pe(3).FaceColor = 'k';
p= 1-binocdf(sum(pcells2(armsig(:,1)<.05)<.05)-1,sum(~isnan(pcells2(armsig(:,1)<.05))),.05);
title(['Cells Significantly Modulated' newline ' by Direction of Theta Sweep'])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])
set(gcf,'Position',[  680   388   816   590])

set(gca,'FontName','Arial')
set(gcf,'renderer','Painters')
helper_savefig([savefolder '\Theta\NonLocalThetaDifference_pfc_SummaryFigs3_SigOnly_' num2str(cutoff) '_' savelab])
helper_saveandclosefig([savefolder '\Figure3\NonLocalThetaDifference_pfc_SummaryFigs3_SigOnly_' num2str(cutoff) '_' savelab])

%%
figure; hold on

% realdat = celldat(:,2)-celldat(:,1);
plot(1,celldat(:,1),'+','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',20)
plot(2,celldat(:,2),'+','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',20)
errorbar(1,mean(celldat(:,1)),std(celldat(:,1))./sqrt(sum(~isnan(celldat(:,1)))),'k','LineWidth',3)
errorbar(2,mean(celldat(:,2)),std(celldat(:,2))./sqrt(sum(~isnan(celldat(:,2)))),'k','LineWidth',3)
p1 = signrank(celldat(:,1),celldat(:,2));
plot([ones(size(celldat,1),1) 2*ones(size(celldat,1),1)]', celldat','k')
if p1<.001
    yl = get(gca,'ylim');
    plot([1 2],[yl(2)*1.05 yl(2)*1.05],'k')
    plot([1.43 1.5 1.57],[yl(2)*1.15 yl(2)*1.15 yl(2)*1.15],'r*','MarkerSize',10)
end
ylim([yl(1) yl(2)*1.7])
xlim([.7 2.3])
title(['ranksum p = ' num2str(p1)])

axes('Position',[.7 .72 .2 .2])
box on
hold on
pie([sum(~isnan(pcells1)) sum(pcells1<.05)])
p= 1-binocdf(sum(pcells1<.05)-1,sum(~isnan(pcells1)),.05);
% title(['p = ' num2str(p)])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])

axes('Position',[.3 .32 .2 .2])
box on
hold on
pie([sum(~isnan(pcells2)) sum(pcells2<.05)])
p= 1-binocdf(sum(pcells2<.05)-1,sum(~isnan(pcells2)),.05);
% title(['p = ' num2str(p)])
xlim([-1.7 1.7])
ylim([-1.5 1.5])
set(gca,'xtick',[],'ytick',[])


helper_saveandclosefig([savefolder '\Theta\NonLocalThetaDifference_pfc_SummaryFigs1_' num2str(cutoff) '_' savelab])


%% is the change in FR of each cell to different armed theta cycles in hp related to whether it gets arm input? yes
celldat = []; %pcells = []; 
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prospectiveTheta_pcells_perms','prospectiveTheta_celldat')
    if isempty(prospectiveTheta_celldat)
        celldat1 = squeeze(nanmean(NaN(1,2,size(prospectiveTheta_pcells_perms,1)),1));
    else
%     pcells = cat(1,pcells,prospectiveTheta_pcells);
    celldat1 = squeeze(nanmean(prospectiveTheta_celldat,1));        
    end
    celldat = cat(2,celldat,celldat1);
end
celldatdiff = abs(diff(celldat,[],1)); %./sum(celldat);
diffperc = celldatdiff';

touse = true(size(armsig,1),1); %armsig(:,5)==1;
armsig2 = armsig(touse,:);
diffdat = nanmean(diffperc(touse,:),2); %,[],2); %,2);
figure; hold on;
set(gcf,'Position',[ 2073         146        1202         744])
subplot(3,1,1), hold on
h(1) = histogram(diffdat(armsig2(:,1)>=.05),20,'Normalization','probability');
h(2) = histogram(diffdat(armsig2(:,1)<.05),20,'Normalization','probability');
yl = get(gca,'ylim');
plot([nanmedian(diffdat(armsig2(:,1)>=.05)) nanmedian(diffdat(armsig2(:,1)>=.05))],yl,'b-','LineWidth',3)
plot([nanmedian(diffdat(armsig2(:,1)<.05)) nanmedian(diffdat(armsig2(:,1)<.05))],yl,'r-','LineWidth',3)

p = ranksum(diffdat(armsig2(:,1)<.05),diffdat(armsig2(:,1)>=.05),'tail','right');
legend(h,'Not Mod','Arm Mod')
tt = text(.2,.2,['p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('HP Non-local Theta Modulation')
ylabel('Probability')
set(gca,'FontSize',18)


subplot(3,1,2)
plot(armsig2(:,2),diffdat,'.k','MarkerSize',15)
[r,p] = corr(armsig2(:,2),diffdat,'rows','complete');
tt = text(15,.8,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('armSSD')
ylabel('Non-local Theta Mod')
set(gca,'FontSize',18)


subplot(3,1,3)
plot(log(armsig2(:,2)),diffdat,'.k','MarkerSize',15)
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
[r,p] = corr(log(armsig2(:,2)),diffdat,'rows','complete');
tt = text(-10,.9,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('log(armSSD)')
ylabel('Non-local Theta Mod')
set(gca,'FontSize',18)
set(gcf,'Position',[2028         -73         624        1064])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Theta\NonLocalThetaDifference_pfc_Vs_Armsig_' num2str(cutoff) '_' savelab])

%%


figure; hold on
plot(log(armsig2(:,2)),diffdat,'.k','MarkerSize',15)
x1 = [min(log(armsig2(:,2))) max(log(armsig2(:,2)))];
touse = ~isnan(armsig2(:,2)) & ~isnan(diffdat);
b = polyfit(log(armsig2(touse,2)),diffdat(touse),1);
y2 = polyval(b,x1);
plot(x1,y2,'r','LineWidth',3)
[r,p] = corr(armsig2(:,2),diffdat,'rows','complete');
tt = text(-9,.9,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
[r,p] = corr(log(armsig2(:,2)),diffdat,'rows','complete');
tt = text(-9,.8,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('log(armSSD)')
ylabel('Non-local Theta Mod')
set(gca,'FontSize',18)
set(gcf,'Position',[ 2211          47        1087         796])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Theta\NonLocalThetaDifference_pfc_Vs_Armsig2_' num2str(cutoff) '_' savelab])

%%


figure; hold on
plot(log(armsig2(:,end)),diffdat,'.k','MarkerSize',15)
x1 = [min(log(armsig2(:,end))) max(log(armsig2(:,end)))];
touse = ~isnan(armsig2(:,end)) & ~isnan(diffdat);
b = polyfit(log(armsig2(touse,2)),diffdat(touse),1);
y2 = polyval(b,x1);
plot(x1,y2,'r','LineWidth',3)
[r,p] = corr(armsig2(:,end),diffdat,'rows','complete');
tt = text(-9,.9,['Lin, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
[r,p] = corr(log(armsig2(:,end)),diffdat,'rows','complete');
tt = text(-9,.8,['Log, r = ' num2str(round(r,2,'significant')) ' p = ' num2str(round(p,2,'significant'))]);
if p<.05
    tt.Color = 'r';
    tt.FontSize = 14;
end
xlabel('log(normalized armSSD)')
ylabel('Non-local Theta Mod')
set(gca,'FontSize',18)
set(gcf,'Position',[ 2211          47        1087         796])
helper_saveandclosefig([savefolder '\Theta\NonLocalThetaDifference_pfc_Vs_Armsig3_normalized_' num2str(cutoff) '_' savelab])