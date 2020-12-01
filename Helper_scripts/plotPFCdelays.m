function plotPFCdelays(savefolder)

EstBin = .04;
velcutoff = 5;
cellcutoff = 4;
spikecutoff = 5;
numshuff = 1000;
errcutoff = 0.5;
binsize = 20;
PFCdelays = [.01:.01:.07 .09 .11];
diffscore = NaN(2,length(PFCdelays));
pval = NaN(length(PFCdelays),1);
f = 'F:\XY_matdata\Figures\ForPaperReviews\AllCells_AllArmEvents\Figure3';
for idelay = 1:length(PFCdelays)
    PFCdelay = PFCdelays(idelay);
    figlab = ['_BinSize' num2str(EstBin) '_VelCutoff' num2str(velcutoff) '_cellcutoff' num2str(cellcutoff) '_spikecutoff' num2str(spikecutoff) '_errcutoff' num2str(errcutoff) '_numshuff' num2str(numshuff) '_binsize' num2str(binsize) '_PFCdelay' num2str(PFCdelay)];
    figname = [ f '\BarPlots_ReplayVSLocalNonLocal_NL_Shuffle_CellByCell_Sig Cells_Norm_High Vel' figlab '.fig'];
    if exist(figname,'file')
        open(figname)
        fg = gcf;
        diffscore(1,idelay) = fg.Children(1).Children(1).YData(2)-fg.Children(1).Children(1).YData(1);
        neurons = NaN(size(fg.Children(1).Children,1)-1,1);
        for ineuron = 2:size(fg.Children(1).Children,1)
            neurons(ineuron-1) = fg.Children(1).Children(ineuron).YData(2)-fg.Children(1).Children(ineuron).YData(1);
        end
        diffscore(2,idelay) = std(neurons)./sqrt(sum(~isnan(neurons)));
        pval(idelay) = str2num(fg.Children(1).Title.String(end-4:end));
        close all
    end
end
%%

figure; hold on;
plot(PFCdelays,diffscore(1,:),'ok','MarkerSize',15)
errorbar(PFCdelays,diffscore(1,:),diffscore(2,:),'k','LineWidth',1)
xl = get(gca,'xlim');
plot(xl,[0 0],'k--')
yl = get(gca,'ylim');

for idelay = 1:length(PFCdelays)
    if pval(idelay)<.01
        text(PFCdelays(idelay),yl(2),'*','FontSize',30)
    elseif pval(idelay)<.05
        text(PFCdelays(idelay),yl(2),'*','FontSize',20)
    elseif pval(idelay)<.1
        text(PFCdelays(idelay),yl(2),'~','FontSize',20)
    end
end
set(gca,'ylim',[yl(1) yl(2)*1.1],'xlim',[0 max(PFCdelays)+.01])
ylabel('Real-Shuffle')
xlabel('Delay of PFC (sec)')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure3\AcrossDelays_All'])
%%
PFCdelays = PFCdelays(1:7);
diffscore = diffscore(:,1:7);
figure; hold on;
plot(PFCdelays,diffscore(1,:),'ok','MarkerSize',15)
errorbar(PFCdelays,diffscore(1,:),diffscore(2,:),'k','LineWidth',1)
xl = get(gca,'xlim');
plot(xl,[0 0],'k--')
yl = get(gca,'ylim');

for idelay = 1:length(PFCdelays)
    if pval(idelay)<.01
        text(PFCdelays(idelay),yl(2),'*','FontSize',30)
    elseif pval(idelay)<.05
        text(PFCdelays(idelay),yl(2),'*','FontSize',20)
    elseif pval(idelay)<.1
        text(PFCdelays(idelay),yl(2),'~','FontSize',20)
    end
end
set(gca,'ylim',[yl(1) yl(2)*1.1],'xlim',[0 max(PFCdelays)+.01])
ylabel('Real-Shuffle')
xlabel('Delay of PFC (sec)')
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure3\AcrossDelays_Short'])