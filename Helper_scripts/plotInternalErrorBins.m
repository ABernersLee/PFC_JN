function plotInternalErrorBins(savefolder)

EstBin = .04;
velcutoff = 5;
cellcutoff = 4;
spikecutoff = 5;
numshuff = 1000;
% errcutoff = 0.5;
% binsize = 20;
binerror = [0.5 20 8 4; 1.5 20 8 12; 2.5 20 8 20; 0.5 15 6 3; 1.5 15 6 9; 0.5 12 5 2.5; 1.5 12 5 7; 4.5 20 8 36]; 
PFCdelay = 0.04;
diffscore = NaN(2,size(binerror,1));
pval = NaN(size(binerror,1),1);
f = 'F:\XY_matdata\Figures\ForPaperReviews\AllCells_AllArmEvents\Figure3';
for ibinerror = 1:size(binerror,1)
    errcutoff = binerror(ibinerror,1);
    binsize = binerror(ibinerror,2);
    figlab = ['_BinSize' num2str(EstBin) '_VelCutoff' num2str(velcutoff) '_cellcutoff' num2str(cellcutoff) '_spikecutoff' num2str(spikecutoff) '_errcutoff' num2str(errcutoff) '_numshuff' num2str(numshuff) '_binsize' num2str(binsize) '_PFCdelay' num2str(PFCdelay)];
    figname = [ f '\BarPlots_ReplayVSLocalNonLocal_NL_Shuffle_CellByCell_Sig Cells_Norm_High Vel' figlab '.fig'];
    if exist(figname,'file')
        open(figname)
        ag = gca;
        diffscore(1,ibinerror) = ag.Children(1).YData(2)-ag.Children(1).YData(1);        
        neurons = NaN(size(ag.Children,1)-1,1);        
        for ineuron = 2:size(ag.Children,1)
            neurons(ineuron-1) = ag.Children(ineuron).YData(2)-ag.Children(ineuron).YData(1);
        end
        diffscore(2,ibinerror) = std(neurons)./sqrt(sum(~isnan(neurons)));
        pval(ibinerror) = str2num(ag.Title.String(end-4:end));
        close all
    end
end
%%

[~,ord] = sort(binerror(:,4));
xtlab = cell(size(binerror,1),1);
for ibinerror = 1:size(binerror,1)
   xtlab{ibinerror,1} = [num2str(binerror(ord(ibinerror),4))]; % '-' num2str(binerror(ord(ibinerror),3)) ''];  
end

%need errorbars
figure; hold on;
plot(binerror(ord,4),diffscore(1,ord),'ok','MarkerSize',15)
errorbar(binerror(ord,4),diffscore(1,ord),diffscore(2,ord),'k','LineWidth',1);
% xl = get(gca,'xlim');
plot([min((binerror(:,4)))-1 max(binerror(:,4))+.01],[0 0],'k--')
yl = get(gca,'ylim');

for idelay = 1:size(binerror,1)
    if pval(idelay)<.01
        text(binerror(idelay,4),yl(2),'*','FontSize',30)
    elseif pval(idelay)<.05
        text(binerror(idelay,4),yl(2),'*','FontSize',20)
    elseif pval(idelay)<.1
        text(binerror(idelay,4),yl(2),'~','FontSize',20)
    end
end
set(gca,'ylim',[yl(1) yl(2)*1.1],'xlim',[min((binerror(:,4)))-1 max(binerror(:,4))+.01])
set(gca,'xtick',binerror(ord,4),'xticklabel',xtlab)
% set(gca,'FontSize',18)
ylabel('Real-Shuffle')
xlabel('Error cutoff (cm)')
set(gcf,'Position',[ 1722         186        1473         545])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure3\AcrossBinErrorSizes'])