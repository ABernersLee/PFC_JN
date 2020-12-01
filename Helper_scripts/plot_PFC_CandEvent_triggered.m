function plot_PFC_CandEvent_triggered(thisdir,label,withmodpatch,smoothsize,savefolder,id)

load(thisdir,[label '_CandEventAssign'],[label '_replay_singlebothjoint'],[label '_replay_replayarm'],[label '_PFCcandspikes_binned'],[label '_PFCcandspikes_list'],...
    'other_cells',[label '_Cand_sig_modu_include'],[label '_Cand_p'],[label '_pSSDreplay']) %[label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list']

dday = [2;4;5;1;2;3;6;1;3;4;1];
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
eval(['p_SSD1 = ' label '_pSSDreplay;'])
eval(['p_SSD2 = ' label '_Cand_p;'])
eval(['PFCcandspikes_binned = ' label '_PFCcandspikes_binned;'])
eval(['PFCcandspikes_list = ' label '_PFCcandspikes_list;'])
% eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
% eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['CandEventAssign = ' label '_CandEventAssign;'])
clear([label '_PFCcandspikes_list'],[label '_PFCcandspikes_binned'],[label '_Cand_sig_modu_include'])
pfc = other_cells(Cand_sig_modu_include(:,3)==1); clear other_cells;
PFCcandspikes_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];
% PFCreplayspikes_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];
p_SSD1(Cand_sig_modu_include(:,3)==0,:) = [];
p_SSD2(Cand_sig_modu_include(:,3)==0,:) = [];
CandEventAssign(end) = [];

%load  p_SSD2

if strcmp(label,'SD')
    xlab = 'Hippocampal Spike Density Event';
    ylab = 'Spike Density Event Number';
elseif strcmp(label,'RP')
    xlab = 'Hippocampal Ripple';
    ylab = 'Ripple Number';
end
% binsize = .002;
% window = [-.5 .5];
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
ind = [window(1):binsize:window(2)];
% smoothsize = 20;
% withmodpatch = true;
% tousereplays = ~isnan(replayarm) & singlebothjoint~=2;
touseR22 = CandEventAssign>0;
sigma = 1.5;
plotwindow = [-.6 .6];       
newcol = [0 0 0;1 0 0];


modind = [ind(ind>=0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<=-.1)];
for icell = 1:length(pfc)
    toadd = 0;
    figure; hold on
    for ievent = 1:2
        
        
        if ievent==1
%             touseR = ~tousereplays; 
            touseR2 = ~touseR22;
        elseif ievent==2
%             touseR = tousereplays; touseR2 = false(size(touseR22));
            touseR2 = touseR22;
        end
%         dl1 = PFCreplayspikes_list;
%         dl1(ismember(PFCreplayspikes_list(:,3),find(~touseR)),:) = [];        
%         dl2 = PFCcandspikes_list;
%         dl2(ismember(PFCcandspikes_list(:,3),find(~touseR2)),:) = [];        
%         dl = [dl1;dl2];
        dl = PFCcandspikes_list;
        dl(ismember(PFCcandspikes_list(:,3),find(~touseR2)),:) = [];        
        
%         pb = cat(3,PFCreplayspikes_binned(icell,:,touseR),PFCcandspikes_binned(icell,:,touseR2));
        pb = PFCcandspikes_binned(icell,:,touseR2);
        dat1 = squeeze(mean(pb,3))./binsize;
%         sem1 = squeeze(std(pb./binsize,[],3))./sqrt(sum(touseR)+sum(touseR2));
        sem1 = squeeze(std(pb./binsize,[],3))./sqrt(sum(touseR2));
        sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
        gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
        dat = conv(dat1',gaussFilter,'same')';
        sem = conv(sem1',gaussFilter,'same')';
        rev = dat+sem;
        fwd = dat-sem;     


    %     rev = dat+sem;
    %     sm = smoothts([dat dat dat dat],'b',smoothsize);

%         modind = [ind(ind>0 & ind<=.2)];
%         baseind = [ind(ind>=-.5 & ind<-.1)];

         dl2 = dl(dl(:,2)==pfc(icell),:);
%         [~,~,jj] = unique(dl2(:,3));
        jj = unique(dl2(:,3));
        test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
        test(test<0) = NaN;
        [~,newj] = min(test);
        if max(newj)>0
            subplot(2,1,1), hold on
        %     line([PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),1) PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),1)]',[PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),3)-.5 PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
        %     plot(PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),1)+window(1),PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),3),'.k','MarkerSize',3)
    %         plot(dl(dl(:,2)==pfc(icell),1),dl(dl(:,2)==pfc(icell),3),'.','Color',newcol(ievent,:),'MarkerSize',3)
            plot(dl2(:,1),newj+toadd,'.','Color',newcol(ievent,:),'MarkerSize',5)


            xlim(plotwindow)        
    %         toadd = toadd+sum(touseR);
            toadd = toadd+max(newj);
            ylim([0 toadd])
            yl = get(gca,'ylim');
            plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            xlabel(['Time From Start of ' xlab ' (ms)'])
            ylabel(ylab)
            set(gca,'FontSize',18,'FontName','Helvetica')    
            title(['Rat ' thisdir(1) ', Day ' num2str(dday(id)) ', Cell ' num2str(pfc(icell))])
        end
        subplot(2,1,2), hold on
    %     patch([ind ind(end:-1:1)],[dat-sem rev(end:-1:1)],'black','FaceAlpha',.1,'EdgeAlpha',.1)
    %     plot(ind,sm(length(dat)*2:length(dat)*3-1)./binsize,'k','LineWidth',2)
        p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p.FaceColor= newcol(ievent,:);
        p.EdgeColor= newcol(ievent,:);
        p.FaceAlpha=.1;
        p.EdgeAlpha=0;
        plot(ind,dat,'Color',newcol(ievent,:),'LineWidth',2)
        xlim(plotwindow)
        if ievent==2
            yl = get(gca,'ylim');
            plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            if withmodpatch
                patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*yl(1)+((yl(2)-yl(1))*binsize)],'red','FaceAlpha',.3,'EdgeAlpha',0)
                patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*yl(1)+((yl(2)-yl(1))*binsize)],'black','FaceAlpha',.3,'EdgeAlpha',0)
            end
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            
            ttt = text(plotwindow(2)*.3,yl(1)+(range(yl)*.95),['p = ' num2str(round(p_SSD1(icell),2,'significant'))]);
            if p_SSD1(icell)<.05
                ttt.Color = 'r';
            end
            ylim(yl)
        end
    end
    xlabel(['Time From Start of ' xlab ' (ms)'])
    ylabel('Firing Rate (Hz)')
    set(gcf,'Position',[2271        -103         730        1099])
    set(gca,'FontSize',18,'FontName','Ariel')
    
    if ~isfolder([savefolder '\CandEventTriggered\'])
        mkdir([savefolder '\CandEventTriggered\'])
    end
    set(gcf,'renderer','Painters')
    if p_SSD1(icell)<.05
        labaddto = '_sig';
    elseif p_SSD1(icell)<.1
        labaddto = '_trend';
    else
        labaddto = '_ns';
    end
    helper_saveandclosefig([savefolder '\CandEventTriggered\Replays_' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) labaddto])
    
    toadd = 0;
    figure; hold on
    
    ievent=1;
%      touseR = singlebothjoint~=3;
    touseR = true(size(touseR22));
      
    dl = PFCcandspikes_list;
    dl(ismember(PFCcandspikes_list(:,3),find(~touseR)),:) = [];

    dat1 = squeeze(mean(PFCcandspikes_binned(icell,:,touseR),3))./binsize;
    sem1 = squeeze(std(PFCcandspikes_binned(icell,:,touseR)./binsize,[],3))./sqrt(sum(touseR));
    sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
    dat = conv(dat1',gaussFilter,'same')';
    sem = conv(sem1',gaussFilter,'same')';
    rev = dat+sem;
    fwd = dat-sem;     


%     rev = dat+sem;
%     sm = smoothts([dat dat dat dat],'b',smoothsize);

%         modind = [ind(ind>0 & ind<=.2)];
%         baseind = [ind(ind>=-.5 & ind<-.1)];

     dl2 = dl(dl(:,2)==pfc(icell),:);
%         [~,~,jj] = unique(dl2(:,3));
    jj = unique(dl2(:,3));
    test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
    test(test<0) = NaN;
    [~,newj] = min(test);

    subplot(2,1,1), hold on
%     line([PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),1) PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),1)]',[PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),3)-.5 PFCcandspikes_list(PFCcandspikes_list(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
%     plot(PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),1)+window(1),PFCcandspikes_list(PFCcandspikes_list(:,2)==pfc(icell),3),'.k','MarkerSize',3)
%         plot(dl(dl(:,2)==pfc(icell),1),dl(dl(:,2)==pfc(icell),3),'.','Color',newcol(ievent,:),'MarkerSize',3)
    plot(dl2(:,1),newj+toadd,'.','Color',newcol(ievent,:),'MarkerSize',5)
    figure; hold on; line([dl2(:,1) dl2(:,1)]',[newj'+toadd-1 newj'+toadd+1]','Color',newcol(ievent,:),'LineWidth',1)



    xlim(plotwindow)        
%         toadd = toadd+sum(touseR);
    toadd = toadd+max(newj);
    ylim([0 toadd])
    yl = get(gca,'ylim');
    plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
    set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
    xlabel(['Time From Start of ' xlab ' (ms)'])
    ylabel(ylab)
    set(gca,'FontSize',18,'FontName','Helvetica')    
    title(['Rat ' thisdir(1) ', Day ' num2str(dday(id)) ', Cell ' num2str(pfc(icell))])

    subplot(2,1,2), hold on
%     patch([ind ind(end:-1:1)],[dat-sem rev(end:-1:1)],'black','FaceAlpha',.1,'EdgeAlpha',.1)
%     plot(ind,sm(length(dat)*2:length(dat)*3-1)./binsize,'k','LineWidth',2)
    p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
    p.FaceColor= newcol(ievent,:);
    p.EdgeColor= newcol(ievent,:);
    p.FaceAlpha=.1;
    p.EdgeAlpha=0;
    plot(ind,dat,'Color',newcol(ievent,:),'LineWidth',2)
    xlim(plotwindow)
    yl = get(gca,'ylim');
    plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
    if withmodpatch
        patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*yl(1)+((yl(2)-yl(1))*binsize)],'red','FaceAlpha',.3,'EdgeAlpha',0)
        patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*yl(1)+((yl(2)-yl(1))*binsize)],'black','FaceAlpha',.3,'EdgeAlpha',0)
    end
    set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
     ttt = text(plotwindow(2)*.3,yl(1)+(range(yl)*.95),['p = ' num2str(round(p_SSD2(icell),2,'significant'))]);
            if p_SSD2(icell)<.05
                ttt.Color = 'r';
            end
    ylim(yl)
    xlabel(['Time From Start of ' xlab ' (ms)'])
    ylabel('Firing Rate (Hz)')
    set(gcf,'Position',[2271        -103         730        1099])
    set(gca,'FontSize',18,'FontName','Ariel')
    
    set(gcf,'renderer','Painters')
    if p_SSD2(icell)<.05
        labaddto = '_sig';
    elseif p_SSD2(icell)<.1
        labaddto = '_trend';
    else
        labaddto = '_ns';
    end
    helper_saveandclosefig([savefolder '\CandEventTriggered\EventOnly_' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) labaddto])

%     helper_saveandclosefig(['D:\XY_matdata\Figures\ForPaper\CandEventTriggered\' label '_' dirname '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell))])
end