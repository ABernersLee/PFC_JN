function plot_PFC_FwdRevReplay_triggered(thisdir,label,withmodpatch,smoothsize)
%%
load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_corr'],[label '_replay_dirbias'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'other_cells','dirname')


eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Rcorr = ' label '_replay_corr;'])
eval(['dirbias = ' label '_replay_dirbias;'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])

pfc = other_cells; clear other_cells

%only single replays
PFCreplayspikes_binned(:,:,singlebothjoint~=1) = [];

% pos Rcorr is out, neg Rcorr is in ??? DOUBLE CHECK
%-1 bias is out, 1 is in
fwd = (Rcorr>0 & dirbias<0) | (Rcorr<0 & dirbias>0);
rev = (Rcorr>0 & dirbias>0) | (Rcorr<0 & dirbias<0);
include = NaN(size(PFCreplayspikes_binned,3),1);
include(fwd) = 1;
include(rev) = 2;
include1 = include(singlebothjoint==1);

PFCreplayspikes_binned(:,:,isnan(include1)) = [];
PFCreplayspikes_list(ismember(PFCreplayspikes_list(:,3),find(singlebothjoint~=1 | isnan(include))),:) = [];
include1(isnan(include1)) = [];
fwdrev = include1;
numcat   = 2;

binsize = .002;
window = [-.5 .5];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

modind = [ind(ind>0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<-.1)];
newcol = [1 0 0; 0 0 0];
        
for icell = 1:length(pfc)
    figure; hold on
    a = subplot(2,1,1); hold on
    b = subplot(2,1,2); hold on
    toadd = 0;
    for idir = 1:numcat
        db = PFCreplayspikes_binned(icell,:,fwdrev==idir);
        dl = PFCreplayspikes_list;
        dl(ismember(PFCreplayspikes_list(:,3),find(include~=idir)),:) = [];
        
        dat = squeeze(mean(db(1,:,:),3));
        sem = squeeze(std(db(1,:,:),[],3))./sqrt(size(db,3));
        rev = dat+sem;
        sm = smoothts([dat dat dat dat],'b',smoothsize);
        
        
        [~,~,jj] = unique(dl(dl(:,2)==pfc(icell),3));
        subplot(a)
    % %     line([dl(dl(:,2)==pfc(icell),1) dl(dl(:,2)==p(icell),1)]',[dl(dl(:,2)==p(icell),3)-.5 dl(dl(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
        plot(dl(dl(:,2)==pfc(icell),1)+window(1),jj+toadd,'.','Color',newcol(idir,:),'MarkerSize',5)
        
        toadd = toadd+max(jj);
        if idir == numcat
            xlim(window)
            ylim([0 toadd])    
            yl = get(gca,'ylim');
            plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            xlabel(['Time From Start of Replay (ms)'])
            ylabel('Replay #')
            set(gca,'FontSize',18,'FontName','Helvetica')
        end
        
        subplot(b)
        p = patch([ind ind(end:-1:1)],[dat-sem rev(end:-1:1)],'black');
        p.FaceColor= newcol(idir,:);
        p.FaceAlpha=.1;
        p.EdgeAlpha=.1;
        plot(ind,sm(length(dat)*2:length(dat)*3-1)./binsize,'Color',newcol(idir,:),'LineWidth',2)

        if idir == numcat
            yl = get(gca,'ylim');
            yl(1) = max([yl(1) 0]);
            plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            if withmodpatch                
                patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*(yl(2)-yl(1))*.02],'red','FaceAlpha',.3,'EdgeAlpha',0)
                patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*(yl(2)-yl(1))*.02],'black','FaceAlpha',.3,'EdgeAlpha',0)
            end
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            ylim(yl)    
            xlabel(['Time From Start of Replay (ms)'])
            ylabel('Firing Rate (Hz)')    
            set(gca,'FontSize',18,'FontName','Helvetica')
        end
    end
    set(gcf,'Position',[2271        -103         730        1099])
    helper_saveandclosefig(['D:\XY_matdata\Figures\ForPaper\FwdRevReplayTriggered\' label '_' dirname '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell))])
end