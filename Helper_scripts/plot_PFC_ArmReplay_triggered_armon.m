function plot_PFC_ArmReplay_triggered_armon(thisdir,label,withmodpatch,smoothsize)

load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_replayarm'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'other_cells',...
    [label '_Cand_sig_modu_include'],[label '_replay_stnd'],[label '_pSSDarm'],[label  '_replay_shuffle_p'],'pos','armpos')

eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['p_SSD = ' label '_pSSDarm;'])
eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])
pfc = other_cells; %(Cand_sig_modu_include(:,3)==1); clear other_cells;
% p_SSD(Cand_sig_modu_include(:,3)~=1) = [];
% PFCreplayspikes_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];
%only single replays
% PFCreplayspikes_binned(:,:,singlebothjoint==3) = [];
% PFCreplayspikes_list(ismember(PFCreplayspikes_list(:,3),find(singlebothjoint==3)),:) = [];
% replayarm1 = replayarm;
% replayarm = replayarm1(singlebothjoint~=3);
replayarm(singlebothjoint==3) = NaN;
sigma = 1.5;


        
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

ind = [window(1):binsize:window(2)];

modind = [ind(ind>=0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<=-.1)];
% newcol = [0 1 .9; 1 0 0; .8 .2 1];

newcol = [75 0 130;255 130 0;34 139 34]/255;
plotwindow = [-.6 .6];        

[~,~,i] = histcounts(Event,pos(:,1));
% Event(i==0) = [];
% replayarm(i==0) = [];
armon = armpos(i);

for icell = 1:length(pfc)
    for iarmon = 1:3
        oths = setdiff(1:3,iarmon);
        figure; hold on
        a = subplot(2,1,1); hold on
        b = subplot(2,1,2); hold on
        toadd = 0;
        for ioth = 1:2
            iarm = oths(ioth);
            db = PFCreplayspikes_binned(icell,:,replayarm==iarm & armon==iarmon);
            dl = PFCreplayspikes_list;
            dl(ismember(PFCreplayspikes_list(:,3),find(replayarm~=iarm | armon~=iarmon)),:) = [];
            if isempty(dl); continue; end
            dat1 = squeeze(mean(db(1,:,:),3))./binsize;
            sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        

            sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
            dat = conv(dat1',gaussFilter,'same')';
            sem = conv(sem1',gaussFilter,'same')';
            rev = dat+sem;
            fwd = dat-sem;        

    %         sm = smoothts([dat1 dat1 dat1 dat1],'b',smoothsize);        
    %         dat = sm(length(dat1)*2:length(dat1)*3-1);
    %         rev1 = dat1+sem1;
    %         fwd1 = dat1-sem1;        
    %         rev2 = smoothts([rev1 rev1 rev1 rev1],'b',smoothsize); 
    %         rev = rev2(length(dat1)*2:length(dat1)*3-1);
    %         fwd2 = smoothts([fwd1 fwd1 fwd1 fwd1],'b',smoothsize); 
    %         fwd = fwd2(length(dat1)*2:length(dat1)*3-1);

            dl2 = dl(dl(:,2)==pfc(icell),:);
    %         [~,~,jj] = unique(dl2(:,3));
            jj = unique(dl2(:,3));
            test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
            test(test<0) = NaN;
            [~,newj] = min(test);
            if isempty(test); continue; end
            subplot(a)
        % %     line([dl(dl(:,2)==pfc(icell),1) dl(dl(:,2)==p(icell),1)]',[dl(dl(:,2)==p(icell),3)-.5 dl(dl(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
    %         plot(dl(dl(:,2)==pfc(icell),1),jj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
            plot(dl2(:,1),newj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
            toadd = toadd+max(newj);
            if ioth == 2
                xlim(plotwindow)
                if toadd~=0
                    ylim([0 toadd])   
                end                
                yl = get(gca,'ylim');
                plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
                xlabel(['Time From Start of Replay (ms)'])
                ylabel('Replay #')
                set(gca,'FontSize',18,'FontName','Helvetica')
            end

            subplot(b), hold on
            p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
            p.FaceColor= newcol(iarm,:);
            p.EdgeColor= newcol(iarm,:);
            p.FaceAlpha=.1;
            p.EdgeAlpha=0;
            plot(ind,dat,'Color',newcol(iarm,:),'LineWidth',2)

            if ioth == 2
                set(gcf,'Position',[2271        -103         730        1099])
                xlim(plotwindow)


                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])

                xlabel(['Time From Start of Replay (ms)'])
                ylabel('Firing Rate (Hz)')    
                set(gca,'FontSize',18,'FontName','Helvetica')

                yl = get(gca,'ylim');
                plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
    %             yl(1) = max([yl(1) 0]);

                if withmodpatch                
                    patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*yl(1)+((yl(2)-yl(1))*.02)],'red','FaceAlpha',.3,'EdgeAlpha',0)
                    patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*yl(1)+((yl(2)-yl(1))*.02)],'black','FaceAlpha',.3,'EdgeAlpha',0)
                end
              

                ylim(yl) 
            end
        end
        if ~isfolder('E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\')
            mkdir('E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\')
        end
        helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmReplayTriggered\' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) '_armon' num2str(iarmon)])
    end
end