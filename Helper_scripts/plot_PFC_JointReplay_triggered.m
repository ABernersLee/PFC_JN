function plot_PFC_JointReplay_triggered(thisdir,label,withmodpatch,smoothsize,timetrig,beta,armdiff1,armdiff2)

tosave = false;
[~,PFCreplayspikes_binned] = make_PFC_jointreplayeventtriggeredmat(thisdir,label,timetrig,tosave);

load(thisdir,[label  '_replay_jointreplayarm'],[label '_Cand_sig_modu_include'],'other_cells','dirname')
eval(['JointArm = ' label  '_replay_jointreplayarm;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
clear([label '_jointreplayarm'],[label '_Cand_sig_modu_include'])
JointArm(sum(~isnan(JointArm),2)==0,:) = [];
replayarm = JointArm(:,1);
otharms = JointArm(:,2);

pfc = other_cells(Cand_sig_modu_include(:,3)==1); clear other_cells;
PFCreplayspikes_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];
sigma = 1.5;
        
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

ind = [window(1):binsize:window(2)];

modind = [ind(ind>=0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<=-.1)];
newcol = [0 1 .9; 1 0 0; .8 .2 1;.5 .5 .5; 1 0 0; .8 .2 1];

plotwindow = [-.6 .6];        

for icell = 1:length(pfc)
    figure; hold on
%     a = subplot(2,1,1); hold on
%     b = subplot(2,1,2); hold on
%     toadd = 0;
    for iarm = 1:3
        otharm = setdiff(1:3,iarm);
        for ioth = 1:2
            
            db = PFCreplayspikes_binned(icell,:,replayarm==iarm & otharms==otharm(ioth));
%             dl = PFCreplayspikes_list;
%             dl(ismember(PFCreplayspikes_list(:,3),find(replayarm~=iarm | otharms~=otharm(ioth))),:) = [];

            dat1 = squeeze(mean(db(1,:,:),3))./binsize;
            sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        

            sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
            dat = conv(dat1',gaussFilter,'same')';
            sem = conv(sem1',gaussFilter,'same')';
            rev = dat+sem;
            fwd = dat-sem;        


%             dl2 = dl(dl(:,2)==pfc(icell),:);
%     %         [~,~,jj] = unique(dl2(:,3));
% %             jj = unique(dl2(:,3));
%             test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
%             test(test<0) = NaN;
%             [~,newj] = min(test);
%             subplot(a)
%             plot(dl2(:,1),newj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
%             toadd = toadd+max(newj);
%             if iarm == 3 && ioth==2
%                 xlim(plotwindow)
%                 ylim([0 toadd])    
%                 yl = get(gca,'ylim');
%                 plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
%                 set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
%                 xlabel(['Time From Start of Replay (ms)'])
%                 ylabel('Replay #')
%                 set(gca,'FontSize',18,'FontName','Helvetica')
%             end

%             subplot(b), hold on
            subplot(3,2,(iarm-1)*2+ioth); hold on
            p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
            p.FaceColor= newcol(iarm,:);
            p.EdgeColor= newcol(iarm,:);
            p.FaceAlpha=.1;
            p.EdgeAlpha=0;
            plot(ind,dat,'Color',newcol(iarm,:),'LineWidth',2)

            if iarm == 3 && ioth==2
                set(gcf,'Position',[ 2271        -103        1083        1099])
            end
                xlim(plotwindow)


                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
                    
            if iarm == 3
                if timetrig==1
                    xlabel(['Time From Start of Replay (ms)'])
                elseif timetrig==2
                    xlabel(['Time From Middle of Replay (ms)'])
                end
            end
            if ioth == 1
                ylabel(['c: ' num2str(round(beta(icell,iarm),3,'significant')) ', d: ' num2str(round(armdiff1(icell,iarm),3,'significant'))...
                    ', a: ' num2str(round(armdiff2(icell,iarm),3,'significant'))])    
            end
            set(gca,'FontSize',18,'FontName','Helvetica')
            title(['JR: ' num2str(iarm) '-> ' num2str(otharm(ioth))])

                yl = get(gca,'ylim');
                plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
    %             yl(1) = max([yl(1) 0]);

                if withmodpatch                
                    patch([modind modind(end:-1:1)],[ones(size(modind))*yl(1) ones(size(modind))*yl(1)+((yl(2)-yl(1))*.02)],'red','FaceAlpha',.3,'EdgeAlpha',0)
                    patch([baseind baseind(end:-1:1)],[ones(size(baseind))*yl(1) ones(size(baseind))*yl(1)+((yl(2)-yl(1))*.02)],'black','FaceAlpha',.3,'EdgeAlpha',0)
                end
%                 ttt = text(plotwindow(2)*.3,yl(1)+(range(yl)*.95),['p = ' num2str(round(p_SSD(icell),2,'significant'))]);
%                 if p_SSD(icell)<.05
%                     ttt.Color = 'r';
%                 end

                ylim(yl) 
            end
    end
        
    suptitle([dirname ', Cell ' num2str(pfc(icell)) ', mean coeff: ' num2str(nanmean(beta(icell,:)))])
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\JointArmReplayTriggered\')
        mkdir('E:\XY_matdata\Figures\ForPaper\JointArmReplayTriggered\')
    end
    

    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\JointArmReplayTriggered\' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) '_timetrig' num2str(timetrig)])
    
    
    
end