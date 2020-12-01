function plot_PFC_ArmTheta_triggered(thisdir,withmodpatch,smoothsize,cutoff,st)
label = num2str(cutoff);
load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','other_cells')
if st 
    load(thisdir,'PFCthetaspikes_binned_st','PFCthetaspikes_list_st');
    PFCthetaspikes_binned = PFCthetaspikes_binned_st;
    PFCthetaspikes_list = PFCthetaspikes_list_st;
    stlab = '_fromTHstart';
else
    load(thisdir,'PFCthetaspikes_binned','PFCthetaspikes_list');
    stlab = ''; %from Theta end
end
pfc = other_cells; clear other_cells
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
clear times_armon_thetaof_headingarm_lap_thetahalf_all

sigma = smoothsize;
thetaarm = Th(:,4);
thetaon = Th(:,3);
thetaarm(thetaon==thetaarm | isnan(thetaon) | isnan(thetaarm) | Th(:,end)<cutoff | isnan(Th(:,end))) = NaN;
        
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

ind = [window(1):binsize:window(2)];

modind = [ind(ind>=0 & ind<=.2)];
baseind = [ind(ind>=-.5 & ind<=-.1)];
% newcol = [0 1 .9; 1 0 0; .8 .2 1];

newcol = [75 0 130;34 139 34;255 130 0]/255;
% plotwindow = [-.6 .6];        
plotwindow = [-.05 .1];        

for icell = 1:length(pfc)
    
    figure; hold on
    a = subplot(2,1,1); hold on
    b = subplot(2,1,2); hold on
    toadd = 0;
    for iarm = 1:3
        
        db = PFCthetaspikes_binned(icell,:,thetaarm==iarm);
        dl = PFCthetaspikes_list;
        dl(ismember(PFCthetaspikes_list(:,3),find(thetaarm~=iarm)),:) = [];
        if isempty(dl); continue; end
        
        dat1 = squeeze(mean(db(1,:,:),3))./binsize;
        sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        
        
        sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
        gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
        dat = conv(dat1',gaussFilter,'same')';
        sem = conv(sem1',gaussFilter,'same')';
        rev = dat+sem;
        fwd = dat-sem;        
        
        dl2 = dl(dl(:,2)==pfc(icell),:);
%         [~,~,jj] = unique(dl2(:,3));
        jj = unique(dl2(:,3));
        test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
        test(test<0) = NaN;
        [~,newj] = min(test);
        subplot(a)
    % %     line([dl(dl(:,2)==pfc(icell),1) dl(dl(:,2)==p(icell),1)]',[dl(dl(:,2)==p(icell),3)-.5 dl(dl(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
%         plot(dl(dl(:,2)==pfc(icell),1),jj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
        plot(dl2(:,1),newj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',10)
        if isempty(newj); continue; end
        toadd = toadd+max(newj);
        if iarm == 3
            xlim(plotwindow)
            if toadd~=0
                ylim([0 toadd])    
            end
            yl = get(gca,'ylim');
            plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            xlabel(['Time From Start of Theta (ms)'])
            ylabel('Theta #')
            set(gca,'FontSize',18,'FontName','Helvetica')
        end
        
        subplot(b), hold on
        p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
        p.FaceColor= newcol(iarm,:);
        p.EdgeColor= newcol(iarm,:);
        p.FaceAlpha=.1;
        p.EdgeAlpha=0;
        plot(ind,dat,'Color',newcol(iarm,:),'LineWidth',2)

        if iarm == 3
            set(gcf,'Position',[2271        -103         730        1099])
            xlim(plotwindow)
            
            
            set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
            
            xlabel(['Time From Start of Theta (ms)'])
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
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\ArmThetaTriggered\')
        mkdir('E:\XY_matdata\Figures\ForPaper\ArmThetaTriggered\')
    end
    suptitle([ label ' ' thisdir(1) thisdir(3:end-4) ' Cell' num2str(pfc(icell)) stlab])
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmThetaTriggered\' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) stlab])    
    
    if 0
    
        for iarmon = 1:3
        oths = setdiff(1:3,iarmon);
        figure; hold on
        a = subplot(2,1,1); hold on
        b = subplot(2,1,2); hold on
        toadd = 0;
        for ioth = 1:2
            iarm = oths(ioth);
            db = PFCthetaspikes_binned(icell,:,thetaarm==iarm & thetaon==iarmon);
            dl = PFCthetaspikes_list;
            dl(ismember(PFCthetaspikes_list(:,3),find(thetaarm~=iarm | thetaon~=iarmon)),:) = [];
            if isempty(dl); continue; end
            
            dat1 = squeeze(mean(db(1,:,:),3))./binsize;
            sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        

            sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
            gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
            dat = conv(dat1',gaussFilter,'same')';
            sem = conv(sem1',gaussFilter,'same')';
            rev = dat+sem;
            fwd = dat-sem;        

            dl2 = dl(dl(:,2)==pfc(icell),:);
    %         [~,~,jj] = unique(dl2(:,3));
            jj = unique(dl2(:,3));
            test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
            test(test<0) = NaN;
            [~,newj] = min(test);
            if isempty(newj); continue; end
            subplot(a)
        % %     line([dl(dl(:,2)==pfc(icell),1) dl(dl(:,2)==p(icell),1)]',[dl(dl(:,2)==p(icell),3)-.5 dl(dl(:,2)==p(icell),3)+.5]','Color','k','LineWidth',3)
    %         plot(dl(dl(:,2)==pfc(icell),1),jj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
            plot(dl2(:,1),newj+toadd,'.','Color',newcol(iarm,:),'MarkerSize',5)
            toadd = toadd+max(newj);
            if ioth == 2 & iarmon==3
                xlim(plotwindow)
                if toadd~=0
                    ylim([0 toadd])    
                end
%                 yl = get(gca,'ylim');
                plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
                xlabel(['Time From Start of Theta (ms)'])
                ylabel('Theta #')
                set(gca,'FontSize',18,'FontName','Helvetica')
            end

            subplot(b), hold on
            p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
            p.FaceColor= newcol(iarm,:);
            p.EdgeColor= newcol(iarm,:);
            p.FaceAlpha=.1;
            p.EdgeAlpha=0;
            plot(ind,dat,'Color',newcol(iarm,:),'LineWidth',2)

            if ioth == 2 & iarmon==3
                set(gcf,'Position',[2271        -103         730        1099])
                xlim(plotwindow)


                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])

                xlabel(['Time From Start of Theta (ms)'])
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
        helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmThetaTriggered\' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) '_armon' num2str(iarmon) stlab])    
        end
        
    
   
        figure; hold on    
        for iarm = 1:3
            oths = setdiff(1:3,iarm);
               clear dat2 sem2
            for ioth = 1:2
                iarmon = oths(ioth);
                db = PFCthetaspikes_binned(icell,:,thetaarm==iarm & thetaon==iarmon);
                dl = PFCthetaspikes_list;
                dl(ismember(PFCthetaspikes_list(:,3),find(thetaarm~=iarm | thetaon~=iarmon)),:) = [];
                if isempty(dl); continue; end

                dat1 = squeeze(mean(db(1,:,:),3))./binsize;
                sem1 = squeeze(std(db(1,:,:)./binsize,[],3))./sqrt(size(db,3));        

                sz = length(dat1); x = linspace(-sz / 2, sz / 2, sz);
                gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2)); gaussFilter = gaussFilter / sum (gaussFilter);
                dat2(:,ioth) = conv(dat1',gaussFilter,'same')';
                sem2(:,ioth) = conv(sem1',gaussFilter,'same')';
            end
            dat = nanmean(dat2,2)';
            sem = nanmean(sem2,2)';
            rev = dat+sem;
            fwd = dat-sem;        

            dl2 = dl(dl(:,2)==pfc(icell),:);
    %         [~,~,jj] = unique(dl2(:,3));
            jj = unique(dl2(:,3));
            test = repmat(dl2(:,3)', [size(jj,1) 1])-repmat(jj,[1 size(dl2,1)]);
            test(test<0) = NaN;
            [~,newj] = min(test);
            if isempty(newj); continue; end


            p = patch([ind ind(end:-1:1)],[fwd rev(end:-1:1)],'black');
            p.FaceColor= newcol(iarm,:);
            p.EdgeColor= newcol(iarm,:);
            p.FaceAlpha=.1;
            p.EdgeAlpha=0;
            plot(ind,dat,'Color',newcol(iarm,:),'LineWidth',2)

            if iarm==3
                xlim(plotwindow)
                set(gca,'xtick',[-.5 -.25 0 .25 .5],'xticklabel',[-500 -250 0 250 500])
                xlabel(['Time From Start of Theta (ms)'])
                ylabel('Firing Rate (Hz)')    
                set(gca,'FontSize',18,'FontName','Helvetica')
                yl = get(gca,'ylim');
                plot([0 0], yl,':','Color',[.5 .5 .5],'LineWidth',3)
            end
        end                
        %         helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmThetaTriggered\' label '_' thisdir(1:end-4) '_smooth' num2str(smoothsize) '_' num2str(withmodpatch) 'modpatch_Cell' num2str(pfc(icell)) '_armonall'  stlab])    
    end
    
    
end

