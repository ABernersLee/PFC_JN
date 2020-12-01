function plot_lapbylap_PFC_all(label,laptype)   

%across laps get for each PFC cell:
    %spatial information 
    
    %FR
    %theta locking mrv
    %theta phase preference
    %theta preceission correlation x dir
    %theta preceission slope x dir
        %and across each arm and heading arm
        
    %replay event modulation
    %arm-replay modulation SSD
        %modu across each arm
    
% lapdata = Cell Num X lapnum X 10
% laparmdata = Cell Num X lapnum X 8 X 3
cd('E:\XY_matdata\AllDays\v1\')
d2 = dir('*.mat');
lapdata = []; laparmdata = []; datalapP = []; datalapPZ = []; datahalf = []; parm = [];
for id = [1:size(d2,1)]
    thisdir = d2(id).name;
    load(thisdir,[label 'laptype' num2str(laptype) '_lapdata'],[label 'laptype' num2str(laptype) '_laparmdata'],[label '_pSSDarm'])
    
    eval(['lapdata1 = ' label 'laptype' num2str(laptype) '_lapdata;'])
    eval(['pSSD = ' label '_pSSDarm;'])
    eval(['laparmdata1 = ' label 'laptype' num2str(laptype) '_laparmdata;'])
    clear([label 'laptype' num2str(laptype) '_lapdata'],[label 'laptype' num2str(laptype) '_laparmdata'])
    
    for itype = 1:size(lapdata1,2)
       tmp = squeeze(lapdata1(:,itype,:));
       tmp(:,:,2) = repmat(((itype)-1)./(size(lapdata1,2)-1),[size(tmp,1) size(tmp,2)]);
       datalapP = cat(1,datalapP,tmp);        
       tmpz = tmp; tmpz(:,:,1) = nanzscore(tmp(:,:,1),[],2);
       datalapPZ = cat(1,datalapPZ,tmpz);
    end
    n = floor(size(lapdata1,2)/2);
    half = NaN(size(lapdata1,1),2,size(lapdata1,3));
    for icell = 1:size(lapdata1,1)
        half(icell,1,:) = nanmean(lapdata1(icell,1:n,:));
        half(icell,2,:) = nanmean(lapdata1(icell,end+1-n:end,:));
    end
    datahalf = cat(1,datahalf,half);
%     if size(lapdata,2)>size(lapdata1,2) && id>1
%         lapdata1 = cat(2,lapdata1,NaN(size(lapdata1,1),size(lapdata,2)-size(lapdata1,2),size(lapdata1,3)));
%         laparmdata1 = cat(2,laparmdata1,NaN(size(laparmdata1,1),size(laparmdata,2)-size(laparmdata1,2),size(laparmdata1,3),size(laparmdata1,4)));
%     elseif size(lapdata1,2)>size(lapdata,2) && id>1
%         lapdata = cat(2,lapdata,NaN(size(lapdata,1),size(lapdata1,2)-size(lapdata,2),size(lapdata,3)));
%         laparmdata = cat(2,laparmdata,NaN(size(laparmdata,1),size(laparmdata1,2)-size(laparmdata,2),size(laparmdata,3),size(laparmdata,4)));
%     end
        if size(lapdata,2)>size(lapdata1,2) && id>1
            lapdata = lapdata(:,1:size(lapdata1,2),:);
            laparmdata = laparmdata(:,1:size(laparmdata1,2),:,:);
        elseif size(lapdata1,2)>size(lapdata,2) && id>1
            lapdata1 = lapdata1(:,1:size(lapdata,2),:);
            laparmdata1 = laparmdata1(:,1:size(laparmdata,2),:,:);
        end

    lapdata = cat(1,lapdata,lapdata1);
    laparmdata = cat(1,laparmdata,laparmdata1);
    parm = cat(1,parm,pSSD);
    
end

laplab = {'Spatial Information';'Firing Rate';'HP Theta MRV'; ...
    'HP Theta Precession Corr In'; 'HP Theta Precession Corr Out';...
    'HP Theta Precession Slope In'; 'HP Theta Precession Slope Out'; ...
    [label ' Modulation'];['Arm-Specific ' label ' Normalized SSD'];['Arm-Specific ' label ' Raw Modu SSD']};
    
% armlaplab = {'Firing Rate';'HP Theta MRV';'HP Theta Phase Pref'; ...
%     'HP Theta Precession Corr In'; 'HP Theta Precession Corr Out';...
%     'HP Theta Precession Slope In'; 'HP Theta Precession Slope Out'; ['Arm-Specific ' label ' Modulation']};

laplab2 = {'laps_coverspace','laps_twoarms','laps_singlepass'};

lapdata(:,:,[4 10]) = [];
datahalf(:,:,[4 10]) = [];
if 0
    
for iz = 1:2
    figure; hold on
    spn = ceil(sqrt(size(lapdata,3)));
    for icol = 1:size(lapdata,3)
        subplot(spn,ceil(size(lapdata,3)/spn),icol), hold on
        if iz == 2
            dat = nanzscore(lapdata(:,:,icol),[],2);
            if icol==9
                dat = nanzscore(abs(lapdata(:,:,icol)),[],2);
            end
        elseif iz == 1
            dat = lapdata(:,:,icol);
            if icol == 9
                dat = abs(lapdata(:,:,icol));
            end
        end
        
        pl = nanmean(dat,1);
        sem = nanstd(dat)./sqrt(sum(~isnan(dat),1));
        rev = pl+sem;
        patch([1:size(dat,2) size(dat,2):-1:1],[pl-sem rev(end:-1:1)],'black','FaceAlpha',.3,'EdgeAlpha',.3)
        plot(1:size(dat,2),pl,'k-','LineWidth',3)
        xlabel('Laps')
        ylabel(laplab{icol})
        [r,p] = corr([1:size(dat,2)]',pl','rows','complete');
        ll = legend(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))],'Location','Best');
        if p<.05
            ll.TextColor = 'r';
        end
        xlim([0 size(lapdata,2)+1])
    end
    subplot(spn,ceil(size(lapdata,3)/spn),icol+1)
    text(0,.4,laplab2{laptype})
    text(0,.8,label)
    if iz ==1
        text(0,0,['Raw'])
    elseif iz == 2
        text(0,0,['Z-scored'])
    end
    axis off
    set(gcf,'Position',[  2006          -9        1103         914])
    helper_saveandclosefig(['D:\XY_matdata\Figures\ForPaper\OverLaps\' label '_' laplab2{laptype} '_z' num2str(iz) '_cutoff'])
end




for iz = 1:2
    figure; hold on
    spn = ceil(sqrt(size(lapdata,3)));
    for icol = 1:size(lapdata,3)
        subplot(spn,ceil(size(lapdata,3)/spn),icol), hold on

        if iz == 1
            dat1 = squeeze(datalapP(:,icol,1));
        else
            dat1 = squeeze(datalapPZ(:,icol,1));
        end

        if icol == 9
            dat1 = abs(dat1);
        end
        dat2 = squeeze(datalapP(:,icol,2));
        plot(dat2,dat1,'*k','MarkerSize',10);
        xlabel('Percent of Epoch (each lap)')
        ylabel(laplab{icol})
        [r,p] = corr(dat2,dat1,'rows','complete');
        ll = legend(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))],'Location','Best');
        if p<.05
            ll.TextColor = 'r';
        end
        xlim([0 1])
    end
    subplot(spn,ceil(size(lapdata,3)/spn),icol+1)
    text(0,.4,laplab2{laptype})
    text(0,.8,label)
    if iz ==1
        text(0,0,['Raw'])
    elseif iz == 2
        text(0,0,['Z-scored'])
    end
    axis off
    set(gcf,'Position',[  2006          -9        1103         914])
    helper_saveandclosefig(['D:\XY_matdata\Figures\ForPaper\OverLaps\percent' label '_' laplab2{laptype} '_z' num2str(iz)])
end

end


for iz = 1 %:2
    figure; hold on
    spn = ceil(sqrt(size(lapdata,3)));
    
    for icol = 1:size(lapdata,3)
        subplot(spn,ceil(size(lapdata,3)/spn),icol), hold on

         if iz == 2
%             dat = nanzscore(lapdata(:,:,icol),[],2);
%             if icol==9
%                 dat = nanzscore(abs(lapdata(:,:,icol)),[],2);
%             end
        elseif iz == 1
            dat = datahalf(:,:,icol);
%             if icol == 9
%                 dat = abs(datahalf(:,:,icol));
            if icol>=size(lapdata,3)-2
                 dat = datahalf(parm<.05,:,icol);
            end
         end
         
        
        
        mdat1 = nanmean(dat(:,1));
        sem1 = nanstd(dat(:,1))./sqrt(sum(~isnan(dat(:,1))));
        errorbar(1,mdat1,sem1,'k-','LineWidth',2)
        
        mdat2 = nanmean(dat(:,2));
        sem2 = nanstd(dat(:,2))./sqrt(sum(~isnan(dat(:,2))));
        errorbar(2,mdat2,sem2,'k-','LineWidth',2)
        
%         h = boxplot2(dat1(i==1),1);
%         make_box_plot2_color(h,[0 0 0])
%         h = boxplot2(dat1(i==2),2);
%         make_box_plot2_color(h,[0 0 0])
        
%         xlabel('Percent of Epoch (each lap)')
        title(laplab{icol})
%         [r,p] = corr(c',meansem(:,1),'rows','complete');
%         p = ranksum(dat1(i==1),dat1(i==2));
        p = signrank(dat(:,1),dat(:,2));
        
        yl = get(gca,'ylim');
        if yl(2)<0
            ll = text(1.3,yl(2)*1.05,['p = ' num2str(round(p,2,'significant'))]);
        elseif yl(2)>=0
            ll = text(1.3,yl(2)*.95,['p = ' num2str(round(p,2,'significant'))]);
        end
        if p<.05
            ll.Color = 'r';
        end
        ll.FontSize = 12;
        xlim([.8 2.2])
        set(gca,'xtick',1:2,'xticklabel',{'First Half','Second Half'})
        set(gca,'ylim',yl)
        set(gca,'FontSize',12)
    end
    
    subplot(spn,ceil(size(lapdata,3)/spn),icol+1), hold on
     dat = cat(1,datahalf(:,:,5),datahalf(:,:,6));
        mdat1 = nanmean(dat(:,1));
        sem1 = nanstd(dat(:,1))./sqrt(sum(~isnan(dat(:,1))));
        errorbar(1,mdat1,sem1,'k-','LineWidth',2)
        
        mdat2 = nanmean(dat(:,2));
        sem2 = nanstd(dat(:,2))./sqrt(sum(~isnan(dat(:,2))));
        errorbar(2,mdat2,sem2,'k-','LineWidth',2)
        
%         xlabel('Percent of Epoch (each lap)')
        title('Phase Precession r (concat dir)')
        p = signrank(dat(:,1),dat(:,2));
        
        yl = get(gca,'ylim');        
        if yl(2)<0
            ll = text(1.3,yl(2)*1.05,['p = ' num2str(round(p,2,'significant'))]);
        elseif yl(2)>=0
            ll = text(1.3,yl(2)*.95,['p = ' num2str(round(p,2,'significant'))]);        
        end
        if p<.05
            ll.Color = 'r';
        end
        ll.FontSize = 12;
        xlim([.8 2.2])
        set(gca,'ylim',yl)
        set(gca,'xtick',1:2,'xticklabel',{'First Half','Second Half'})
           set(gca,'FontSize',12)
    
    
        subplot(spn,ceil(size(lapdata,3)/spn),icol+2), hold on
        dat = cat(1,datahalf(:,:,7),datahalf(:,:,8));
        mdat1 = nanmean(dat(:,1));
        sem1 = nanstd(dat(:,1))./sqrt(sum(~isnan(dat(:,1))));
        errorbar(1,mdat1,sem1,'k-','LineWidth',2)
        
        mdat2 = nanmean(dat(:,2));
        sem2 = nanstd(dat(:,2))./sqrt(sum(~isnan(dat(:,2))));
        errorbar(2,mdat2,sem2,'k-','LineWidth',2)
        
%         xlabel('Percent of Epoch (each lap)')
        title('Phase Precession slope (concat dir)')
        p = signrank(dat(:,1),dat(:,2));
        
        yl = get(gca,'ylim');
        if yl(2)<0
            ll = text(1.3,yl(2)*1.05,['p = ' num2str(round(p,2,'significant'))]);
        elseif yl(2)>=0
            ll = text(1.3,yl(2)*.95,['p = ' num2str(round(p,2,'significant'))]);
        end
        if p<.05
            ll.Color = 'r';
        end
        ll.FontSize = 12;
        xlim([.8 2.2])
        set(gca,'ylim',yl)

        set(gca,'xtick',1:2,'xticklabel',{'First Half','Second Half'})
        set(gca,'FontSize',12)
    
%     text(0,.4,laplab2{laptype})
%     text(0,.8,label)
%     if iz ==1
%         text(0,0,['Raw'])
%     elseif iz == 2
%         text(0,0,['Z-scored'])
%     end
%     axis off
    
    set(gcf,'Position',[  2006          -9        1103         914])
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\OverLaps\')        
        mkdir('E:\XY_matdata\Figures\ForPaper\OverLaps\')
    end
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\OverLaps\' laplab2{laptype} '_z' num2str(iz) '_' label])
end
end