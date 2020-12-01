function plot3PFC_arm_modulation(label,igroup,savefolder)
tlab = {'Modu';'SSD'};
for itype = 1:2
    cd('F:\XY_matdata\AllDays')
    d2 = dir('*.mat');
    modu = []; p_SSD = [];
    for id = 1:size(d2,1)
        thisdir = d2(id).name;
        load(thisdir,[label '_pSSDarm'],[label '_moduarm'],[label '_SSDarm'],[label '_moduarm2'])
        eval(['p_SSD2 = ' label '_pSSDarm;'])
        if itype ==1
            eval(['modu2 = ' label '_moduarm'';'])
        else
            eval(['modu2 = ' label '_moduarm2'';'])
        end
        load(thisdir,'other_cells_touse')        
        modu = cat(1,modu,modu2(other_cells_touse(:,igroup),:));
        p_SSD = cat(1,p_SSD,p_SSD2(other_cells_touse(:,igroup),:));
    end
        modu2 = modu';
    
    % p_SSD = p_SSD2;
    ind = sum(isnan(modu2))==0 & ~isnan(p_SSD)';
    modu2(:,~ind) = [];
    toplot = NaN(3,size(modu2,2));
    for iarm = 1:3
        toplot(iarm,:) = [modu2(iarm,p_SSD(ind)>=.05)';modu2(iarm,p_SSD(ind)<.05)'];
    end

    jnk = toplot>0;
    allcombs = [1 1 1; 1 1 0; 1 0 1; 0 1 1; 0 1 0; 1 0 0; 0 0 1; 0 0 0];
%     allcombs_new = [1 1 1; 1 1 0; 1 0 1; 0 1 1; 0 1 0; 1 0 0; 0 0 1; 0 0 0];
    allcombs2 = allcombs; allcombs2(allcombs==0) = -1;
    wh = squeeze(sum(repmat(allcombs,[1 1 size(jnk,2)])==permute(repmat(jnk,[1 1 8]),[3 1 2]),2));
    [~,mm] = max(wh);
    %%
    figure; hold on
% %     histogram(mm(1:sum(p_SSD(ind)<.05)),.5:8.5)
    histogram(mm(end-sum(p_SSD(ind)<.05)+1:end),.5:8.5,'FaceColor','k')
% %     set(gca,'xtick',1:8,'xticklabel',num2str(allcombs2))
%     h = hist(mm(end-sum(p_SSD(ind)<.05)+1:end),.5:8.5);
%     plot(1:4,h,'.k','MarkerSize',20)
    set(gca,'xtick',1:8,'xticklabel',{'+ + +';'+ + -';'+ - +';'- + +';'- + -';'+ - -';'- - +';'- - -'})
    xlabel('Direction of Modulation on Each Arm')
    ylabel('Number of Cells Significantly Differentially Modulated')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    if ~isfolder([savefolder '\ArmReplayTriggered\'])
        mkdir([savefolder '\ArmReplayTriggered\'])
    end
    set(gcf,'renderer','Painters')
    if ~isfolder([savefolder '\ArmReplayTriggered\'])
        mkdir([savefolder '\ArmReplayTriggered\'])
    end
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\HistArmPlotAllSig_' label '_' tlab{itype}])
    
        figure; hold on
%     histogram(mm(1:sum(p_SSD(ind)<.05)),.5:8.5)
    mm2 = sum(jnk);
%     histogram(-mm2(end-sum(p_SSD(ind)<.05)+1:end),-3.5:1:1,'FaceColor','k')
% %     set(gca,'xtick',1:8,'xticklabel',num2str(allcombs2))
    h = histc(-mm2(end-sum(p_SSD(ind)<.05)+1:end),-3.5:1:1);
    plot(1:4,h(1:end-1),'ok','MarkerSize',20)
    set(gca,'xtick',1:4,'xticklabel',{'+ + +';'+ + -';'+ - -';'- - -'})
    set(gca,'xlim',[.5 4.5])
    set(gca,'ylim',[0 max(h)+1])
    xlabel('Direction of Modulation on Each Arm')
    ylabel('Number of Cells Significantly Differentially Modulated')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    if ~isfolder([savefolder '\ArmReplayTriggered\'])
        mkdir([savefolder '\ArmReplayTriggered\'])
    end
    set(gcf,'renderer','Painters')
    if ~isfolder([savefolder '\ArmReplayTriggered\'])
        mkdir([savefolder '\ArmReplayTriggered\'])
    end
    helper_savefig([savefolder '\Figure2\HistArmPlotAllSigNew_' label '_' tlab{itype} '_nobar'])
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\HistArmPlotAllSigNew_' label '_' tlab{itype} '_nobar'])
    %%
      figure; hold on
%     histogram(mm(1:sum(p_SSD(ind)<.05)),.5:8.5)
    histogram(mm,.5:8.5,'FaceColor','k')
%     set(gca,'xtick',1:8,'xticklabel',num2str(allcombs2))
    set(gca,'xtick',1:8,'xticklabel',{'+ + +';'+ + -';'+ - +';'- + +';'- + -';'+ - -';'- - +';'- - -'})
    xlabel('Direction of Modulation on Each Arm')
    ylabel('Number of Cells Significantly Differentially Modulated')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    if ~isfolder([savefolder '\ArmReplayTriggered\'])
        mkdir([savefolder '\ArmReplayTriggered\'])
    end
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\HistArmPlotAllShort_' label '_' tlab{itype}])
    %%
    hh = bigfigure; 
    toplot2 = toplot;
    h = scatter3(toplot2(1,:),toplot2(2,:),toplot2(3,:),60*ones(size(toplot,2),1),...
        [repmat([0 0 0],[sum(p_SSD(ind)>=.05) 1]);repmat([1 0 0],[sum(p_SSD(ind)<.05) 1])],'LineWidth',2);
    % h.MarkerFaceColor = [.5 .5 .5];
    hold on
    axis tight
%     set(gca,'xlim',[-.4 .4],'ylim',[-.4,.4],'zlim',[-.4 .4])
    plot3(get(gca,'XLim'),[0 0],[0 0],'color',[.5 .5 .5],'LineWidth',3);
    plot3([0 0],[0 0],get(gca,'ZLim'),'color',[.5 .5 .5],'LineWidth',3);
    plot3([0 0],get(gca,'YLim'),[0 0],'color',[.5 .5 .5],'LineWidth',3);
%     scatter3(toplot(1,:),toplot(2,:),toplot(3,:),60*ones(size(toplot,2),1),[repmat([0 0 0],[sum(p_SSD(ind)>=.05) 1]);repmat([1 0 0],[sum(p_SSD(ind)<.05) 1])],'LineWidth',2);
    xlabel('Modulation to Left Arm')
    ylabel('Modulation to Right Arm')
    zlabel('Modulation to Center Arm')
    set(gca,'FontSize',18)
    set(gcf,'Position',[1.0891    0.0926    0.3297    0.4769])
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\3DPlotAll_' label '_' tlab{itype}])
    %%
    l = sum(p_SSD(ind)>=.05)+1;
    figure; hold on
    quads(2) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(1) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    quads(3) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(4) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(6) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    quads(5) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    bar(1:6,quads,'k')
    set(gca,'xtick',1:6,'xticklabel',{'L>C>R';'L>R>C';'R>L>C';'R>C>L';'C>R>L';'C>L>R'})
    xlabel('Order of Modulation of Arms (Left, Right, Center)')
    ylabel('Number of Significantly Differentially Modulated Cells')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\3DPlotQuadrantsAllSig_' label '_' tlab{itype}])
    %%
    l = 1;
    figure; hold on
    quads(2) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(1) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    quads(3) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)>toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(4) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)>toplot(3,l:end));
    quads(6) = sum((toplot(1,l:end) > toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    quads(5) = sum((toplot(1,l:end) < toplot(2,l:end)) & toplot(1,l:end)<toplot(3,l:end) & toplot(2,l:end)<toplot(3,l:end));
    bar(1:6,quads,'k')
    set(gca,'xtick',1:6,'xticklabel',{'L>C>R';'L>R>C';'R>L>C';'R>C>L';'C>R>L';'C>L>R'})
    xlabel('Order of Modulation of Arms (Left, Right, Center)')
    ylabel('Number of Significantly Differentially Modulated Cells')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\3DPlotQuadrantsAll_' label '_' tlab{itype}])
    %%
    l = sum(p_SSD(ind)>=.05)+1;
    ind2 = l:size(toplot,2);
    ind2 = ind2(mm(end-sum(p_SSD(ind)<.05)+1:end)==1);
    figure; hold on
    quads(2) = sum((toplot(1,ind2) > toplot(2,ind2)) & toplot(1,ind2)>toplot(3,ind2) & toplot(2,ind2)>toplot(3,ind2));
    quads(1) = sum((toplot(1,ind2) > toplot(2,ind2)) & toplot(1,ind2)>toplot(3,ind2) & toplot(2,ind2)<toplot(3,ind2));
    quads(3) = sum((toplot(1,ind2) < toplot(2,ind2)) & toplot(1,ind2)>toplot(3,ind2) & toplot(2,ind2)>toplot(3,ind2));
    quads(4) = sum((toplot(1,ind2) < toplot(2,ind2)) & toplot(1,ind2)<toplot(3,ind2) & toplot(2,ind2)>toplot(3,ind2));
    quads(6) = sum((toplot(1,ind2) > toplot(2,ind2)) & toplot(1,ind2)<toplot(3,ind2) & toplot(2,ind2)<toplot(3,ind2));
    quads(5) = sum((toplot(1,ind2) < toplot(2,ind2)) & toplot(1,ind2)<toplot(3,ind2) & toplot(2,ind2)<toplot(3,ind2));
    bar(1:6,quads,'k')
    set(gca,'xtick',1:6,'xticklabel',{'L>C>R';'L>R>C';'R>L>C';'R>C>L';'C>R>L';'C>L>R'})
    xlabel('Order of Modulation of Arms (Left, Right, Center)')
    ylabel('Number of Significantly Differentially Modulated Cells')
    set(gca,'FontSize',18)
    set(gcf,'Position',[ 2056          37        1175         832])
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\3DPlotQuadrantsPluses_' label '_' tlab{itype}])
    %%
    updown = [sum(toplot(:,sum(p_SSD(ind)>=.05):end)>0,2) sum(toplot(:,sum(p_SSD(ind)>=.05):end)<0,2)];
    figure; hold on
    subplot(3,1,1)
    bar(updown,'stacked')
    title(num2str(updown))
    [~,m] = max(toplot(:,sum(p_SSD(ind)>=.05):end));
    subplot(3,1,2); 
    histogram(m);
    h = hist(m,1:3);
    title(num2str(h))
    [~,n] = min(toplot(:,sum(p_SSD(ind)>=.05):end));
    subplot(3,1,3); 
    histogram(n);
    h = hist(n,1:3);
    title(num2str(h))
    set(gcf,'Position',[2319          64         472         687])
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\NumberModu_' label '_' tlab{itype}])
    %%
    
    figure; hold on
    f = scatter(toplot(1,:),toplot(2,:));
    f.CData = toplot(3,:); 
    c = colorbar;
    c.Label.String = 'Modulation to Center Arm';
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    plot(xl,[0 0],'k--')
    plot([0 0],yl,'k--')
    xlabel('Modulation to Left Arm')
    ylabel('Modulation to Right Arm')
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([savefolder '\ArmReplayTriggered\3DPlotColor_' label '_' tlab{itype}])
end