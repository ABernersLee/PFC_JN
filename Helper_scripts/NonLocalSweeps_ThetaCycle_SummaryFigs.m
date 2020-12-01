function NonLocalSweeps_ThetaCycle_SummaryFigs(dirs,cutoff,savefolder)

cd(dirs.homedir)
d2 = dir('*.mat');
allevents = NaN(size(d2,1),2,2,2);
%first dim is day, 2nd is which half of theta cycle, 3rd is global zero or
%shifted by sequence score, 4th is all events or nonlocal>.4 events
Thall = [];
for id = 1:size(d2,1)    
    thisdir = d2(id).name;        
    load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all') 
    
    Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
    Thall = cat(1,Thall,Th);
    clear times_armon_thetaof_headingarm_lap_thetahalf_all
    %Th 7 is from centered theta sequence, 8 is from centered global zero
    %9 is the non-local fraction
    
    ind = ~isnan(Th(:,9));
    allevents(id,1,1,1) = sum(Th(ind,7)==1);
    allevents(id,1,2,1) = sum(Th(ind,8)==1);
    allevents(id,2,1,1) = sum(Th(ind,7)==2);
    allevents(id,2,2,1) = sum(Th(ind,8)==2);
    
    ind = Th(:,9)>cutoff;
    allevents(id,1,1,2) = sum(Th(ind,7)==1);
    allevents(id,1,2,2) = sum(Th(ind,8)==1);
    allevents(id,2,1,2) = sum(Th(ind,7)==2);
    allevents(id,2,2,2) = sum(Th(ind,8)==2);
end
%%

if ~isfolder([savefolder '\Theta\'])
    mkdir([savefolder '\Theta\'])
end
% [r,p] = corr(Thall(:,9),Thall(:,8),'type','Spearman','rows','complete');
typelab = {'Theta Sequence Shifted';'Global Zero'};

for itype = 1:2
    figure; hold on
    dat1 = Thall(Thall(:,6+itype)==1,9);
    errorbar(1,nanmean(dat1),nanstd(dat1)./sqrt(sum(~isnan(dat1))),'k','LineWidth',3)
    dat2 = Thall(Thall(:,6+itype)==2,9);
    errorbar(2,nanmean(dat2),nanstd(dat2)./sqrt(sum(~isnan(dat2))),'k','LineWidth',3)
    xlim([.5 2.5])
%     if itype == 1
%         ylim([.18 .2])
%         set(gca,'ytick',[.18:.005:.2])
%     end
%     h = hist(dat1,5,'Normalization','probability');
%     [hh,c] = hist(dat2,5,'Normalization','probability');
%     figure; hold on
%     h1 = histcounts(dat1,0:.1:1)'; 
%     [h2,c] = histcounts(dat2,0:.1:1); 
%     h2 = h2';
%     b = bar([h1./sum(h1) h2./sum(h2)],'LineWidth',2,'BarWidth',1);
%     b(2).FaceColor = [0 0 0];
%     b(1).FaceColor = [1 1 1];    
%     h = boxplot2(dat1,1);
%     make_box_plot2_color(h,[0 0 0])
%     hh = boxplot2(dat2,2);
%     make_box_plot2_color(hh,[0 0 0])
    p = ranksum(dat1,dat2);
    xlabel('Half of Theta')
    ylabel('Proportion Non-local')
    set(gca,'xtick',1:2,'xticklabel',{'1st Half';'2nd Half'})
    title([typelab{itype} ' p = ' num2str(round(p,2,'significant'))])
    helper_saveandclosefig([savefolder '\Theta\NonLocalSweeps_ThetaCycle_new_'  typelab{itype}])
end

%%
siglab = {'All Theta Sweeps';['Non-local cutoff ' num2str(cutoff)]};
typelab = {'Theta Sequence Shifted';'Global Zero'};
for itype = 1:2
    for isig = 1:2
        figure; hold on
        for ihalf = 1:2
            plot(ihalf,allevents(:,ihalf,itype,isig),'+','MarkerEdgeColor',[.5 .5 .5],'MarkerSize',20)
            errorbar(ihalf,mean(allevents(:,ihalf,itype,isig)),std(allevents(:,ihalf,itype,isig))./sqrt(size(allevents,1)),'k','LineWidth',3)
        end
        plot((ones(size(allevents,1),1)*[1 2])',allevents(:,:,itype,isig)','k')
        p = signrank(allevents(:,1,itype,isig),allevents(:,2,itype,isig),'tail','left');
        if p<.05
            yl = get(gca,'ylim');
            plot([1 2],[yl(2)*1.1 yl(2)*1.1],'k')
            plot(1.5,yl(2)*1.15,'r*','MarkerSize',10)
        end
%         ylim([yl(1) yl(2)*1.2])
        xlim([.7 2.3])
        title([siglab{isig} ' ' typelab{itype} ' p = ' num2str(p)])
        helper_saveandclosefig([savefolder '\Theta\NonLocalSweeps_ThetaCycle_new_ ' siglab{isig} ' ' typelab{itype}]) 
    end
end
    
    