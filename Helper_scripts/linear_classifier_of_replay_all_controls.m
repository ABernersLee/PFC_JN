function linear_classifier_of_replay_all_controls(label,savelab,figlab,igroup,savefolder)

cd('E:\XY_matdata\AllDays')
d2 = dir('*.mat');

type = 1;

p = []; dist = []; shuff = []; mdist = []; mshuff = [];
for id =  1:size(d2,1) %[1 3 6:size(d2,1)] % should change this?
    thisdir = d2(id).name;
    
    load(thisdir,'other_cells_touse')
    if sum(other_cells_touse(:,igroup))==0
        continue
    end
    
    load(thisdir,[label '_' savelab '_Classifier_type' num2str(type) '_p'],...
        [label '_' savelab '_Classifier_type' num2str(type) '_r'],...
        [label '_' savelab '_Classifier_type' num2str(type) '_s'])
    eval(['dist1 = ' label '_' savelab '_Classifier_type' num2str(type) '_r;'])
    eval(['sdist1 = ' label '_' savelab '_Classifier_type' num2str(type) '_s;'])
    eval(['p1 = ' label '_' savelab '_Classifier_type' num2str(type) '_p;'])
    clear([label '_' savelab '_Classifier_type' num2str(type) '_p'],...
        [label '_' savelab '_Classifier_type' num2str(type) '_r'],...
        [label '_' savelab '_Classifier_type' num2str(type) '_s'])
    p = cat(1,p,p1);
    dist = cat(1,dist,dist1(:));
    mdist = cat(1,mdist,nanmean(dist1(:)));
    shuff = cat(1,shuff,sdist1);
    mshuff = cat(1,mshuff,nanmean(sdist1(:)));
end

savelab = [savelab figlab];
lab = {'linear discriminant analysis';'multinomial logistic regression';'psuedolinear discriminant analysis'};
c = [0 0 0];

% p1 = ranksum(dist,shuff,'tail','right');
p1 = (sum(mean(shuff)>=mean(dist))+1)/(size(shuff,2)+1);
% p2 = signrank(mdist,mean(shuff,2),'tail','right');
p2 = signrank(mdist,mshuff,'tail','right');
% [~,p2] = ttest(mdist,mshuff,'tail','right');

figure; hold on;

subplot(3,1,1), hold on
histogram(mean(shuff),'FaceColor','k')
yl = get(gca,'ylim'); xl = get(gca,'xlim');
plot([mean(dist) mean(dist)],yl,'r','LineWidth',5)
title('Average Across All Data')
t1 = text(xl(2)*.99,yl(1)+(range(yl)*.8),['p = ' num2str(round(p1,2,'significant'))]);
if p1<.05
    t1.Color = 'r';
end


subplot(3,1,2), hold on
plot([dist mean(shuff,2)]','k.-','LineWidth',2,'MarkerSize',20)
title('Mean Within Each Day')
yl = get(gca,'ylim');
ylim([yl(1)*.85 yl(2)*1.04])
xlim([.8 2.2])
yl = get(gca,'ylim'); xl = get(gca,'xlim');
t2 = text(1.7,yl(1)+(range(yl)*.8),['p = ' num2str(round(p2,2,'significant'))]);
if p2<.05
    t2.Color = 'r';
end
set(gca,'xtick',1:2,'xticklabel',{'Data';'Shuffle'})

subplot(3,1,3)
p3 = (1-binocdf((sum(p<.05)-1),length(p),.05));
pp = pie([sum(p>=.05) sum(p<.05)],{['Days Not Significant (' num2str(sum(p>=.05)) ')'],['Days Significant, (' num2str(sum(p<.05)) '), p = ' num2str(round(p3,2,'significant'))]});
pp(1).FaceColor = 'k';
if length(pp)>2
pp(3).FaceColor = [.5 .5 .5];
if p3<.05
    pp(4).Color = 'r';
    pp(3).FaceColor = 'r';
end
end
if sum(p>=.05)==0
    pp(1).FaceColor = 'r';
end
suptitle([lab{type} ' ' savelab ' ' label])

set(gcf,'Position',[680   244   614   734])
if ~isfolder([savefolder '\Classify\'])
    mkdir([savefolder '\Classify\'])
end
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Classify\' savelab '_' lab{type} ' ' label])


    