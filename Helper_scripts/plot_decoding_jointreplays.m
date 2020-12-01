function plot_decoding_jointreplays(dirs,label)

cd(dirs.homedir)
d2 = dir('*.mat');
PFCreal = [];PFCshuff = [];HPreal = [];HPshuff = [];
for id = 1:size(d2,1)    
    load(d2(id).name,[label  '_decodejointreplays_PFCreal'],[label  '_decodejointreplays_PFCshuff'],...
        [label  '_decodejointreplays_HPreal'],[label  '_decodejointreplays_HPshuff'])
    eval(['PFCreal1 = ' label  '_decodejointreplays_PFCreal;'])
    eval(['HPreal1 = ' label  '_decodejointreplays_HPreal;'])
    eval(['PFCshuff1 = ' label  '_decodejointreplays_PFCshuff;'])
    eval(['HPshuff1 = ' label  '_decodejointreplays_HPshuff;'])
    PFCreal = [PFCreal;PFCreal1];
    PFCshuff = [PFCshuff; PFCshuff1];
    HPreal = [HPreal;HPreal1];
    HPshuff = [HPshuff; HPshuff1];
end

if ~isfolder('E:\XY_matdata\Figures\ForPaper\JointReplayDecoding\')
    mkdir('E:\XY_matdata\Figures\ForPaper\JointReplayDecoding\')
end
% Joint Replay Prediction Significance
figure; hold on; histogram(nanmean(PFCshuff),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(PFCshuff)>=nanmean(PFCreal))+1)/(size(PFCshuff,2)+1);
title([label ' Joint Replay Prediction: PFC vs Shuffle, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\JointReplayDecoding\PFCvsSHuff_' label])

figure; hold on; histogram(nanmean(HPreal),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(HPreal)>=nanmean(PFCreal))+1)/(size(HPreal,2)+1);
title([label ' Joint Replay Prediction: PFC vs HP, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\JointReplayDecoding\PFCvsHP_' label])

figure; hold on; 
histogram(nanmean(HPshuff),'FaceColor','k'); 
yl = get(gca,'ylim');
plot([nanmean(nanmean(HPreal)) nanmean(nanmean(HPreal))],yl,'r','LineWidth',3)
p = (sum(nanmean(HPshuff)>=nanmean(nanmean(HPreal)))+1)/(size(HPshuff,2)+1);
title([label ' Joint Replay Prediction: HP vs Shuffle, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\JointReplayDecoding\HPvsSHuff_' label])