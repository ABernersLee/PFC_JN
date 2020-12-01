function plot_decoding_nextarm_usingFR(dirs)

cd(dirs.homedir)
d2 = dir('*.mat');
PFCreal = [];PFCshuff = [];HPreal = [];HPshuff = [];
for id = 1:3 %size(d2,1)    
    load(d2(id).name,'decodenextarmFR_PFCreal','decodenextarmFR_PFCshuff', ...
        'decodenextarmFR_HPreal','decodenextarmFR_HPshuff')
    PFCreal1 = decodenextarmFR_PFCreal;
    HPreal1 = decodenextarmFR_HPreal;
    PFCshuff1 = decodenextarmFR_PFCshuff;
    HPshuff1 = decodenextarmFR_HPshuff;
    PFCreal = [PFCreal;PFCreal1];
    PFCshuff = [PFCshuff; PFCshuff1];
    HPreal = [HPreal;HPreal1];
    HPshuff = [HPshuff; HPshuff1];
end

if ~isfolder('E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\')
    mkdir('E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\')
end
figure; hold on; histogram(nanmean(PFCshuff),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(PFCshuff)>=nanmean(PFCreal))+1)/(size(PFCshuff,2)+1);
title(['Prospective FR Prediction: PFC vs Shuffle, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\PFCvsSHuff'])

figure; hold on; histogram(nanmean(HPreal),'FaceColor','k'); yl = get(gca,'ylim');plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r','LineWidth',3)
p = (sum(nanmean(HPreal)>=nanmean(PFCreal))+1)/(size(HPreal,2)+1);
title(['Prospective FR Prediction: PFC vs HP, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\PFCvsHP'])

figure; hold on; 
histogram(nanmean(HPshuff),'FaceColor','k'); 
yl = get(gca,'ylim');
plot([nanmean(nanmean(HPreal)) nanmean(nanmean(HPreal))],yl,'r','LineWidth',3)
p = (sum(nanmean(HPshuff)>=nanmean(nanmean(HPreal)))+1)/(size(HPshuff,2)+1);
title(['Prospective FR Prediction: HP vs Shuffle, p = ' num2str(p)])
set(gca,'FontSize',18)
set(gcf,'Position',[2069         208         960         572])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\HPvsSHuff'])