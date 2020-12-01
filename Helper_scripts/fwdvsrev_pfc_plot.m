function fwdvsrev_pfc_plot(dirs,igroup,savefolder)
%%
cd(dirs.homedir)
d2 = dir('*.mat');
AA = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    alldatz = fwdvsrev_pfc(thisdir,'RP',igroup);
    AA = [AA;alldatz];
end
nanmedian(AA);
nanmean(AA);
p = signrank(AA(:,1),AA(:,2),'tail','right');
signrank(AA(:,1),AA(:,2));

figure; hold on
plot(1,AA(:,1),'.')
errorbar(1,mean(AA(:,1)),std(AA(:,1)/sqrt(size(AA,1))),'k','LineWidth',3)
plot(2,AA(:,2),'.')
errorbar(2,mean(AA(:,2)),std(AA(:,2)/sqrt(size(AA,1))),'k','LineWidth',3)
ylim([min(min(AA)) max(max(AA))])
xlim([.7 2.3])
set(gca,'xtick',1:2,'xticklabel',{'Away';'Towards'})
ylabel('Z-Scored FR')
title([' p = ' num2str(round(p,3)) ' N = ' num2str(size(AA,1)) ''])
if ~isfolder([savefolder '\Supplement\'])
    mkdir([savefolder '\Supplement\'])
end
helper_saveandclosefig([savefolder '\Supplement\Away_Vs_Towards_pfc_N'  num2str(size(AA,1)) '_Allsig_ripplesigcells'])