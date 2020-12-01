function plot_prospective_theta_modu(thisdir)

load(thisdir,'prospectiveTheta_toplot','prospectiveTheta_armcombs','other_cells')

dat = nanmean(prospectiveTheta_toplot,4);
dat2 = reshape(prospectiveTheta_armcombs,[size(prospectiveTheta_armcombs,1) 3 2]);
dat(isnan(dat2)) = NaN;
datsem = nanstd(prospectiveTheta_toplot,[],4)./sqrt(sum(~isnan(prospectiveTheta_toplot),4));

cs = varycolor(6);
cs(cs==1)= .7;
cs = cs([1 4 3 5 2 6],:);

armlab = {'Left';'Center';'Right'};
for icell = 1:size(dat,1)
   figure; hold on
   yl = NaN(3,2);
   for iarm = 1:3
       otharms = setdiff(1:3,iarm);
       subplot(3,1,iarm); hold on
       a1 = errorbar(1,dat(icell,iarm,1),datsem(icell,iarm,1),'.','color',cs((iarm-1)*2+1,:),'LineWidth',3,'MarkerSize',30);
       a2 = errorbar(2,dat(icell,iarm,2),datsem(icell,iarm,2),'.','color',cs((iarm-1)*2+2,:),'LineWidth',3,'MarkerSize',30);
       legend([a1 a2],{['Theta Sweeps of ' armlab{otharms(1)} ' Arm'],['Theta Sweeps of ' armlab{otharms(2)} ' Arm']},'Location','EastOutside')       
       ylabel(['FR on ' armlab{iarm} ' Arm'])
       set(gca,'xlim',[.5 2.5],'xtick',[])
       set(gca,'FontSize',18)
       yl(iarm,:) = get(gca,'ylim');
   end         
   for iarm = 1:3
       subplot(3,1,iarm)
       set(gca,'ylim',[min(yl(:,1)) max(yl(:,2))])
   end

   suptitle(['Cell ' num2str(other_cells(icell)) ', Rat ' num2str(thisdir(1)) ' ' thisdir(3:end-4)])
   set(gca,'FontSize',18)
   set(gcf,'Position',[  680    91   980   887])
%    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\PFCmodulation_' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell))])

    if ~isfolder(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\'])
        mkdir(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\'])
    end
   helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\PFCmodulation'])
end

% b = dat;
% bb = reshape(b,[size(b,1) 2*size(b,2)]);
% p = kruskalwallis(bb);
% 
% for icell = 1:size(prospectiveTheta_toplot,1)
%     b = permute(squeeze(prospectiveTheta_toplot(icell,:,:,:)),[3 1 2]);
%     bb = reshape(b,[size(b,1) 2*size(b,2)]);
%     ps(icell) = kruskalwallis(bb);
% end


