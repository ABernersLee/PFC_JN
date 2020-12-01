function plot_prospective_FR(thisdir,velcutoff,numbins)

%formerly prospective_coding_fields
load(thisdir,'other_cells','spikedata','vel','behave_change_log','laps_singlepass','linposcat','armpos')


% numbins = 1; %was 10 (changed 10/22/18) change back probably, then was 1
s = spikedata(ismember(spikedata(:,2),other_cells),:);
fr = zeros(numbins,max(laps_singlepass),length(other_cells));       
thearms = NaN(max(laps_singlepass),2);
for ilap = 1:max(laps_singlepass)
    
    %normal
   leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
   lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first')-1; % 6 is arrive middle
   
%    lapind = find(laps_singlepass==ilap,1,'first'):find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first')-1; % 6 is arrive middle
   
    %testing
%    leave = find(behave_change_log(:,4) & laps_singlepass == ilap,1,'first'); %4 leave lick
%    lapind = find(laps_singlepass==ilap,1,'first'):leave;
   
   nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
   thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
%    lapind = leave:find(armpos(1:find(laps_singlepass == ilap,1,'last')==thisarm,1,'last')); %thought about changing to on 10/22

   lapind(vel(lapind)<velcutoff) = []; %testing
   
   lappos = linposcat(lapind);
   test = ((lappos-min(linposcat(armpos==thisarm)))/range(linposcat(armpos==thisarm)))*100;
   posbin = ceil(test./(100/numbins));
   
   
   for icell = 1:length(other_cells)
       for ipos = 1:max(posbin)
           fr(ipos,ilap,icell) = sum(s(:,2)==other_cells(icell) & ismember(s(:,3),lapind(posbin==ipos)))./(sum(posbin==ipos)/60);
       end
   end
   thearms(ilap,:) = [thisarm nextarm];
end


toplotfd = NaN(size(fr,3),size(fr,1),3,2);
toplotfd_sem = NaN(size(fr,3),size(fr,1),3,2);
for iarm = 1:3
   otharms = setdiff(1:3,iarm);
   if sum(thearms(:,1)==iarm)>0 %sum(thearms(:,1)==iarm & thearms(:,2)==otharms(1))>0 && sum(thearms(:,1)==iarm & thearms(:,2)==otharms(2))>0             
           
       armslabel = thearms(thearms(:,1)==iarm,2);
       fr2 = fr(:,thearms(:,1)==iarm,:);
       
       toplotfd(:,:,iarm,1) = squeeze(nanmean(fr2(:,armslabel==otharms(1),:),2))';
       toplotfd(:,:,iarm,2) = squeeze(nanmean(fr2(:,armslabel==otharms(2),:),2))';
       toplotfd_sem(:,:,iarm,1) = squeeze(nanstd(fr2(:,armslabel==otharms(1),:),[],2))'./sqrt(sum(armslabel==otharms(1)));
       toplotfd_sem(:,:,iarm,2) = squeeze(nanstd(fr2(:,armslabel==otharms(2),:),[],2))'./sqrt(sum(armslabel==otharms(2)));       
   end
end


% save(thisdir,'prospectiveFR_percdiff','prospectiveFR_perchalf',...
%     'prospectiveFR_whicharm','prospectiveFR_toplotfd','prospectiveFR_armcombs','-append')

cs = varycolor(6);
cs(cs==1)= .7;
cs = cs([1 4 3 5 2 6],:);

armlab = {'Left';'Center';'Right'};

for icell = 1:size(fr,3)
   forxl = (nansum(nansum(~isnan(toplotfd(icell,:,:,:)),3),4))~=0;
   figure; hold on
   yl = NaN(3,2);
   for iarm = 1:3
       otharms = setdiff(1:3,iarm);
       subplot(3,1,iarm); hold on
       a1 = plot(toplotfd(icell,:,iarm,1),'color',cs((iarm-1)*2+1,:),'LineWidth',3);
       fwd = toplotfd(icell,:,iarm,1)+toplotfd_sem(icell,:,iarm,1);
       rev = toplotfd(icell,end:-1:1,iarm,1)-toplotfd_sem(icell,end:-1:1,iarm,1);
       frind = find(~isnan(fwd));
       patch([frind frind(end:-1:1)],[fwd(~isnan(fwd)) rev(~isnan(rev))], ...
           cs((iarm-1)*2+1,:),'FaceAlpha',.3,'EdgeAlpha',0)
       
       a2 = plot(toplotfd(icell,:,iarm,2),'color',cs((iarm-1)*2+2,:),'LineWidth',3);       
       fwd = toplotfd(icell,:,iarm,2)+toplotfd_sem(icell,:,iarm,2);
       rev = toplotfd(icell,end:-1:1,iarm,2)-toplotfd_sem(icell,end:-1:1,iarm,2);
       frind = find(~isnan(fwd));
       patch([frind frind(end:-1:1)],[fwd(~isnan(fwd)) rev(~isnan(rev))], ...
           cs((iarm-1)*2+2,:),'FaceAlpha',.3,'EdgeAlpha',0)
%        legend([a1 a2],{['Arm ' num2str(otharms(1)) ' bound'],['Arm ' num2str(otharms(2)) ' bound']})       
       legend([a1 a2],{['Going to the ' armlab{otharms(1)} ' Arm'],['Going to the ' armlab{otharms(2)} ' Arm']},'Location','eastoutside')       
       ylabel(['FR on ' armlab{iarm} ' Arm'])
       set(gca,'xlim',[find(forxl,1,'first') find(forxl,1,'last')])
       set(gca,'FontSize',18)
       yl(iarm,:) = get(gca,'ylim');
   end      
   for iarm = 1:3
       subplot(3,1,iarm)
       set(gca,'ylim',[min(yl(:,1)) max(yl(:,2))])
   end
   xlabel('Position Bin Across Arm')   
   suptitle(['Cell ' num2str(other_cells(icell)) ', Rat ' num2str(thisdir(1)) ' ' thisdir(3:end-4)])
   set(gca,'FontSize',18)
   set(gcf,'Position',[1950          69        1035         781])
%    set(gcf,'Position',[ 1950         -11        1544         861])
%    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell))])

     if ~isfolder(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\'])
        mkdir(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\'])
    end
   helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\Examples\' num2str(thisdir(1:end-4)) '_Cell' num2str(other_cells(icell)) '\ProspectiveCoding'])
end

disp('Done with plot_prospective_FR')