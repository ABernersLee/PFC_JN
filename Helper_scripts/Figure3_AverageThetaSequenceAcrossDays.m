function Figure3_AverageThetaSequenceAcrossDays(dirs,savefolder)
cd(dirs.homedir)
load([dirs.homedir '/matfiles/theta_sequence_zero.mat'],'theta_sequence_zero')
d2 = dir('*.mat');
%average theta sequence across days
SpeedCutoff = 5;
newmat = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'theta_globalzero') %,'times_armon_thetaof_headingarm_lap_thetahalf_all');
%     shift = draft_thetaseq_again(thisdir,ach,armposindex,curp,oc,pc,otherplots);
    
    load(thisdir,'HP_Theta')
    HP_Theta_globalzero = HP_Theta;
    HP_Theta_globalzero(:,4) = mod(HP_Theta_globalzero(:,4)-theta_globalzero,360);
%     theta = HP_Theta_globalzero;
    HP_Theta_globalzero2 = HP_Theta_globalzero;
    HP_Theta_globalzero2(:,4) = mod(HP_Theta_globalzero2(:,4)-theta_sequence_zero,360);
    theta = HP_Theta_globalzero2;
    

    ind = find(theta(:,4)>mod(359,360) | theta(:,4)<1);
    ind(diff(ind)==1) = [];
    anglechange = theta(ind,1);
    as = [anglechange(1:end-1) anglechange(2:end)];
    ind2 = (as(:,2)-as(:,1))<.1 | (as(:,2)-as(:,1))>.2;    
    as(ind2,:) = [];  
    achange = as;    
    if ~isfolder([savefolder '\Figure3\'])
        mkdir([savefolder '\Figure3\'])
    end
    if ~isfolder([savefolder '\Theta\'])
        mkdir([savefolder '\Theta\'])
    end

    load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos','RP_CandEventTimes')
    excl = false(size(pos,1),1);
    for icd = 1:size(RP_CandEventTimes(:,1))
        excl = excl | (pos(:,1)>=RP_CandEventTimes(icd,1) & pos(:,1)<=RP_CandEventTimes(icd,2));
    end

    [~,armp] = max(armposindex,[],2); 
    lappro = NaN(size(laps_singlepass,1),3);
    lapout = NaN(size(laps_singlepass,1),2);
    % lapdir = NaN(size(laps_singlepass,1),1);
    for ilap = 1:max(laps_singlepass)
       leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
       lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first')-1; % 6 is arrive middle      
       lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform

       nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
       thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
       lappro(lapind,:) = [ilap*ones(length(lapind),1) thisarm*ones(length(lapind),1) nextarm*ones(length(lapind),1)];
       lapout(lapout1,:) = [ilap*ones(length(lapout1),1) nextarm*ones(length(lapout1),1)];
    end
    lappro(vel<SpeedCutoff | dirdat==1 | excl==1,:) = NaN; %dirdat in is 0, out is 1
    lapout(vel<SpeedCutoff | dirdat==0 | excl==1,:) = NaN; %getting rid of high velocity and ripple periods


    [~,~,i2] = histcounts(achange(:,1),pos(:,1));
    ind = [i2(1:end-1)+1 i2(2:end)];
    ind = [ind; ind(end)+1 size(pos,1)];
    achange(i2==0,:) = [];

    prospectivecycles = NaN(size(achange,1),3);
    outcycle = NaN(size(achange,1),2);
    currpos = NaN(size(achange,1),1);
    for itheta = 1:size(achange,1)
        prospectivecycles(itheta,:) = mode(lappro(ind(itheta,1):ind(itheta,2),:));
        currpos(itheta,:) = binpos(ind(itheta,1),1);
        outcycle(itheta,:) = mode(lapout(ind(itheta,1):ind(itheta,2),:));
    end

    torun = ~isnan(outcycle(:,1)) | ~isnan(prospectivecycles(:,1)); 
    ach = achange(torun,:);
    oc = outcycle(torun,:);
    pc = prospectivecycles(torun,:); %lap and arm
    curp = currpos(torun,:);

    %%%%%%%%% makes theta sequence figures and calculates shift to center
    %%%%%%%%% sequence (figures for average theta sequence and shifted theta
    %%%%%%%%% sequence)

    newmat1 = draft_thetaseq_again2(thisdir,ach,armposindex,curp,oc,pc);
    
    newmat = cat(3,newmat,newmat1);
    figure; imagesc(nanmean(newmat1,3)); axis xy
    title([thisdir(1) thisdir(3:end-4)])
    colormap gray; axis xy
    cm = colormap; cm2 = cm(end:-1:1,:);
    set(gca,'colormap',cm2)
%     cl = get(gca,'clim');
    set(gca,'clim',[.01 .03])
    ylabel('Distance from Current Position (cm)')
    xlabel('Theta Phase (degrees)')
    colorbar 
    set(gca,'FontSize',18)
    set(gcf,'Position',[1998          73         635         504])
    helper_saveandclosefig([savefolder '\Theta\Average_Theta_Sequence_Across_Days_' [thisdir(1) thisdir(3:end-4)]])
    disp(id)
end

if 0 % for making the theta_sequence_zero
matmat = nanmean(newmat,3);
wm = NaN(size(matmat,2),1);
tm = 1:size(matmat,2);
for ij = tm
    wm(ij) = [1:size(newmat,1)]*(matmat(:,ij)/sum(matmat(:,ij)));
end

bins = 1:size(matmat,1);
sl = NaN(size(matmat,1),1);
for i = bins
   ind = mod(tm+i,size(matmat,1))+1;   
%    pp = polyfit(tm(1:size(Ind,2)/2)',wm(ind(1:size(Ind,2)/2)),1);    
   pp = polyfit(tm',wm(ind),1);
   sl(i) = pp(1);   
end
[~,mx] = max(sl);
theta_sequence_zero = mod((360/125)*5*mx,360);
cd matfiles
save('theta_sequence_zero','theta_sequence_zero')
cd(dirs.homedir)

shifted_mat = newmat(:,mod(tm+bins(mx),size(matmat,1))+1,:);
end
shifted_mat = newmat;
% cd matfiles
% save('AverageThetaPhase.mat','newmat')
%%
figure; hold on

imagesc(nanmean(shifted_mat,3));  
colormap gray; axis xy
cm = colormap; cm2 = cm(end:-1:1,:);
set(gca,'colormap',cm2)
cl = get(gca,'clim');
set(gca,'clim',[0.01 .04])
set(gca,'xtick',4:5:size(shifted_mat,2),'xticklabel',[25:25:125]*360/125)
set(gca,'ytick',[1:2:size(shifted_mat,1)],'yticklabel',[1:2:size(shifted_mat,1)]-floor(size(shifted_mat,1)/2),'FontSize',18)
xlim([1 size(shifted_mat,2)])
ylim([1 size(shifted_mat,1)])
ylabel('Distance from Current Position (cm)')
xlabel('Theta Phase (degrees)')
colorbar 
set(gca,'FontSize',18)
set(gcf,'Position',[1998          73         635         504])
suptitle(['Average Theta Sequence'])

helper_saveandclosefig([savefolder '\Figure3\Average_Theta_Sequence_Across_Days_04'])


