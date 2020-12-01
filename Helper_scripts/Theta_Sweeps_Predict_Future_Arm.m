function Theta_Sweeps_Predict_Future_Arm(dirs,cutoff,savefolder,fromend,igroup,numshuff)

cd(dirs.homedir)
d2 = dir('*.mat');
if igroup==3
    dduse = NaN(size(d2,1),1);
    for id = 1:size(d2,1)    
        load(d2(id).name,'other_cells_touse')
        if sum(other_cells_touse(:,igroup))>0
            dduse(id) = true;
        elseif sum(other_cells_touse(:,igroup))==0
            dduse(id) = false;
        end
    end
    if sum(isnan(dduse))>0
        error('day ind NaN')
    end
    d2 = d2(dduse==1);
end
    
lapinfo = cell(size(d2,1),1);
for id = 1:size(d2,1)    
    thisdir = d2(id).name;

    load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos')
    [~,armp] = max(armposindex,[],2); clear armposindex
    lappro = NaN(size(laps_singlepass,1),3);
    lapout = NaN(size(laps_singlepass,1),2);
    % lapdir = NaN(size(laps_singlepass,1),1);
    lapi = NaN(max(laps_singlepass),2);
    for ilap = 1:max(laps_singlepass)
       leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
       lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle      
       lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform

       nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
       thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
       lapi(ilap,:) = [thisarm nextarm];
    end
    lapinfo{id} = lapi;
end

%% Does HP tend to have theta sweeps of where the animal is about to go? yes
cd(dirs.homedir)

pid = NaN(size(d2,1),1);
realall = []; shuffall = [];
realall2 = []; shuffall2 = [];
% numshuff= 2000;
% fromend = [10 180];
for id = 1:size(d2,1)    
    lapi = lapinfo{id};
    thisdir = d2(id).name;
%     load(thisdir,'times_armon_thetaof_headingarm_lap')
    load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all')     
            
    times_armon_thetaof_headingarm_lap = times_armon_thetaof_headingarm_lap_thetahalf_all...
        (times_armon_thetaof_headingarm_lap_thetahalf_all(:,end)>cutoff,1:8);
    aa = times_armon_thetaof_headingarm_lap;
%     if fromend>0
        load(thisdir,'armposindex','armpos','pos','binpos')
        include = true(size(aa,1),1);
        [~,armp] = max(armposindex,[],2); clear armposindex
        [~,ii] = histc(mean([aa(:,1) aa(:,2)],2),pos(:,1)); % index of the middle of the theta sweep

        co = NaN(3,1);
        for iarm = 1:3
            fm = find(armp==iarm,1,'first');
            cutoff2 = fm+fromend(2);
            co(iarm) = cutoff2;
            cutoff1 = fm+fromend(1);
            include(armpos(ii)==iarm & (binpos(ii)<cutoff1 | binpos(ii)>cutoff2)) = false;            
        end
        
        aa = aa(include,:);
%     end
%         figure; hold on; plot(pos(:,2),pos(:,3),'k.'); 
%         plot(pos(ismember(binpos,co),2),pos(ismember(binpos,co),3),'r.')
%         helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\Sweeps_Predict_Future_Arm_' num2str(cutoff) '_cutoff' num2str(fromend(2)) ' ' d2(id).name(1:end-4)])
        
%         sum(include)
        
            
    
    realprop = NaN(size(lapi,1),1);
    for ilap = 1:max(aa(:,6))
        dat = aa(aa(:,6)==ilap,1:5);
        a = sum(dat(:,4)==unique(dat(:,5))); b = sum(dat(:,4)~=unique(dat(:,5)));
        if (a+b)>0
            realprop(ilap,1) = (a)./(a+b);
        end
    end
    shufflap = [];
    for iarm = 1:3
        shufflaps = lapi(lapi(:,1)==iarm,2); 
        h = histc(shufflaps,unique(shufflaps));
        if nchoosek(sum(h),min(h))<numshuff
            shufflap1 = NaN(size(shufflaps,1),numshuff);
            shufflap = [shufflap;shufflap1];
            disp(['Day ' num2str(id) ', Skipping arm ' num2str(iarm) ', Only ' num2str(nchoosek(sum(h),min(h))) ' perms'])
            continue
        end
%         disp(['Ran with ' num2str(nchoosek(sum(h),min(h))) ' perms'])
        currlaps = unique(aa(aa(:,3)==iarm,6));
        shufflap1 = NaN(length(shufflaps),numshuff);
        for is = 1:numshuff
            sslaps = shufflaps(randperm(length(shufflaps)));
            while sum(sslaps~=shufflaps)==0
                sslaps = shufflaps(randperm(length(shufflaps)));
            end
            if is>1; while any(sum(sslaps~=shuffsave)==0); sslaps = shufflaps(randperm(length(shufflaps))); end; shuffsave = [shuffsave sslaps]; elseif is == 1;shuffsave = sslaps; end
            for ilap = 1:length(currlaps)
                a = sum(aa(aa(:,6)==currlaps(ilap),4)==sslaps(ilap)); b = sum(aa(aa(:,6)==currlaps(ilap),4)~=sslaps(ilap));
                if (a+b)>0
                    shufflap1(ilap,is) =  (a)./(a+b);
                end
            end
        end
        clear shuffsave
        shufflap = [shufflap;shufflap1];
    end
    if isnan(nanmean(nanmean(shufflap)))
        disp(['Skipping Day ' num2str(id)])
        continue
    end
    pid(id) = (sum(nanmean(shufflap)>=nanmean(realprop))+1)/(size(shufflap,2)+1);
    1-binocdf(sum(realprop>nanmean(nanmean(shufflap)))-1,size(realprop,1),.05);
    realall = [realall; nanmean(realprop)];
    shuffall = [shuffall;nanmean(shufflap)];
    realall2 = [realall2; realprop(~isnan(realprop))];
    shuffall2 = [shuffall2;shufflap(~isnan(realprop),:)];
end

p = (sum(nanmean(shuffall)>=nanmean(realall))+1)/(size(shuffall,2)+1);
figure; histogram(nanmean(shuffall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realall) nanmean(realall)],yl,'r-','LineWidth',3)
if p<.05
    text(nanmean(realall)*.98,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','r')
else
    text(nanmean(realall)*.98,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','k')
end
ylabel('Shuffle Counts')
xlabel('Proportion going to Heading Arm')
set(gca,'FontSize',18)
% xlim([.5 .7])
p2 = 1-binocdf(nansum(pid<.05)-1,sum(~isnan(pid)),.05);
title(['Days - Difference of ' num2str(round((nanmean(realall)-nanmean(nanmean(shuffall)))*100,2)) '%, ' num2str(nansum(pid<.05)) '/' num2str(sum(~isnan(pid))) ' days sig (p=' num2str(round(p2,2,'significant')) ')'])
if ~isfolder([savefolder '\Figure4\'])
    mkdir([savefolder '\Figure4\'])
end
set(gcf,'renderer','Painters')
set(gcf,'Position',[  2037         255         650         522])
helper_saveandclosefig([savefolder '\Figure4\Sweeps_Predict_Future_Arm_' num2str(cutoff) '_from' num2str(fromend(1)) '_to' num2str(fromend(2)) '_Days' '_numshuff' num2str(numshuff)])

1-binocdf(sum(realall>.5)-1,11,.05)


p = (sum(nanmean(shuffall)>=nanmean(realall))+1)/(size(shuffall,2)+1);
figure; 
% histogram(nanmean(shuffall),'FaceColor','k'); 
hold on; 
[h,c] = hist(nanmean(shuffall),20); 
plot(c,h,'.-k','LineWidth',2,'MarkerSize',20)

yl = get(gca,'ylim');
plot([nanmean(realall) nanmean(realall)],yl,'r-','LineWidth',3)
if p<.05
    text(nanmean(realall)*.96,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','r')
else
    text(nanmean(realall)*.96,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','k')
end
ylabel('Shuffle Counts')
xlabel('Proportion going to Heading Arm')
set(gca,'FontSize',18)
xlim([min(c)-.01 nanmean(realall)+.01])
p2 = 1-binocdf(nansum(pid<.05)-1,sum(~isnan(pid)),.05);
title(['Days - Difference of ' num2str(round((nanmean(realall)-nanmean(nanmean(shuffall)))*100,2)) '%, ' num2str(nansum(pid<.05)) '/' num2str(sum(~isnan(pid))) ' days sig (p=' num2str(round(p2,2,'significant')) ')'])
if ~isfolder([savefolder '\Figure4\'])
    mkdir([savefolder '\Figure4\'])
end
set(gcf,'renderer','Painters')
set(gcf,'Position',[  2037         255         650         522])
helper_saveandclosefig([savefolder '\Figure4\Sweeps_Predict_Future_Arm_' num2str(cutoff) '_from' num2str(fromend(1)) '_to' num2str(fromend(2)) '_Days' '_numshuff' num2str(numshuff) '_nobar'])


if 0
p = (sum(nanmean(shuffall2)>=nanmean(realall2))+1)/(size(shuffall2,2)+1);
figure; histogram(nanmean(shuffall2),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realall2) nanmean(realall2)],yl,'r-','LineWidth',3)
if p<.05
    text(nanmean(realall2)*.99,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','r')
else
    text(nanmean(realall2)*.99,yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))],'Color','k')
end
ylabel('Shuffle Counts')
xlabel('Proportion going to Heading Arm')
set(gca,'FontSize',18)
% xlim([.5 .7])
title(['Laps - Difference of ' num2str(round((nanmean(realall)-nanmean(nanmean(shuffall)))*100,2)) '%'])
if ~isfolder([savefolder '\Figure4\'])
    mkdir([savefolder '\Figure4\'])
end
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure4\Sweeps_Predict_Future_Arm_' num2str(cutoff) '_from' num2str(fromend(1)) '_to' num2str(fromend(2)) '_Laps' '_numshuff' num2str(numshuff)])
end

% else
%     helper_saveandclosefig([savefolder '\Theta\Sweeps_Predict_Future_Arm_' num2str(cutoff)])
% end
%%
if 0 
%% ridge regression, prospective coding is a better predictor until you raise the threshold of theta sweep up, then isnt anymore (still right direction)
% label = 'RP'; (just using sig mod cells doesnt change it much)
cd(dirs.homedir)
d2 = dir('*.mat');
k = 1e-5; %0:1e-5:5e-3;
bsall = [];
modind1 = [0 .2];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap','other_cells','spikedata')
    load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','RP_armsig') 
    times_armon_thetaof_headingarm_lap = times_armon_thetaof_headingarm_lap_thetahalf_all...
        (times_armon_thetaof_headingarm_lap_thetahalf_all(:,end)>.5,1:8);
    %     pfcfr = prsopectiveTheta_pfcfr;    
        aa = times_armon_thetaof_headingarm_lap;
        label = 'RP';
        load(thisdir, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'])
        eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])
    %     load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap',[label '_Cand_sig_modu_include'])
    %     eval(['Cand = ' label '_Cand_sig_modu_include;'])
    %     pfcfr = prsopectiveTheta_pfcfr(:,Cand(:,1)==1);    

        spikes = spikedata(ismember(spikedata(:,2),other_cells),:);
        modind = [aa(:,2)+modind1(1) aa(:,2)+modind1(2)];                
        ms = NaN(size(modind,1),length(other_cells));
        for ievent = 1:size(modind,1)
            modspikes = spikes(spikes(:,1)>=modind(ievent,1) & spikes(:,1)<modind(ievent,2),2);
            ms(ievent,:) = histc(modspikes,other_cells);
        end
%                 pfcfr = ms(:,armsig(:,1)<.05);
                pfcfr = ms;
%         other_cells = other_cells(armsig(:,1)<.05);
    
%     bs = NaN(size(pfc,2),3,length(k),3);
%     bs = NaN(size(pfc,2),3,3); 
    bs = NaN(size(pfc,2),2,3); 
    for iarm = 1:3
        
        Y = aa(aa(:,3)==iarm,4:5); % 4 is theta of, 5 is heading arm
        D = x2fx(Y,'interaction');
        pfc = pfcfr(aa(:,3)==iarm,:); 
        pfc2 = zscore(pfc);
        for icell = 1:size(pfc,2)
            if sum(pfc(:,icell))>0                
%                 [b,~,~] = glmfit(Y,pfc2(:,icell));        
%                 b = ridge(pfc(:,icell),D(:,2:end),k);                        
                b = ridge(pfc(:,icell),Y,k);        
%                 bs(icell,:,:,iarm) = b;
                bs(icell,:,iarm) = b;
            end
        end
        disp(['Done with arm ' num2str(iarm)])
    end    
    bsall = cat(1,bsall,bs);
    disp(['Done with day ' num2str(id)])
end
bsall2 = bsall;
bsall = abs(bsall);
nanmean(nanmean(bsall,3),1)
signrank(nanmean(bsall(:,1,:),3),nanmean(bsall(:,2,:),3))



%%
[p,tbl,stats]= kruskalwallis(nanmean(bsall,3))
p = multcompare(stats) 
%%
[p,tbl,stats]= kruskalwallis(bsall(:,:,1))
p = multcompare(stats) 
[p,tbl,stats]= kruskalwallis(bsall(:,:,2))
p = multcompare(stats) 
[p,tbl,stats]= kruskalwallis(bsall(:,:,3))
p = multcompare(stats) 
%%


%% results from single factor test
 [p,tbl,stats]= kruskalwallis(bsall(:,1,:)-bsall(:,2,:)); % sig difference across arms
 p = multcompare(stats) %center arm is sig more than other two arms
 
 p1 = ranksum(bsall(:,1,1),bsall(:,2,1))
 nanmedian(bsall(:,1,1))-nanmedian(bsall(:,2,1)) % left arm is non-significantly more heading than theta
 
 p3 = ranksum(bsall(:,1,3),bsall(:,2,3))
 nanmedian(bsall(:,1,3))-nanmedian(bsall(:,2,3)) % right arm is non-significantly more heading than theta (switches but not sig when abs(beta))
 
 p1 = ranksum(bsall(:,1,2),bsall(:,2,2))
 nanmedian(bsall(:,1,2))-nanmedian(bsall(:,2,2)) % center arm is significantly more about theta than heading

 %% results from interaction test

 [p,tbl,stats]= kruskalwallis(nanmean(bsall,3)); % no difference when average over arms
 
[p,tbl,stats]= kruskalwallis(bsall(:,:,1)); % no difference
 
 [p,tbl,stats]= kruskalwallis(bsall(:,:,2)); % on center arm, significantly more theta than heading or interaction
 p = multcompare(stats) 
 
 [p,tbl,stats]= kruskalwallis(bsall(:,:,3)); % no difference
 
 [p,tbl,stats]= kruskalwallis(bsall(:,1,:)); % significantly more theta on center arm than other two arms (not when I take abs value)
 p = multcompare(stats) 
end