function [times_armon_thetaof_headingarm_lap_thetahalf_all,armsig_mx] =  prospective_theta(thisdir,SpeedCutoff,exampleplots,otherplots)
% now outdated!!!! use Theta_Shifts_Sequences_Xcorr
%formerly make_Theta_Seq_forXWdat, 

% SpeedCutoff = 10;
load(thisdir,'HP_Theta')
theta = HP_Theta;
load(thisdir,'hp_cells','hpinterneurons','spikedata')
hpcells = hp_cells(~ismember(hp_cells,hpinterneurons));
spikedata = spikedata(ismember(spikedata(:,2),hpcells),:);

[~,c] = histc(spikedata(:,1),theta(:,1));

[h,cc] = hist(theta(c(c~=0),4),1:2:360);
% [~,mm] = max(h);
test = smoothts([h h h],'b',10);
smthed = test(size(h,2)+1:size(h,2)*2);
[~,mm] = max(smthed);
mx = cc(mm);

if otherplots
figure; hold on
hist(theta(c(c~=0),4),180);
plot(cc,smthed,'g','LineWidth',3)
yl = get(gca,'ylim');
plot([mx mx],yl,'k-','LineWidth',3)
text(mx*1.01,yl(2)*.99,num2str(mx))
title(['HP population activity to HP theta ' num2str(thisdir(1:end-4))])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\HPpopulationactivity_' num2str(thisdir(1:end-4))])
end

% ind = find(theta(:,4)>358 & theta(:,4)<360);
ind = find(theta(:,4)>mod((mx-2),360) & theta(:,4)<mx);
ind(diff(ind)==1) = [];
anglechange = theta(ind,1);
as = [anglechange(1:end-1) anglechange(2:end)];
ind2 = (as(:,2)-as(:,1))<.1 | (as(:,2)-as(:,1))>.2;    
as(ind2,:) = [];  
achange = as;    

%save out if going to use again
HP_Theta_globalzero = HP_Theta;
HP_Theta_globalzero(:,4) = mod(HP_Theta_globalzero(:,4)-mx,360);
% theta_globalzero = mx;
% save(thisdir,'HP_Theta_globalzero','theta_globalzero','-append')

%NOW adjusted to "global zero"
theta = HP_Theta_globalzero;

% armsig_mx = crosscov_thetatimes(theta,thisdir);

if otherplots
    [h,cc] = hist(theta(c(c~=0),4),1:2:360);
    test = smoothts([h h h],'b',10);
    smthed = test(size(h,2)+1:size(h,2)*2);
%     [~,mm] = max(smthed);
%     mx = cc(mm);
    figure; hold on
    hist(theta(c(c~=0),4),180);
    plot(cc,smthed,'g','LineWidth',3)
    yl = get(gca,'ylim');
%     plot([mx mx],yl,'k-','LineWidth',3)
%     text(mx*1.01,yl(2)*.99,num2str(mx))
    title(['Adjusted to global zero ' num2str(thisdir(1:end-4))])
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\HPpopulationactivity_adjusted_' num2str(thisdir(1:end-4))])
end

load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos')
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
% lapall(~isnan(lapall(:,1)),2) = armpos(~isnan(lapall(:,1)));
% lapall = [lapall lapdir];
lappro(vel<SpeedCutoff | dirdat==1,:) = NaN; %dirdat in is 0, out is 1
lapout(vel<SpeedCutoff | dirdat==0,:) = NaN;
% lapall(vel<SpeedCutoff,:) = NaN;


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
% torun = ~isnan(outcycle(:,1));
% torun = ~isnan(prospectivecycles(:,1));
ach = achange(torun,:);
oc = outcycle(torun,:);
pc = prospectivecycles(torun,:); %lap and arm
curp = currpos(torun,:);


shift = draft_thetaseq_again(thisdir,ach,armposindex,curp,oc,pc,otherplots);

HP_Theta_globalzero2 = theta;
HP_Theta_globalzero2(:,4) = mod(HP_Theta_globalzero2(:,4)+shift,360);

armsig_mx = crosscov_thetatimes(HP_Theta_globalzero2,thisdir,false);


%Run  draft_thetaseq_again.m
[Mat,Index,binspike] = make_PosteriorfromThetaCandEvents(thisdir,ach);

rat = NaN(size(ach,1),3);
for c = 1:size(ach,1)
    Ind = 4*Index(c)+1:4*Index(c+1)-3;
    b = binspike(:,Ind);
    rat(c,2) = sum(sum(b,2)>0);
    if ~isnan(pc(c,2)) && rat(c,2)>4
        
        sumarm = NaN(3,1);
        for iarm = 1:3
            sumarm(iarm)= nansum(nansum(Mat(armp==iarm,Ind)./(sum(armp==iarm)/length(armp))));
        end
        sa = sumarm;
        sumarm(pc(c,2)) = NaN;
        [~,inda] = max(sumarm);
    %     [~,indi] = min(sumarm);
        rat(c,1) = (sumarm(inda)/sum(sa));    
        rat(c,3) = inda;
        
        midp = NaN(length(Ind),1);
        midp(1:floor(length(Ind)/2)) = 1;
        midp(ceil(length(Ind)/2)+1:end) = 0;
        midp2 = circshift(midp,round((shift*(length(Ind)/360))));        
        sumarm2 = NaN(3,2);
        for iarm = 1:3               
            sumarm2(iarm,1)= (nansum(nansum(Mat(armp==iarm,midp==1)./(sum(armp==iarm)/length(armp)))))./sum(midp==1);
            sumarm2(iarm,2)= (nansum(nansum(Mat(armp==iarm,Ind(midp==0))./(sum(armp==iarm)/length(armp)))))./sum(midp==0);
        end
        sumarm3 = sumarm2(inda,:)./sum(sumarm2,1);
        if length(unique(sumarm3))>1; [~,b] = max(sumarm3); else; b = NaN; end
        rat(c,4) = b;
        sumarm2 = NaN(3,2);
        for iarm = 1:3               
            sumarm2(iarm,1)= (nansum(nansum(Mat(armp==iarm,midp2==1)./(sum(armp==iarm)/length(armp)))))./sum(midp2==1);
            sumarm2(iarm,2)= (nansum(nansum(Mat(armp==iarm,Ind(midp2==0))./(sum(armp==iarm)/length(armp)))))./sum(midp2==0);
        end
        sumarm3 = sumarm2(inda,:)./sum(sumarm2,1);
        if length(unique(sumarm3))>1; [~,b2] = max(sumarm3); else; b = NaN; end
        rat(c,5) = b2;
        
        if rat(c,1)>.6  && exampleplots    
%             [~,~,i]= histcounts(ach(c,:),pos(:,1));
%             figure; subplot(2,1,1); plot(linposcat,'k.'); hold on;  plot(i(1):i(end),linposcat(i(1):i(end)),'r*')
%             subplot(2,1,2); imagesc(Mat(:,4*Index(c)+1:4*Index(c+1))); axis xy
%             disp('stop')
            c1 = find(achange(:,1)==ach(c,1));
            [~,~,i]= histcounts(mean(achange(c1-20:c1+20,:),2),pos(:,1));
            [Mat2,Index2,~] = make_PosteriorfromThetaCandEvents(thisdir,achange(c1-20:c1+20,:));
            
            figure; imagesc(Mat2); colormap hot; axis xy; hold on
            plot([4*Index2(21)+1 4*Index2(21)+1],[0 163],'w-','LineWidth',3)
            plot([4*Index2(22) 4*Index2(22)],[0 163],'w-','LineWidth',3)
            plot(Index2(1:end-1)*4,binpos(i),'w*')            
%             xt = get(gca,'xtick');
%             set(gca,'xticklabel',(xt*5)/1000)
            xtimes = achange(c1-20:c1+20,1)-achange(c1-20,1);
            set(gca,'xtick',Index2(2:4:end)*4,'xticklabel',xtimes(1:4:end))
            xlabel('Time (Seconds)')
            ylabel('Binned Position')
            set(gca,'FontSize',18)
            set(gcf,'Position',[ 1921         -79        1600        1083])
            helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\Example_' num2str(thisdir(1:end-4)) '_' num2str(c)])
        end
    end
end

ind = find(rat(:,1)>.4); %was .34, then .4
ttimes = ach(ind,:);
tprosp = rat(ind,3);
tarmon = pc(ind,2);
theadarm = pc(ind,3);
tlap = pc(ind,1);
times_armon_thetaof_headingarm_lap = [ttimes tarmon tprosp theadarm tlap rat(ind,4) rat(ind,5)]; % added which half of the theta cycle the representation was strongest
times_armon_thetaof_headingarm_lap_thetahalf_all = [ach pc(:,2) rat(:,3) pc(:,3) pc(:,1) rat(:,4) rat(:,5) rat(:,1)]; % added which half of the theta cycle the representation was strongest
% use this to then see if the theta sweeps predict the behavior

disp('Half Theta:')
[sum(rat(~isnan(rat(:,1)),4)==2)./sum(rat(~isnan(rat(:,1)),4)==1) sum(rat(rat(:,1)>.4,4)==2)./sum(rat(rat(:,1)>.4,4)==1)]
disp('Half ThetaSeq:')
[sum(rat(~isnan(rat(:,1)),5)==2)./sum(rat(~isnan(rat(:,1)),5)==1) sum(rat(rat(:,1)>.4,5)==2)./sum(rat(rat(:,1)>.4,5)==1)]


% thdiff = NaN(3,2);
% aa = times_armon_thetaof_headingarm_lap;
% for iarm = 1:3
%     otharms = setdiff(1:3,iarm);
%     for ihead = 1:2        
%         a = sum(aa(:,3)==iarm & aa(:,4)==otharms(ihead) & aa(:,5)==otharms(ihead));
%         b = sum(aa(:,3)==iarm & aa(:,4)==otharms(setdiff(1:2,ihead)) & aa(:,5)==otharms(ihead));
%         thdiff(iarm,ihead) = (a-b)./(a+b);
%     end
% end

load(thisdir,'spikedata','other_cells')



spks = spikedata(ismember(spikedata(:,2),other_cells),1:2);
spks2 = [];
pfcfr = NaN(length(ttimes),size(other_cells,1));
for it = 1:size(ttimes,1)
    spks2 = [spks2; spks(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2),:)...
        tprosp(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        tarmon(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        it*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        theadarm(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        tlap(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)];    
    h = histc(spks(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2),2),other_cells);
    pfcfr(it,:) = h./diff(ttimes(it,:),[],2);
end

pfc3 = NaN(size(other_cells,1),3);
for iarm = 1:3
    tm = sum(diff(ttimes(tprosp==iarm,:),[],2));    
    for icell = 1:length(other_cells)
        pfc3(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==iarm);
    end
    pfc3(:,iarm) = pfc3(:,iarm)./tm; 
end

pfc3b = NaN(size(other_cells,1),3);
for iarm = 1:3
    tm = sum(diff(ttimes(tarmon==iarm,:),[],2));    
    for icell = 1:length(other_cells)
        pfc3b(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,4)==iarm);
    end
    pfc3b(:,iarm) = pfc3b(:,iarm)./tm; 
end

pfc6 = NaN(size(other_cells,1),3,2);
pfc6alt = NaN(size(other_cells,1),3,2,size(ttimes,1));
for iarm = 1:3
    otharm = setdiff(1:3,iarm);
    for ioth = 1:2
        tm = sum(diff(ttimes(tprosp==otharm(ioth) & tarmon==iarm,:),[],2));    
        ttimes2 = find(tprosp==otharm(ioth) & tarmon==iarm);
        for icell = 1:length(other_cells)
            pfc6(icell,iarm,ioth) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm);
            for it = 1:length(ttimes2)
                pfc6alt(icell,iarm,ioth,ttimes2(it)) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm & spks2(:,5)==ttimes2(it))...
                    ./sum(diff(ttimes(ttimes2(it),:),[],2));
            end
        end
        pfc6(:,iarm,ioth) = pfc6(:,iarm,ioth)./tm;         
    end
end
pfc6v1 = pfc6;
pfc6 = reshape(pfc6,[size(pfc6,1) 6]);
prospectiveTheta_armcombs = pfc6;
prospectiveTheta_toplot = pfc6alt;
clear pfc6 pfc6alt

%for broken up by heading arm too

pfc6head = NaN(size(other_cells,1),3,2,2);

% pfc6althead = NaN(size(other_cells,1),3,2,max(tlap),2);
testtm = NaN(3,2,2);
for iarm = 1:3
    otharm = setdiff(1:3,iarm);
    for ioth = 1:2
        for ihead = 1:2
            tm = sum(diff(ttimes(tprosp==otharm(ioth) & tarmon==iarm & theadarm==otharm(ihead),:),[],2));    
%             ttimes2 = find(tprosp==otharm(ioth) & tarmon==iarm & theadarm==otharm(ihead));
            for icell = 1:length(other_cells)
                pfc6head(icell,iarm,ioth,ihead) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm & spks2(:,6)==otharm(ihead));
                
%                 ulaps = unique(tlap(ttimes2));
%                 for it = 1:length(ulaps)
%                     pfc6althead(icell,iarm,ioth,ulaps(it),ihead) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm & spks2(:,7)==ulaps(it) & spks2(:,6)==otharm(ihead))...
%                         ./sum(diff(ttimes(ttimes2(tlap(ttimes2)==ulaps(it)),:),[],2));
%                     %%%check if this makes sense, changed
%                 end
            end
            pfc6head(:,iarm,ioth,ihead) = pfc6head(:,iarm,ioth,ihead)./tm; 
            testtm(iarm,ioth,ihead) = tm;
        end
    end
end
pfc6headv1 = pfc6head;
pfc6head = reshape(pfc6head,[size(pfc6head,1) 6 2]); % check this
    
prospectiveTheta_armcombs_head = pfc6head;
% prospectiveTheta_toplot_head = pfcfr;
prsopectiveTheta_pfcfr = pfcfr;
%% again for outward


% rat = NaN(size(ach,1),3);
% for c = 1:size(ach,1)
%     b = binspike(:,4*Index(c)+1:4*Index(c+1));
%     rat(c,2) = sum(sum(b,2)>0);
%     if ~isnan(oc(c,2)) && rat(c,2)>4
%         sumarm = NaN(3,1);
%         for iarm = 1:3
%             sumarm(iarm)= sum(sum(Mat(armp==iarm,4*Index(c)+1:4*Index(c+1))./(sum(armp==iarm)/length(armp))));
%         end
%         sa = sumarm;
%         sumarm(oc(c,2)) = NaN;
%         [~,inda] = max(sumarm);
%     %     [~,indi] = min(sumarm);
%         rat(c,1) = (sumarm(inda)/sum(sa));    
%         rat(c,3) = inda;
%     end
% end
% 
% ind = find(rat(:,1)>.4); %was .34
% ttimes = ach(ind,:);
% tprosp = rat(ind,3);
% tarmon = oc(ind,2);
% 
% spks = spikedata(ismember(spikedata(:,2),other_cells),1:2);
% spks2 = [];
% for it = 1:length(ttimes)
%     spks2 = [spks2; spks(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2),:)...
%         tprosp(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
%         tarmon(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)];
% end
% 
% Opfc3 = NaN(size(other_cells,1),3);
% for iarm = 1:3
%     tm = sum(diff(ttimes(tprosp==iarm,:),[],2));    
%     for icell = 1:length(other_cells)
%         Opfc3(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==iarm);
%     end
%     Opfc3(:,iarm) = Opfc3(:,iarm)./tm; 
% end
% 
% Opfc3b = NaN(size(other_cells,1),3);
% for iarm = 1:3
%     tm = sum(diff(ttimes(tarmon==iarm,:),[],2));    
%     for icell = 1:length(other_cells)
%         Opfc3b(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,4)==iarm);
%     end
%     Opfc3b(:,iarm) = Opfc3b(:,iarm)./tm; 
% end
% 
% Opfc6 = NaN(size(other_cells,1),3,2);
% for iarm = 1:3
%     otharm = setdiff(1:3,iarm);
%     for ioth = 1:2
%         tm = sum(diff(ttimes(tprosp==otharm(ioth) & tarmon==iarm,:),[],2));    
%          for icell = 1:length(other_cells)
%             Opfc6(icell,iarm,ioth) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm);
%         end
%         Opfc6(:,iarm,ioth) = Opfc6(:,iarm,ioth)./tm; 
%     end
% end
% Opfc6 = reshape(Opfc6,[size(Opfc6,1) 6]);
    
%%

save(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','armsig_mx','prospectiveTheta_armcombs','prospectiveTheta_toplot','prospectiveTheta_armcombs_head','prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap','-append');


disp('Done with prospective_theta')


% load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all')
%     rat = times_armon_thetaof_headingarm_lap_thetahalf_all(times_armon_thetaof_headingarm_lap_thetahalf_all(:,4)>.6,7);
% test(id) = sum(rat==2)./sum(rat==1);