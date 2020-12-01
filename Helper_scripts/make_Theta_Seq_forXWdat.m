function [pfc3,pfc6,pfc3b,Opfc3,Opfc6,Opfc3b] = make_Theta_Seq_forXWdat(thisdir)

SpeedCutoff = 10;
load(thisdir,'HP_Theta')
theta = HP_Theta;

ind = find(theta(:,4)>358 & theta(:,4)<360);
ind(diff(ind)==1) = [];
anglechange = theta(ind,1);
as = [anglechange(1:end-1) anglechange(2:end)];
ind2 = (as(:,2)-as(:,1))<.1 | (as(:,2)-as(:,1))>.2;    
as(ind2,:) = [];  
achange = as;    


load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','linposcat')
[~,armp] = max(armposindex,[],2); clear armposindex
lappro = NaN(size(laps_singlepass,1),2);
lapout = NaN(size(laps_singlepass,1),2);
% lapdir = NaN(size(laps_singlepass,1),1);
for ilap = 1:max(laps_singlepass)
   leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
   lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle      
   lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform
   
   nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
   thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
   lappro(lapind,:) = [ilap*ones(length(lapind),1) thisarm*ones(length(lapind),1)];
   lapout(lapout1,:) = [ilap*ones(length(lapout1),1) nextarm*ones(length(lapout1),1)];
end
% lapall(~isnan(lapall(:,1)),2) = armpos(~isnan(lapall(:,1)));
% lapall = [lapall lapdir];
lappro(vel<SpeedCutoff | dirdat==1,:) = NaN;
lapout(vel<SpeedCutoff | dirdat==2,:) = NaN;
% lapall(vel<SpeedCutoff,:) = NaN;


[~,~,i2] = histcounts(achange(:,1),pos(:,1));
ind = [i2(1:end-1)+1 i2(2:end)];
ind = [ind; ind(end)+1 size(pos,1)];
achange(i2==0,:) = [];

prospectivecycles = NaN(size(achange,1),2);
outcycle = NaN(size(achange,1),2);
for itheta = 1:size(achange,1)
    prospectivecycles(itheta,:) = mode(lappro(ind(itheta,1):ind(itheta,2),:));
    outcycle(itheta,:) = mode(lapout(ind(itheta,1):ind(itheta,2),:));
end

torun = ~isnan(outcycle(:,1)) | ~isnan(prospectivecycles(:,1));
% torun = ~isnan(prospectivecycles(:,1));
ach = achange(torun,:);
oc = outcycle(torun,:);
pc = prospectivecycles(torun,:);
[Mat,Index,binspike] = make_PosteriorfromThetaCandEvents(thisdir,ach);

rat = NaN(size(ach,1),3);
for c = 1:size(ach,1)
    b = binspike(:,4*Index(c)+1:4*Index(c+1)-3);
    rat(c,2) = sum(sum(b,2)>0);
    if ~isnan(pc(c,2)) && rat(c,2)>4
        sumarm = NaN(3,1);
        for iarm = 1:3
            sumarm(iarm)= sum(sum(Mat(armp==iarm,4*Index(c)+1:4*Index(c+1)-3)./(sum(armp==iarm)/length(armp))));
        end
        sa = sumarm;
        sumarm(pc(c,2)) = NaN;
        [~,inda] = max(sumarm);
    %     [~,indi] = min(sumarm);
        rat(c,1) = (sumarm(inda)/sum(sa));    
        rat(c,3) = inda;
%         if rat(c,1)>.6          
%             [~,~,i]= histcounts(ach(c,:),pos(:,1));
%             figure; subplot(2,1,1); plot(linposcat,'k.'); hold on;  plot(i(1):i(end),linposcat(i(1):i(end)),'r*')
%             subplot(2,1,2); imagesc(Mat(:,4*Index(c)+1:4*Index(c+1))); axis xy
% %             disp('stop')
%         end
    end
end

ind = find(rat(:,1)>.4); %was .34
ttimes = ach(ind,:);
tprosp = rat(ind,3);
tarmon = pc(ind,2);

load(thisdir,'spikedata','other_cells')



spks = spikedata(ismember(spikedata(:,2),other_cells),1:2);
spks2 = [];
for it = 1:length(ttimes)
    spks2 = [spks2; spks(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2),:)...
        tprosp(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        tarmon(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)];
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
for iarm = 1:3
    otharm = setdiff(1:3,iarm);
    for ioth = 1:2
        tm = sum(diff(ttimes(tprosp==otharm(ioth) & tarmon==iarm,:),[],2));    
         for icell = 1:length(other_cells)
            pfc6(icell,iarm,ioth) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm);
        end
        pfc6(:,iarm,ioth) = pfc6(:,iarm,ioth)./tm; 
    end
end
pfc6 = reshape(pfc6,[size(pfc6,1) 6]);
    


% again for outward


rat = NaN(size(ach,1),3);
for c = 1:size(ach,1)
    b = binspike(:,4*Index(c)+1:4*Index(c+1)-3);
    rat(c,2) = sum(sum(b,2)>0);
    if ~isnan(oc(c,2)) && rat(c,2)>4
        sumarm = NaN(3,1);
        for iarm = 1:3
            sumarm(iarm)= sum(sum(Mat(armp==iarm,4*Index(c)+1:4*Index(c+1)-3)./(sum(armp==iarm)/length(armp))));
        end
        sa = sumarm;
        sumarm(oc(c,2)) = NaN;
        [~,inda] = max(sumarm);
    %     [~,indi] = min(sumarm);
        rat(c,1) = (sumarm(inda)/sum(sa));    
        rat(c,3) = inda;
    end
end

ind = find(rat(:,1)>.4); %was .34
ttimes = ach(ind,:);
tprosp = rat(ind,3);
tarmon = oc(ind,2);

spks = spikedata(ismember(spikedata(:,2),other_cells),1:2);
spks2 = [];
for it = 1:length(ttimes)
    spks2 = [spks2; spks(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2),:)...
        tprosp(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)...
        tarmon(it,1)*ones(sum(spks(:,1)>ttimes(it,1) & spks(:,1)<=ttimes(it,2)),1)];
end

Opfc3 = NaN(size(other_cells,1),3);
for iarm = 1:3
    tm = sum(diff(ttimes(tprosp==iarm,:),[],2));    
    for icell = 1:length(other_cells)
        Opfc3(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==iarm);
    end
    Opfc3(:,iarm) = Opfc3(:,iarm)./tm; 
end

Opfc3b = NaN(size(other_cells,1),3);
for iarm = 1:3
    tm = sum(diff(ttimes(tarmon==iarm,:),[],2));    
    for icell = 1:length(other_cells)
        Opfc3b(icell,iarm) = sum(spks2(:,2)==other_cells(icell) & spks2(:,4)==iarm);
    end
    Opfc3b(:,iarm) = Opfc3b(:,iarm)./tm; 
end

Opfc6 = NaN(size(other_cells,1),3,2);
for iarm = 1:3
    otharm = setdiff(1:3,iarm);
    for ioth = 1:2
        tm = sum(diff(ttimes(tprosp==otharm(ioth) & tarmon==iarm,:),[],2));    
         for icell = 1:length(other_cells)
            Opfc6(icell,iarm,ioth) = sum(spks2(:,2)==other_cells(icell) & spks2(:,3)==otharm(ioth) & spks2(:,4)==iarm);
        end
        Opfc6(:,iarm,ioth) = Opfc6(:,iarm,ioth)./tm; 
    end
end
Opfc6 = reshape(Opfc6,[size(Opfc6,1) 6]);
    






