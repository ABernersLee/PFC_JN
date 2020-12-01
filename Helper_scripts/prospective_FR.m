function prospective_FR(thisdir,velcutoff)
disp('Done with prospective_FR')
%formerly prospective_coding_fields
load(thisdir,'other_cells','spikedata','vel','behave_change_log','laps_singlepass','linposcat','armpos','RP_CandEventTimes','pos')

nS = 10;
numbins = 1; %was 10 (changed 10/22/18) change back probably, then was 1
s = spikedata(ismember(spikedata(:,2),other_cells),:);
fr = zeros(numbins,max(laps_singlepass),length(other_cells));   frlast = fr;    
thearms = NaN(max(laps_singlepass),2);

excl = false(size(pos,1),1);
for icd = 1:size(RP_CandEventTimes(:,1))
    excl = excl | (pos(:,1)>=RP_CandEventTimes(icd,1) & pos(:,1)<=RP_CandEventTimes(icd,2));
end

for ilap = 1:max(laps_singlepass)
%    leave = find(behave_change_log(:,4) & laps_singlepass == ilap,1,'first'); %4 leave licking
%    thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
%    lapind = leave:find(armpos~=thisarm & laps_singlepass == ilap,1,'first')-1;
%    nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
%    
    %normal
   leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
   arrivemid = find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first')-1;
   lapind = leave:arrivemid; % 6 is arrive middle
   lapindlast = arrivemid:find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); 
   
%    lapind = find(laps_singlepass==ilap,1,'first'):find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first')-1; % 6 is arrive middle
   
    %testing
%    leave = find(behave_change_log(:,4) & laps_singlepass == ilap,1,'first'); %4 leave lick
%    lapind = find(laps_singlepass==ilap,1,'first'):leave;
   
   nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
   thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
%    lapind = leave:find(armpos(1:find(laps_singlepass == ilap,1,'last')==thisarm,1,'last')); %thought about changing to on 10/22


   lapind(vel(lapind)<velcutoff | excl(lapind)==1) = []; %exclude low velocity and ripples
   lapindlast(vel(lapind)<velcutoff | excl(lapind)==1) = [];
   
   lappos = linposcat(lapind);
   test = ((lappos-min(linposcat(armpos==thisarm)))/range(linposcat(armpos==thisarm)))*100;
   posbin = ceil(test./(100/numbins));
   lapposlast = linposcat(lapindlast);
   testlast = ((lapposlast-min(linposcat(armpos==nextarm)))/range(linposcat(armpos==nextarm)))*100;
   posbinlast = ceil(testlast./(100/numbins));
   
   
   for icell = 1:length(other_cells)
       for ipos = 1:max(posbin)
           fr(ipos,ilap,icell) = sum(s(:,2)==other_cells(icell) & ismember(s(:,3),lapind(posbin==ipos)))./(sum(posbin==ipos)/60);
           frlast(ipos,ilap,icell) = sum(s(:,2)==other_cells(icell) & ismember(s(:,3),lapindlast(posbinlast==ipos)))./(sum(posbinlast==ipos)/60);
       end
   end
   thearms(ilap,:) = [thisarm nextarm];
end

% thelastarms = [NaN NaN; thearms(2:end,1) thearms(1:end-1,1)];
maxarm = NaN(3,size(fr,3));
percdiff = maxarm;
perchalf = NaN(3,size(fr,3),2);
SSD = NaN(size(fr,3),3);

SSDf = NaN(size(fr,3),nS,3);
allarm = NaN(3,3,size(fr,3));
toplotfd = NaN(size(fr,3),size(fr,1),3,2);
toplotfd_shuff = NaN(size(fr,3),size(fr,1),3,2,nS);
toplotfdlast = toplotfd;
for iarm = 1:3
   otharms = setdiff(1:3,iarm);
   if length(unique(thearms(thearms(:,1)==iarm,2)))>1 
       alaps = find(thearms(:,1)==iarm & thearms(:,2)==otharms(1));
       blaps = find(thearms(:,1)==iarm & thearms(:,2)==otharms(2));
       
       
       armslabel = thearms(thearms(:,1)==iarm,2);
       numperm = (factorial(length(armslabel))/(factorial(length(armslabel)-sum(armslabel==otharms(1)))*factorial(sum(armslabel==otharms(1)))));
       if numperm<nS
            disp(['Skipping arm ' num2str(iarm) ' FR , only ' num2str(numperm) ' possible permutations'])
           continue
          
       end
% %        if length(alaps)>0 && length(blaps)>0
           aind = floor(length(alaps)/2);
           bind = floor(length(blaps)/2);
           am1 = nanmean(squeeze(nanmean(fr(:,alaps(1:aind),:),2)));
           bm1 = nanmean(squeeze(nanmean(fr(:,blaps(1:bind),:),2)));
           
           am2 = nanmean(squeeze(nanmean(fr(:,alaps(aind+1:end),:),2)));           
           bm2 = nanmean(squeeze(nanmean(fr(:,blaps(bind+1:end),:),2)));
           
           [a1,~] = nanmax([am1;bm1],[],1);
           [a12,~] = nanmin([am1;bm1],[],1);           
           perchalf(iarm,:,1) = a1-a12; %((a1-a12)./a1)*100;
           
           [a1,~] = nanmax([am2;bm2],[],1);
           [a12,~] = nanmin([am2;bm2],[],1);           
           perchalf(iarm,:,2) = a1-a12; %((a1-a12)./a1)*100;
%        end
           
       
       
%        armslabellast = thelastarms(thelastarms(:,1)==iarm,2);
       fr2 = fr(:,thearms(:,1)==iarm,:);
       
%        m1 = squeeze(mean(nanmax(fr2(:,armslabel==otharms(1),:)),2));
%        m2 = squeeze(mean(nanmax(fr2(:,armslabel==otharms(2),:)),2));
       m1 = squeeze(mean(nanmean(fr2(:,armslabel==otharms(1),:)),2));  % 10/16 1:24pm did max, changed back 10/22 11:18am
       m2 = squeeze(mean(nanmean(fr2(:,armslabel==otharms(2),:)),2));
       SSD(:,iarm) = (m1-m2).^2;
       SSD(m1==0 & m2==0,iarm) = NaN;
       
%        disp([num2str(length(armslabel)) ' to perm'])
%        if length(armslabel)>10
       jnk = NaN(nS,length(armslabel));
       for iN = 1:nS
           jnk1 = armslabel(randperm(length(armslabel)));
           while sum(jnk1~=armslabel)==0 || (iN>2 && sum(sum(jnk1'~=jnk(1:iN-1,:),2)==0)>0)
               jnk1 = armslabel(randperm(length(armslabel)));
               disp(['Prospective FR Same, nS = ' num2str(iN)])               
           end
           jnk(iN,:) = jnk1;
       end
%        else
%            jnk = perms(armslabel);
%        end
       
%        disp([num2str(size(jnk,1)) ' permutations'])
       for iN = 1:min([size(jnk,1) nS])
%            armp = armslabel(randperm(length(armslabel)));
%            while sum(armp~=armslabel)==0
%                armp = armslabel(randperm(length(armslabel)));
%                disp('armlabel same')
%            end
           armp = jnk(iN,:);
           m1a = squeeze(mean(nanmax(fr2(:,armp==otharms(1),:)),2));
           m2a = squeeze(mean(nanmax(fr2(:,armp==otharms(2),:)),2));
           SSDf(:,iN,iarm) = (m1a-m2a).^2;
           
           
           toplotfd_shuff(:,:,iarm,1,iN) = squeeze(nanmean(fr2(:,armp==otharms(1),:),2))';
           toplotfd_shuff(:,:,iarm,2,iN) = squeeze(nanmean(fr2(:,armp==otharms(2),:),2))';
       end
       
       toplotfd(:,:,iarm,1) = squeeze(nanmean(fr2(:,armslabel==otharms(1),:),2))';
       toplotfd(:,:,iarm,2) = squeeze(nanmean(fr2(:,armslabel==otharms(2),:),2))';
           
       
%        m2 = max(squeeze(max(fr(:,thearms(:,1)==iarm & thearms(:,2)==otharms(2),:),2)));
       [a,b] = max([nanmean(m1) nanmean(m2)],[],2);
       [a2,~] = min([nanmean(m1) nanmean(m2)],[],2);
       maxarm(iarm,:) = otharms(b);
%        a2(a2==0) = .001;
       percdiff(iarm,:) = a-a2; %((a-a2)./a)*100;
       allarm(iarm,otharms(1),:) = m1;
       allarm(iarm,otharms(2),:) = m2;
%        percdiff(iarm,:) = (a-a2);
   end
   if sum(thearms(:,2)==iarm)>0
       armslabellast = thearms(thearms(:,2)==iarm,1);
       fr2last = frlast(:,thearms(:,2)==iarm,:);
   
       toplotfdlast(:,:,iarm,1) = squeeze(nanmean(fr2last(:,armslabellast==otharms(1),:),2))';
       toplotfdlast(:,:,iarm,2) = squeeze(nanmean(fr2last(:,armslabellast==otharms(2),:),2))';
   end
       
end

%changed from on 10/11/18
% [~,ma] = nanmax(nanmean(allarm,1),[],2);
% maxheadingarm = squeeze(ma);
% [~,ma] = nanmax(nanmean(allarm,2),[],1);
% maxonarm = squeeze(ma);

%changed to on 10/11/18
[a,b] = max(abs(min(allarm,[],2)-max(allarm,[],2)));
b(a==0) = NaN;
maxonarm = squeeze(b);
maxheadingarm = NaN(size(maxonarm));
for icell = 1:size(allarm,3)
    if ~isnan(maxonarm(icell))
        [~,maxheadingarm(icell)] = max(allarm(maxonarm(icell),:,icell));     
    end
end
%

whicharm = [maxonarm maxheadingarm];

armsnum = 1:3;
tolook = ~(sum(isnan(SSD),1)==size(SSD,1));
armsnum = armsnum(tolook);
pSSD = NaN(size(fr,3),sum(tolook));
for iarm = 1:sum(tolook)
    pSSD(:,armsnum(iarm)) = (sum(SSD(:,armsnum(iarm))<SSDf(:,:,armsnum(iarm)),2)+1)./(nS+1);
end


toward1 = NaN(size(toplotfd,1),3); toward2 = toward1; 
toward4 = toward1;
armcombs = NaN(size(toplotfd,1),6); 
armcombs_Shuff = NaN(size(toplotfd,1),6,nS); 
armcombslast = armcombs;
for icell = 1:size(toplotfd,1)
    dat = squeeze(max(toplotfd(icell,:,:,:),[],2));
    
    datlast = squeeze(max(toplotfdlast(icell,:,:,:),[],2));
%     toward1(icell,1) = nanmean([dat(2,1)-dat(2,2) dat(3,1)-dat(3,2)]);
%     toward1(icell,2) = nanmean([dat(1,1)-dat(1,2) dat(3,2)-dat(3,1)]);
%     toward1(icell,3) = nanmean([dat(1,2)-dat(1,1) dat(2,2)-dat(2,1)]);
% %         toward1(icell,sum(isnan(dat),2)>0) = NaN;
% 
%     toward2(icell,1) = nanmean([dat(2,1) dat(3,1)]);
%     toward2(icell,2) = nanmean([dat(1,1) dat(3,2)]);
%     toward2(icell,3) = nanmean([dat(1,2) dat(2,2)]);

    armcombs(icell,:) = dat(:); %(dat(:)-nanmin(dat(:)))./(nanmax(dat(:))-nanmin(dat(:)));
    armcombslast(icell,:) = datlast(:);
    toward4(icell,:) = nanmean(dat,2);
    for ii = 1:nS
        datS = squeeze(max(toplotfd_shuff(icell,:,:,:,ii),[],2));
        armcombs_Shuff(icell,:,ii) = datS(:);
    end
end

prospectiveFR_percdiff = percdiff;
prospectiveFR_perchalf = perchalf;
prospectiveFR_whicharm = whicharm;
prospectiveFR_toplotfd = toplotfd;
prospectiveFR_armcombs = armcombs;
prospectiveFR_armcombs_shuff = armcombs_Shuff;
prospectiveFR_armcombslast = armcombslast;
prospectiveFR_meanFR = toward4;
save(thisdir,'prospectiveFR_armcombslast','prospectiveFR_armcombs_shuff','prospectiveFR_meanFR','prospectiveFR_percdiff','prospectiveFR_perchalf',...
    'prospectiveFR_whicharm','prospectiveFR_toplotfd','prospectiveFR_armcombs','-append')

disp('Done with prospective_FR')