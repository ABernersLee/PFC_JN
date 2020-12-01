function [single,joint] = get_PFC_armtriggered_modusig(thisdir,label,ilab)
disp('Start get_PFC_armtriggered_modusig')
load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_replayarm'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'other_cells',[label  '_replay_shuffle_p']...
    ,[label  '_replay_corr'],[label  '_replay_maxjump'],[label  '_replay_armcov'])


eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['corr = ' label '_replay_corr;'])
eval(['armcov1 = ' label  '_replay_armcov;'])

eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
eval(['mj = ' label  '_replay_maxjump;'])
eval(['wc = ' label  '_replay_corr;'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])

pfc = other_cells; clear other_cells

%all single replays %doing this earlier now
% PFCreplayspikes_binned(:,:,singlebothjoint==3) = [];
% PFCreplayspikes_list(ismember(PFCreplayspikes_list(:,3),find(singlebothjoint==3)),:) = [];

% replayarm = replayarm(singlebothjoint~=3 & shuffle_p<.05); 
% replayarm = replayarm(singlebothjoint~=3); %changed 2/4/19

% replaylab = {'AllArmEvents';'All.3Events';'AllSigEvents';'All.3SigEvents';'All.5SigEvents';'All.5Events'}; 'All.3SigEvents_JD.6';'All.5SigEvents';'All.6SigEvents_JD.4'};
% replaylab = {'AllArmEvents';'All.3Events_armcov3';'AllSigEvents_ArmCov3';'All.3SigEventsArmcov5';'All.5SigEventsArmcov5'}; 
% replaylab = {'AllArmEvents';'Wc3ArmCov3';'WC4ArmCov5mj7';'WC5ArmCov5mj4'};
if ilab==1    
    replayarm(singlebothjoint==3) = NaN;
elseif ilab ==2
    replayarm(singlebothjoint==3 | abs(wc)<.3 | armcov1<.3) = NaN; %singlebothjoint==0) = NaN;   
elseif ilab == 3
    replayarm(singlebothjoint==3 | abs(wc)<.4 | armcov1<.5 | mj>.7) = NaN;
%     replayarm(singlebothjoint==3 | shuffle_p>=.05 | armcov1<.3) = NaN;  
elseif ilab == 4
    replayarm(singlebothjoint==3 | abs(wc)<.5 | armcov1<.5 | mj>.4) = NaN;
%     replayarm(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0) = NaN;    
elseif ilab == 5 
    replayarm(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0 | abs(wc)<.5) = NaN;    
elseif ilab == 6
    replayarm(singlebothjoint==3 | singlebothjoint==0 | abs(wc)<.5) = NaN;    
end

single = sum(singlebothjoint~=3 & shuffle_p<.05 | armcov1>=.5);
joint = sum(singlebothjoint==3);
% replayarm(singlebothjoint==3 | singlebothjoint==0) = NaN; 
% replayarm(singlebothjoint==3) = NaN; 

nS = 2000;
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
ind = [window(1):binsize:window(2)];
modind1 = [0 .2];
baseind1 = [-.5 -.1];

modind = ind>=modind1(1) & ind<=(modind1(2));
baseind = ind>=baseind1(1) & ind<=baseind1(2);

modu = NaN(3,length(pfc));
modu2 = modu;
p_SSD = NaN(length(pfc),1);
SSD_obs1 = p_SSD; SSD_obs2 = p_SSD;
combos = nchoosek(1:3,2);

IND = ~isnan(replayarm);
ra = replayarm(IND);
replayarmShSave = replayarm;
for i=1:nS
    replayarmSh = replayarm;
    replayarmSh(IND) = ra(randperm(size(ra,1)));
    while any(sum(replayarmShSave(~isnan(replayarm),:)~=replayarmSh(~isnan(replayarm)))==0) %sum(replayarmSh(~isnan(replayarmSh))~=replayarm(~isnan(replayarm)))==0
        replayarmSh = replayarm(randperm(size(replayarm,1)));
        disp(['triggeredsame i = ' num2str(i)])
    end
    replayarmShSave = cat(2,replayarmShSave,replayarmSh);
end
replayarmShSave = replayarmShSave(:,2:end);
        
for icell = 1:length(pfc)
    
    mm = NaN(3,1);
    mcount = NaN(3,1);
    for iarm = 1:3        
%        m1 = squeeze(sum(PFCreplayspikes_binned(icell,modind,replayarm==iarm),2))./range(modind1);
%         b1 = squeeze(sum(PFCreplayspikes_binned(icell,baseind,replayarm==iarm),2))./range(baseind1);
        m = squeeze(squeeze(sum(PFCreplayspikes_binned(icell,modind,replayarm==iarm),2)./range(modind1)));
        b = squeeze(sum(PFCreplayspikes_binned(icell,baseind,replayarm==iarm),2)./range(baseind1));
        modu(iarm,icell) = nanmean(m-b);
        modu2(iarm,icell) = (nanmean(m)-nanmean(b))/(nanmean(m)+nanmean(b));
        mm(iarm,:) = nanmean(m-b); %nanmean(m-b)/nanmean(m+b); %nanmean(m-b);
%         mm(iarm,:) = (m-b)./(m+b);
    %     mm(iarm,:) = (m1-(mean(b1)))./sum(replayarm==iarm); %baseline subtract?
        mcount(iarm) = sum(squeeze(sum(PFCreplayspikes_binned(icell,modind | baseind,replayarm==iarm),3)));
    end
    
    if sum(mcount)>50  %&& sum(mcount>50)>1
        SSD_obs = 0;
        SSD_obs2a = 0;
        for ic = 1:size(combos,1)
           SSD_obs = SSD_obs+sum((mm(combos(ic,1),:)-mm(combos(ic,2),:)).^2);
           SSD_obs2a = SSD_obs2a+sum((modu2(combos(ic,1),icell)-modu2(combos(ic,2),icell)).^2);
        end

        SSD_sh = zeros(nS,1);
        
       
        for i = 1:nS
            ms = NaN(3,1);
            for iarm = 1:3
                m = sum(PFCreplayspikes_binned(icell,modind,replayarmShSave(:,i)==iarm),2)./range(modind1);      
                b = nanmean(sum(PFCreplayspikes_binned(icell,baseind,replayarmShSave(:,i)==iarm),2)./range(baseind1));
%                 ms(iarm,:) = (m-b)./(m+b);
                ms(iarm,:) = nanmean(m-b); %nanmean(m-b)/nanmean(m+b); %nanmean(m-b);
            end
            for ic = 1:size(combos,1)
               SSD_sh(i,1) = SSD_sh(i,1)+sum((ms(combos(ic,1),:)-ms(combos(ic,2),:)).^2);
            end
        end

        p_SSD(icell,1) = (sum(SSD_sh>=SSD_obs)+1)/(nS+1);
        SSD_obs1(icell,1) = SSD_obs;
        SSD_obs2(icell,1) = SSD_obs2a;
    else
        p_SSD(icell,1) = NaN; SSD_obs1(icell,1) = NaN;
    end
      
end
% if ilab >1
%     label = [label '_' num2str(ilab)];
% end
eval([label '_pSSDarm = p_SSD;'])
eval([label '_moduarm = modu;'])
eval([label '_moduarm2 = modu2;'])
eval([label '_SSDarm = SSD_obs1;'])
eval([label '_SSDarm2 = SSD_obs2;'])
save(thisdir,[label '_pSSDarm'],[label '_moduarm'],[label '_moduarm2'],[label '_SSDarm'],[label '_SSDarm2'],'-append')
disp('Done with get_PFC_armtriggered_modusig')
