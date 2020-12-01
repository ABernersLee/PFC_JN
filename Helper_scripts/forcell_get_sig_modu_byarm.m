function [modu,p_SSD,SSD_obs1,SSD_obs2] = forcell_get_sig_modu_byarm(icell,modind1,baseind1,replayarm,PFCreplayspikes_binned)



binsize = .02;
window = [-2 2];
% window = [-.5 .5];
ind = [window(1):binsize:window(2)];

modind = ind>modind1(1) & ind<=(modind1(2));
baseind = ind>=baseind1(1) & ind<baseind1(2);

modu = NaN(3,1);
nS = 500;
mm = NaN(3,sum(modind)); mm2 = mm;
mcount = NaN(3,1);
for iarm = 1:3        
    m1 = squeeze(sum(PFCreplayspikes_binned(icell,modind,replayarm==iarm),3))./range(modind1);
    b1 = squeeze(sum(PFCreplayspikes_binned(icell,baseind,replayarm==iarm),3))./range(baseind1);
    m = squeeze(squeeze(sum(PFCreplayspikes_binned(icell,modind,replayarm==iarm),2)./range(modind1)));
    b = squeeze(sum(PFCreplayspikes_binned(icell,baseind,replayarm==iarm),2)./range(baseind1));
    modu(iarm) = nanmean((m-b)./(m+b));
    mm(iarm,:) = (m1-mean(b1))./(m1+mean(b1));
    mm2(iarm,:) = m1-mean(b1);
%     mm(iarm,:) = (m1-(mean(b1)))./sum(replayarm==iarm); %baseline subtract?
    mcount(iarm) = sum(squeeze(sum(PFCreplayspikes_binned(icell,:,replayarm==iarm),3)));
end

combos = nchoosek(1:3,2);
if sum(sum(mcount))>50
    SSD_obs1 = 0;
    SSD_obs2 = 0;
    for ic = 1:size(combos,1)
       SSD_obs1 = SSD_obs1+sum((mm(combos(ic,1),:)-mm(combos(ic,2),:)).^2);
       SSD_obs2 = SSD_obs2+sum((mm2(combos(ic,1),:)-mm2(combos(ic,2),:)).^2);
    end

    SSD_sh = zeros(nS,1);
    for i=1:nS
        replayarmSh = replayarm(randperm(size(replayarm,1)));
        ms = NaN(3,sum(modind));
        for iarm = 1:3
                m = sum(PFCreplayspikes_binned(icell,modind,replayarmSh==iarm),3)./range(modind1);      
                b = nanmean(sum(PFCreplayspikes_binned(icell,baseind,replayarmSh==iarm),3)./range(baseind1));
                ms(iarm,:) = (m-b)./(m+b);
        end
        for ic = 1:size(combos,1)
           SSD_sh(i,1) = SSD_sh(i,1)+sum((ms(combos(ic,1),:)-ms(combos(ic,2),:)).^2);
        end

    end
    
    p_SSD = (sum(SSD_sh>=SSD_obs1)+1)/(nS+1);
else
    p_SSD = NaN; SSD_obs1 = NaN; SSD_obs2 = NaN;
end
      