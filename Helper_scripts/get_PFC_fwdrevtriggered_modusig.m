function get_PFC_fwdrevtriggered_modusig(thisdir,label)

load(thisdir,[label '_replay_singlebothjoint'],[label '_replay_corr'],[label '_replay_dirbias'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'other_cells','dirname')


eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Rcorr = ' label '_replay_corr;'])
eval(['dirbias = ' label '_replay_dirbias;'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])

pfc = other_cells; clear other_cells

%only single replays
PFCreplayspikes_binned(:,:,singlebothjoint~=1) = [];



% pos Rcorr is out, neg Rcorr is in ??? DOUBLE CHECK
%-1 bias is out, 1 is in
fwd = (Rcorr>0 & dirbias<0) | (Rcorr<0 & dirbias>0);
rev = (Rcorr>0 & dirbias>0) | (Rcorr<0 & dirbias<0);
include = NaN(size(PFCreplayspikes_binned,3),1);
include(fwd) = 1;
include(rev) = 2;
include1 = include(singlebothjoint==1);

PFCreplayspikes_binned(:,:,isnan(include1)) = [];
PFCreplayspikes_list(ismember(PFCreplayspikes_list(:,3),find(singlebothjoint~=1 | isnan(include))),:) = [];
include1(isnan(include1)) = [];
fwdrev = include1;
numcat   = 2;

nS = 500;
binsize = .002;
window = [-.5 .5];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

modind = ind>0 & ind<=.2;
baseind = ind>=-.5 & ind<-.1;

modu = NaN(numcat,length(pfc));
p_SSD = NaN(length(pfc),1);

for icell = 1:length(pfc)
    
    mm = NaN(numcat,sum(modind));
    mcount = NaN(numcat,1);
    for idir = 1:numcat        
        m1 = squeeze(sum(PFCreplayspikes_binned(icell,modind,fwdrev==idir),3));
        b1 = squeeze(sum(PFCreplayspikes_binned(icell,baseind,fwdrev==idir),3));
        m = nanmean(sum(PFCreplayspikes_binned(icell,modind,fwdrev==idir),2));
        b = nanmean(sum(PFCreplayspikes_binned(icell,baseind,fwdrev==idir),2)./2);
        modu(idir,icell) = (m-b)./(m+b);
        mm(idir,:) = (m1-(mean(b1)./(2*length(m1))))./sum(fwdrev==idir); %baseline subtract?
        mcount(idir) = sum(squeeze(sum(PFCreplayspikes_binned(icell,:,fwdrev==idir),3)));
    end
    
    if sum(mcount<50)==0
        
       SSD_obs = sum(mean(mm(1,:)-mm(2,:)).^2);        

        SSD_sh = NaN(nS,1);
        for i=1:nS
            replayarmSh = fwdrev(randperm(size(fwdrev,1)));
            ms = NaN(numcat,sum(modind));
            for idir = 1:numcat
                m = sum(PFCreplayspikes_binned(icell,modind,replayarmSh==idir),3);      
                b = nanmean(sum(PFCreplayspikes_binned(icell,baseind,replayarmSh==idir),3)./(2*length(m)));
                ms(idir,:) = (m-b)./sum(replayarmSh==idir);
            end
            
           SSD_sh(i,1) = sum(mean(ms(1,:)-ms(2,:)).^2);
            
        end

        p_SSD(icell,1) = (sum(SSD_sh>=SSD_obs)+1)/(nS+1);
    else
        p_SSD(icell,1) = NaN;
    end
      
end

% eval([label '_pSSDfwdrev = p_SSD'])
% eval([label '_modufwdrev = modu'])
% save(thisdir,[label '_pSSDfwdrev'],[label '_modufwdrev'],'-append')
disp('wait')
