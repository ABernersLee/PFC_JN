function make_PFC_candeventtriggeredmat(thisdir,label)
disp('Start make_PFC_candeventtriggeredmat')
load(thisdir,[label '_CandEventTimes'],'other_cells','spikedata')

eval(['Event = ' label '_CandEventTimes(:,1);'])
clear([label '_replay_stnd'])
p = other_cells; clear other_cells


PFCspk = spikedata(ismember(spikedata(:,2),p),1:2); 
nS = 500;
modwin = [0 .2];
basewin = [-.5 -.1];

binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
ind = [window(1):binsize:window(2)];

E = [Event+window(1) Event+window(2)];
PFCcandspikes_list = [];
PFCcandspikes_binned = [];
pfcind = NaN(max(p),1);
pfcind(p) = 1:length(p);
for iE = 1:size(E,1)
    Espikes1 = PFCspk(PFCspk(:,1)>=E(iE,1) & PFCspk(:,1)<=E(iE,2),:);
    [~,~,i] = histcounts(Espikes1(:,1),E(iE,1):binsize:E(iE,2));
    Espikes = [Espikes1 i Espikes1(:,1)-Event(iE,1)'];
    PFCcandspikes_list = cat(1,PFCcandspikes_list,[Espikes(:,[4 2]) ones(size(Espikes,1),1)*iE]);
    tmpind = sub2ind([length(p),length(ind)],pfcind(Espikes(:,2)),Espikes(:,3));
    tmp = false(length(p),length(ind));
    tmp(tmpind) = true;
    PFCcandspikes_binned = cat(3,PFCcandspikes_binned,tmp);    
end

ModShuff = NaN(length(p),size(Event,1),nS);
ModReal = NaN(length(p),size(Event,1));
BaseReal = ModReal;
PFCs = PFCcandspikes_list;
for iE = 1:size(E,1)    
    k = (rand(1,nS)-window(2)/2); %[[0 to 1] - 2] should be -2 to -1 ... should be from -1 to 1. ADDED THE /2 ON 9/17/2020
    j = repmat((PFCs(PFCs(:,3)==iE,1)),[1 nS])+ones(sum(PFCs(:,3)==iE),1)*k; % randomly adding -1 to -2, should be from -1 to 1
    j(j<-window(1)) = window(2)/2+j(j<-window(1)); 
    j(j>(window(2))) = window(2)/2-j(j>(window(2)));
    for icell = 1:length(p)
        ModShuff(icell,iE,:) = sum(j(PFCs(PFCs(:,3)==iE,2)==p(icell),:)>=modwin(1) & j(PFCs(PFCs(:,3)==iE,2)==p(icell),:)<=modwin(2),1);
        q = PFCs(PFCs(:,3)==iE & PFCs(:,2)==p(icell),1);
        ModReal(icell,iE) = sum(q>=modwin(1) & q<=modwin(2));
        BaseReal(icell,iE) = sum(q>=basewin(1) & q<=basewin(2));
    end        
end
meanShuf = squeeze(nanmean(ModShuff,2));
meanReal = nanmean(ModReal,2);
Rdiff = ((meanReal./range(modwin))-(nanmean(meanShuf,2)./range(modwin))).^2;
Sdiff = ((meanShuf./range(modwin))-(nanmean(meanShuf,2)./range(modwin))).^2;

Cand_p = (sum(Sdiff>=Rdiff,2)+1)/(nS+1); %corrected 2/7/19 (didnt have or=)

include = (sum(ModReal,2)+sum(BaseReal,2))>=50;
modu = ((nanmean(ModReal,2)./range(modwin)) - (nanmean(BaseReal,2)./range(basewin)))./ ...
    ((nanmean(ModReal,2)./range(modwin)) + (nanmean(BaseReal,2)./range(basewin)));
sig_modu_include = [Cand_p<.05 modu include];


eval([label '_PFCcandspikes_list = PFCcandspikes_list;'])
eval([label '_PFCcandspikes_binned = PFCcandspikes_binned;'])
eval([label '_Cand_sig_modu_include = sig_modu_include;'])
eval([label '_Cand_p = Cand_p;'])
save(thisdir,[label '_PFCcandspikes_binned'],[label '_PFCcandspikes_list'],[label '_Cand_sig_modu_include'],[label '_Cand_p'],'-append')
disp('Done with make_PFC_candeventtriggeredmat')