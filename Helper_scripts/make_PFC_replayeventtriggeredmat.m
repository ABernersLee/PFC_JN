function make_PFC_replayeventtriggeredmat(thisdir,label)
disp('Start make_PFC_replayeventtriggeredmat')
load(thisdir,[label '_replay_stnd'],'other_cells','spikedata',[label '_replay_singlebothjoint'],[label  '_replay_shuffle_p']...
    ,[label  '_replay_corr'],[label  '_replay_maxjump'])

eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
eval(['mj = ' label  '_replay_maxjump;'])
eval(['wc = ' label  '_replay_corr;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
clear([label '_replay_stnd'])
p = other_cells; clear other_cells

%newly moved up here so list indicies are correct
% Event(singlebothjoint==3) = [];
% Event(singlebothjoint==3 | shuffle_p>=.05) = []; %changed 2/6/19
% Event(singlebothjoint==3) = []; %changed 2/4/19

PFCspk = spikedata(ismember(spikedata(:,2),p),1:2); 
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

ind = [window(1):binsize:window(2)];
E = [Event+window(1) Event+window(2)];

PFCreplayspikes_list = [];
PFCreplayspikes_binned = [];
pfcind = NaN(max(p),1);
pfcind(p) = 1:length(p);
for iE = 1:size(E,1)
    Espikes1 = PFCspk(PFCspk(:,1)>=E(iE,1) & PFCspk(:,1)<=E(iE,2),:);
    [~,~,i] = histcounts(Espikes1(:,1),E(iE,1):binsize:E(iE,2));
    Espikes = [Espikes1 i Espikes1(:,1)-Event(iE,1)'];
    PFCreplayspikes_list = cat(1,PFCreplayspikes_list,[Espikes(:,[4 2]) ones(size(Espikes,1),1)*iE]);
    tmpind = sub2ind([length(p),length(ind)],pfcind(Espikes(:,2)),Espikes(:,3));
    tmp = false(length(p),length(ind));
    tmp(tmpind) = true;
    PFCreplayspikes_binned = cat(3,PFCreplayspikes_binned,tmp);    
end

eval([label '_PFCreplayspikes_list = PFCreplayspikes_list;'])
eval([label '_PFCreplayspikes_binned = PFCreplayspikes_binned;'])
save(thisdir,[label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'-append')

disp('Done with make_PFC_replayeventtriggeredmat')


