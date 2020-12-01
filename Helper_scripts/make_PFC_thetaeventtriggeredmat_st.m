function make_PFC_thetaeventtriggeredmat_st(thisdir)

load(thisdir,'spikedata','times_armon_thetaof_headingarm_lap_thetahalf_all','spikedata','other_cells','pos');
p = other_cells; clear other_cells
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
clear times_armon_thetaof_headingarm_lap_thetahalf_all
Event = Th(:,1);
PFCspk = spikedata(ismember(spikedata(:,2),p),1:2); 
binsize = .02;
window = [-2 2];
% ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

ind = [window(1):binsize:window(2)];
E = [Event+window(1) Event+window(2)];

PFCthetaspikes_list = [];
PFCthetaspikes_binned = [];
pfcind = NaN(max(p),1);
pfcind(p) = 1:length(p);
for iE = 1:size(E,1)
    Espikes1 = PFCspk(PFCspk(:,1)>=E(iE,1) & PFCspk(:,1)<=E(iE,2),:);
    [~,~,i] = histcounts(Espikes1(:,1),E(iE,1):binsize:E(iE,2));
    Espikes = [Espikes1 i Espikes1(:,1)-Event(iE,1)'];
    PFCthetaspikes_list = cat(1,PFCthetaspikes_list,[Espikes(:,[4 2]) ones(size(Espikes,1),1)*iE]);
    tmpind = sub2ind([length(p),length(ind)],pfcind(Espikes(:,2)),Espikes(:,3));
    tmp = false(length(p),length(ind));
    tmp(tmpind) = true;
    PFCthetaspikes_binned = cat(3,PFCthetaspikes_binned,tmp);    
end

PFCthetaspikes_binned_st = PFCthetaspikes_binned;
PFCthetaspikes_list_st = PFCthetaspikes_list;
save(thisdir,'PFCthetaspikes_binned_st','PFCthetaspikes_list_st','-append')

disp('Done with make_PFC_thetaeventtriggeredmat_st')
