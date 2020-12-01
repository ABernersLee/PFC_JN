function make_PFC_replayeventtriggeredmat(thisdir,label)
%this doesnt work, cant use hist like this because the windows overlap
load(thisdir,[label '_replay_stnd'],'other_cells','spikedata',[label '_replay_singlebothjoint'])

eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
clear([label '_replay_stnd'])
p = other_cells; clear other_cells

%newly moved up here so list indicies are correct
Event(singlebothjoint==3) = [];

PFCspk = spikedata(ismember(spikedata(:,2),p),1:2); 
binsize = .002;
window = [-5 5];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

E = [Event+window(1) Event+window(2)];
B = cumsum([0 diff(E')]);
[~,I] = histc(PFCspk(:,1),sortrows(E(:)));


PFCreplayspikes_list = [PFCspk ones(size(PFCspk,1),1)];
for i = 1:2:max(I)
    PFCreplayspikes_list(I==i,1) = PFCspk(I==i,1)-E(ceil(i/2),1);
    PFCreplayspikes_list(I==i,3) = ones(sum(I==i),1)*ceil(i/2);      
    PFCspk(I==i,1) = PFCspk(I==i,1)-E(ceil(i/2),1)+B(ceil(i/2));    
end
PFCspk(~mod(I,2),:)=[];
PFCreplayspikes_list(~mod(I,2),:)=[];

rs = NaN(length(p),length(ind)*size(E,1));
for icell = 1:length(p)    
    rs(icell,:) = histc(PFCspk(PFCspk(:,2)==p(icell),1),B(1)+binsize/2:binsize:round(B(end))-binsize/2);  
end
PFCreplayspikes_binned = reshape(rs,[size(rs,1) length(ind) size(E,1)]);

eval([label '_PFCreplayspikes_list = PFCreplayspikes_list;'])
eval([label '_PFCreplayspikes_binned = PFCreplayspikes_binned;'])
save(thisdir,[label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],'-append')



