function make_slowgama_ofsteps(thisdir,label)

load(thisdir,[label  '_replay_seqtimes'],[label  '_replay_stnd'],'slowgamma_localspikephase','other_cells','hp_cells','hpinterneurons')

eval(['seqtimes = ' label  '_replay_seqtimes;'])
eval(['Cand = ' label  '_replay_stnd;'])
clear([label  '_seqtimes_btwn'])
clear([label  '_seqindex_btwn;'])
clear([label  '_replay_stnd;'])

gspk = slowgamma_localspikephase(ismember(slowgamma_localspikephase(:,2),hp_cells(~ismember(hp_cells,hpinterneurons))),:);
PFCspk = slowgamma_localspikephase(ismember(slowgamma_localspikephase(:,2),other_cells),:);


%HP spikes to get slow gamma and spikes
[~,I] = histc(gspk(:,1),sortrows(Cand(:)));
gspk(~mod(I,2),:)=[];
clear I

[so,o] = sortrows(seqtimes);
[~,I] = histc(gspk(:,1),so);
gspk(I==0,:) = [];
I(I==0) = [];

%PFC spikes
[~,I2] = histc(PFCspk(:,1),sortrows(Cand(:)));
PFCspk(~mod(I2,2),:)=[];
clear I2

[so2,o2] = sortrows(seqtimes);
[~,I2] = histc(PFCspk(:,1),so2);
I2(I2==0) = [];
Is2 = unique(I2);

%get for each seq transition bin
Is = unique(I);
seqphases = NaN(size(seqtimes));
seqHPspikes = zeros(size(seqtimes));
seqPFCspikes = seqHPspikes;
for ii = 1:length(Is)
    seqHPspikes(o(Is(ii))) = sum(I==Is(ii));    
    alpha = deg2rad(gspk(I==Is(ii),5));
    r = sum(exp(1i*alpha),1);
    ph = mod(((angle(r)/pi)*180)+360,360);
    seqphases(o(Is(ii))) = ph;
end

for ii = 1:length(Is2)
    seqPFCspikes(o2(Is2(ii))) = sum(I2==Is2(ii));
end
    
eval([label  '_seqPFCspikes = seqPFCspikes;'])
eval([label  '_seqHPspikes = seqHPspikes;'])
eval([label  '_seqphases = seqphases;'])

save(thisdir,[label '_seqPFCspikes'],[label  '_seqHPspikes'],[label  '_seqphases'],'-append')

