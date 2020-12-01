function make_RPCandEvents_ripple(thisdir,rip_thresh)

%rip_thresh is how many std above the mean
disp('Start make_RPCandEvents_ripple')
VelThresh = 5;
maxlength = 1; 
returnto = 1; 
minlength = .05; 
minsep = .05;

load(thisdir,'HP_Ripple','pos','vel')

% Index=cumsum(hist(pos(:,1),HP_Ripple(:,1)));
% Index(Index<=0)=1;
% Index(Index>size(pos,1))=size(pos,1);
HP_Ripple(HP_Ripple(:,1)<min(pos(:,1)) | HP_Ripple(:,1)>max(pos(:,1)),:) = [];
[~,~,Index] = histcounts(pos(:,1),HP_Ripple(:,1));
[~,~,Index2] = histcounts(HP_Ripple(:,1),pos(:,1));
k = find(Index<=0);
Index(k(k<length(Index)/2))=1;
Index(k(k>length(Index)/2))=max(Index);
k = find(Index2<=0);
Index2(k(k<length(Index2)/2))=1;
Index2(k(k>length(Index2)/2))=max(Index2);



meanRip = mean(HP_Ripple(Index(vel<VelThresh),2));
stdRip = std(HP_Ripple(Index(vel<VelThresh),2));

[~,Locations]=findpeaks(HP_Ripple(:,2),'MinPeakHeight',meanRip+(stdRip*rip_thresh));
% [~,Locations]=findpeaks(HP_Ripple(Index(vel<VelThresh),2),'MinPeakHeight',meanRip+(stdRip*rip_thresh));

% This eliminates peaks that occurred when the rat was moving
Locations=Locations(vel(Index2(Locations))<VelThresh);

% This finds the start and end timepoints for each event.
Ripple_Events=zeros(length(Locations),6);
Ripple_Events(:,3)=HP_Ripple(Locations,1);
for N=1:size(Ripple_Events,1)
    L=Locations(N);
    Ripple_Events(N,6) = L;
    while HP_Ripple(L,2)>(meanRip+(stdRip*returnto)) && L>1 %this finds the closest timepoint prior to the current peak that crosses the mean
        L=L-1;
    end
    Ripple_Events(N,1)=HP_Ripple(L,1);
    Ripple_Events(N,4) = L;
    L=Locations(N);
    while HP_Ripple(L,2)>(meanRip+(stdRip*returnto)) && L<size(HP_Ripple,1) %this finds the closest timepoint after the current peak that crosses the mean
        L=L+1;
    end
    Ripple_Events(N,2)=HP_Ripple(L,1);
    Ripple_Events(N,5) = L;
end

%get rid of weirdly long ones
Ripple_Events((Ripple_Events(:,2)-Ripple_Events(:,1))>maxlength,:) = [];

%combine close ones
for ind = 3:-1:1
    while sum(diff(Ripple_Events(:,ind))<minsep)>0 %ABL added
        merge = [diff(Ripple_Events(:,ind))<minsep;false];
        Ripple_Events(merge,[2 5]) = Ripple_Events(find(merge)+1,[2 5]);
        [~,tallest] = max([HP_Ripple(Ripple_Events(merge,6),2),HP_Ripple(Ripple_Events(find(merge)+1,6),2)],[],2);
        Ripple_Events(merge,[3 6]) = Ripple_Events(find(merge)+tallest-1,[3 6]);
        Ripple_Events(find(merge)+1,:) = [];
    end
end
% Ripple_Events=Ripple_Events(find(Ripple_Events(:,1)>0),:); shouldn't be a thing anymore

%get rid of short ones that couldn't be combined even to make longer ones
Ripple_Events((Ripple_Events(:,2)-Ripple_Events(:,1))<minlength,:) = [];

RP_CandEventTimes = Ripple_Events;

disp(['Number of RP Candidate Events: ' num2str(size(RP_CandEventTimes,1))])

save(thisdir,'RP_CandEventTimes','-append')
disp('Done with make_RPCandEvents_ripple')