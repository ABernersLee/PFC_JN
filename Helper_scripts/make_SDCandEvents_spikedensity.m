function make_SDCandEvents_spikedensity(thisdir,std_cutoff)

Time_Bin_Size=0.001;
VelThresh = 5;

load(thisdir,'hpinterneurons','spikedata','vel','hp_cells','pos')

cellstouse = hp_cells(~ismember(hp_cells,hpinterneurons));
HPSpikes = spikedata(ismember(spikedata(:,2),cellstouse),1);

Time_Bins=min(HPSpikes):Time_Bin_Size:min(HPSpikes)+Time_Bin_Size*ceil(range(HPSpikes)/Time_Bin_Size);
[A,~]=histc(HPSpikes,Time_Bins');
Spike_Density=zeros(length(Time_Bins)-1,3);
Spike_Density(:,2)=Time_Bins(1:end-1)'+Time_Bin_Size/2;
Filter=fspecial('gaussian',[100 1],10); % ref Brad Pfeiffer 2013
Spike_Density(:,1)=filter2(Filter,A(1:end-1));
[~,I]=histc(Spike_Density(:,2),[min([pos(1,1) Time_Bins(1)]); pos(1:end-1,1)+diff(pos(:,1))/2 ; max([pos(end,1) Time_Bins(end)])]);
Spike_Density(:,3)=vel(I);

sMean=mean(Spike_Density(abs(Spike_Density(:,3))<VelThresh,1));
sPeak=sMean+std_cutoff*std(Spike_Density(abs(Spike_Density(:,3))<VelThresh,1));
Check=diff([0 ; Spike_Density(:,1)>sMean & abs(Spike_Density(:,3))<VelThresh ; 0]);
% Check=diff([0 ; Spike_Density(:,1)>sMean ; 0]);
clear VelTresh

DurationBound=[0.1 .5];
SD_CandEventTimes=[];
SS=find(Check==1);
EE=find(Check==-1);
target=find(EE-SS>DurationBound(1)/Time_Bin_Size & EE-SS<DurationBound(2)/Time_Bin_Size);
if ~isempty(target)
    CandSeq1=[Spike_Density(SS(target),2) Spike_Density(EE(target)-1,2)];
    delete=zeros(length(target),1);
    for i=1:length(target)
        if Spike_Density(SS(target(i)):EE(target(i))-1,1)<sPeak
            delete(i)=1;
        end    
    end
    target2=find(delete==1);
    if ~isempty(target2)
        CandSeq1(target2,:)=[];
        target(target2,:) = [];
    end
end
SD_CandEventTimes=[SD_CandEventTimes;CandSeq1]; 


SD_CandEventTimes = [SD_CandEventTimes NaN(size(SD_CandEventTimes,1),1)];
for ievent = 1:length(SD_CandEventTimes)
    [~,m] = max(Spike_Density(SS(target(ievent)):EE(target(ievent))-1,1));
    SD_CandEventTimes(ievent,3) = Spike_Density(SS(target(ievent))+m,2); %add the time of the peak of the spike density
end

disp(['Number of SD Candidate Events: ' num2str(size(SD_CandEventTimes,1))])

save(thisdir,'SD_CandEventTimes','-append')
disp('Done with make_SDCandEvents_spikedensity')