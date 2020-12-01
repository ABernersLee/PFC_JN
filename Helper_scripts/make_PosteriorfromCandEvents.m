function make_PosteriorfromCandEvents(thisdir,label)
disp(['Start make_PosteriorfromCandEvents ' label])
EstBin=0.02;
load(thisdir,'InFR','OutFR','hp_cells','hpinterneurons','spikedata',[label '_CandEventTimes'])

% load('Ripple_Events')
% Cand = [Ripple_Events(:,1)-.02 Ripple_Events(:,2)+.02];

Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
Spike = [spikedata(:,2) spikedata(:,1)];
InFR = InFR(Cell_Number,:);
OutFR = OutFR(Cell_Number,:);

eval(['Cand= ' label  '_CandEventTimes(:,1:2);'])
Cand = [Cand(:,1)-.02 Cand(:,2)+.02];
CandSeq = Cand;
N=round(diff(Cand')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[~,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);

Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
end
Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
% clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for i=1:4
    for CellID=1:length(Cell_Number)   
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
    end                      
end

clear Start_Time End_Time CellID c i 

InFR(sum(InFR~=0,2)==0,:) = 1;
OutFR(sum(OutFR~=0,2)==0,:) = 1;
    
InFR(InFR==0) = .0001;
OutFR(OutFR==0) = .0001;


Number_Of_Bins=size(OutFR,2);
term2=zeros(Number_Of_Bins,TimeBins);
term4=term2;

if Number_Of_Bins>TimeBins    
    for TBin=1:TimeBins
        term2(:,TBin)=prod(OutFR.^binspike(:,TBin),1);    
        term4(:,TBin)=prod(InFR.^binspike(:,TBin),1);
    end
else
    for PosBin=1:Number_Of_Bins   
        term2(PosBin,:)=prod((OutFR(:,PosBin)*ones(1,TimeBins)).^binspike,1);    
        term4(PosBin,:)=prod((InFR(:,PosBin)*ones(1,TimeBins)).^binspike,1);
    end
end

term3=exp(-EstBin*sum(OutFR,1)')*ones(1,TimeBins);

term5=exp(-EstBin*sum(InFR,1)')*ones(1,TimeBins);

Mat=term2.*term3+term4.*term5;

OutMat1=term2.*term3;
OutMatrix=(OutMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';

InMat1=term4.*term5;
InMatrix=(InMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';    

Index=round(B/(EstBin));

eval([label  '_InMatrix =InMatrix;'])
eval([label  '_OutMatrix =OutMatrix;'])
eval([label  '_Index =Index;'])
save(thisdir,[label '_InMatrix'],[label '_OutMatrix'],[label '_Index'],'-append')
disp(['Done with make_PosteriorfromCandEvents ' label])