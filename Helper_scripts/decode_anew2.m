function [Mat,Index] = decode_anew2(cellstouse,spikedata,InFR,OutFR,Cand)
% not using anymore
EstBin = .02;
Cell_Number = cellstouse;
Spike = [spikedata(:,2) spikedata(:,1)];
% InFR2 = repmat(InFR,[1 1 size(cellstouse,1)]);
% OutFR2 = repmat(OutFR,[1 1 size(cellstouse,1)]);
% InFR = InFR(Cell_Number,:);
% OutFR = OutFR(Cell_Number,:);
FR = InFR+OutFR;
FR3 = NaN(size(Cell_Number,1),size(FR,2),size(Cell_Number,2));
% OutFR3 = NaN(size(Cell_Number,1),size(OutFR,2),size(Cell_Number,2));
for icell = 1:size(Cell_Number,2)
    FR2 = FR(Cell_Number(:,icell),:);
    FR2(sum(FR2~=0,2)==0,:) = 1;
    FR2(FR2==0) = .0001;
    FR3(:,:,icell) = FR2;
end

Cand = [Cand(:,1) Cand(:,2)+.005];
% CandSeq = Cand;
% N=round(diff(Cand')/EstBin)';
% t=find(mod(N,2));
% Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
% t=find(~mod(N,2));
% Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[~,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);

% Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
end
Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
% clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;

% times = [times; times(end)+.005]; 
% Ind = Spike(:,2)>=times(1) & Spike(:,2)<=times(end);
% Spike(Ind,2)=Spike(Ind,2)-times(1);
% Spike(~Ind,:)=[]; % only takes spikes within the candidate events
% Start_Time=0;
% End_Time=range(times);
% TimeBins=ceil((End_Time-Start_Time)/EstBin)*4;

% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(size(Cell_Number,1),TimeBins,size(Cell_Number,2));
for i=1:4
    for icell = 1:size(Cell_Number,2)   
        for CellID=1:size(Cell_Number,1)   
            c=histc(Spike(Spike(:,1)==Cell_Number(CellID,icell),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
            binspike(CellID,i:4:TimeBins,icell)=c(1:TimeBins/4);
        end                      
    end
end
% TimeBins=round((End_Time-Start_Time)/EstBin*4);
% binspike = binspike(:,1:TimeBins);

clear Start_Time End_Time CellID c i 
clearvars -except TimeBins binspike FR3 EstBin B


Number_Of_Bins=size(FR3,2);
term2=zeros(Number_Of_Bins,TimeBins,size(binspike,3));
% term4=term2;

if Number_Of_Bins>TimeBins    
    for TBin=1:TimeBins
        term2(:,TBin,:)=prod(FR3.^binspike(:,TBin,:),1);    
%         term4(:,TBin,:)=prod(InFR.^binspike(:,TBin,:),1);
    end
else
    for PosBin=1:Number_Of_Bins   
        term2(PosBin,:,:)=prod((squeeze(FR3(:,PosBin,:))*ones(size(binspike,3),TimeBins)).^binspike,1);    
%         term4(PosBin,:,:)=prod((squeeze(InFR3(:,PosBin,:))*ones(size(binspike,3),TimeBins)).^binspike,1);            
    end
end

term3=exp(-EstBin*squeeze(sum(FR3,1)))*ones(size(binspike,3),TimeBins);

% term5=exp(-EstBin*sum(InFR,1)')*ones(1,TimeBins);
% Mat=term2.*term3+term4.*term5;

Mat=term2.*term3;
clear term2 term3

Index=round(B/(EstBin));