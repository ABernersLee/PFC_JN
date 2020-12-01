function [OpenInFR,OpenOutFR] = make_directionalfields_open(thisdir)

VelThresh = 10;
Bin_Size = 3;
SigmaField = 5;

load(thisdir, 'pos','spikedata','vel','dirdat','hp_cells','other_cells')
Cell_Number = [hp_cells;other_cells];

P1 =[pos(:,1) ceil([pos(:,2)-min(pos(:,2))+.001 pos(:,3)-min(pos(:,3))+.001]./Bin_Size)];
S1=[spikedata(:,2) P1(spikedata(:,3),2) P1(spikedata(:,3),3) vel(spikedata(:,3)) dirdat(spikedata(:,3))];
D2 = dirdat; D2(vel<VelThresh) = [];
P2 = P1; P2(vel<VelThresh,:) = [];
S2 = S1; S2(S1(:,4)<VelThresh,:) = [];

OpenOutFR=zeros(max(P1(:,2)),max(P1(:,3)),max(Cell_Number));
OpenInFR=zeros(max(P1(:,2)),max(P1(:,3)),max(Cell_Number));
Two_D_Filter=fspecial('gaussian',[20 20],SigmaField);

for idir = 1:2    
    P = P2(D2==idir-1,:);
    S = S2(S2(:,5)==idir-1,:);
    
    %Time in each position bin
    Tm=zeros(max(P1(:,2)),max(P1(:,3)));
    for N=1:size(P,1)    
        Tm(P(N,2),P(N,3))=Tm(P(N,2),P(N,3))+(1/30);
    end

    %number of spikes in each position bin is calculated for each cell
    SIP=zeros(max([P1(:,2);S1(:,2)]),max([P1(:,3);S1(:,3)]),max(Cell_Number));
    for i = 1:length(S)        
        SIP(S(i,2),S(i,3),S(i,1))= SIP(S(i,2),S(i,3),S(i,1))+1;    
    end

    FR=zeros(size(SIP,1),size(SIP,2),size(SIP,3));
    for N=1:size(SIP,3)
        FR(:,:,N)=SIP(:,:,N)./Tm;
    end
    FR(isnan(FR))=0;
    FR(isinf(FR))=0;
    
    if idir == 1
        for N=1:size(FR,3)
            OpenInFR(:,:,N)=filter2(Two_D_Filter,FR(:,:,N));
        end
    else
        for N=1:size(FR,3)
            OpenOutFR(:,:,N)=filter2(Two_D_Filter,FR(:,:,N));
        end
    end
end
save(thisdir,'OpenInFR','OpenOutFR','-append')