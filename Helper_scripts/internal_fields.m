function [SSDarm,alldat,alldatS] = internal_fields(thisdir,EstBin,cellcutoff,cmcutoff,toplot,numshuff)
%run inside internal_fields_all
%decoding
if 1
% EstBin=0.04;
load(thisdir,'RP_pSSDarm','RP_SSDarm','InFR','OutFR','hp_cells','hpinterneurons','spikedata','pos','binpos','armpos','armposindex','vel','other_cells')
newcol = [75 0 130;255 130 0;34 139 34]/255;
% load('Ripple_Events')
% Cand = [Ripple_Events(:,1)-.02 Ripple_Events(:,2)+.02];

Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
% Spike = [spikedata(ismember(spikedata(:,2),Cell_Number),2) spikedata(ismember(spikedata(:,2),Cell_Number),1)];
Spike = [spikedata(:,2) spikedata(:,1)];

FR = InFR+OutFR;
FR = FR(Cell_Number,:);

Cand = [min(spikedata(:,1)) max(spikedata(:,1))];
CandSeq = Cand;
N=round(diff(Cand')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[~,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
end
Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
% clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;
times = Cand(1):EstBin/4:Cand(2);
times = times+EstBin/2; times = times(1:end-1);
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for i=1:4
    for CellID=1:length(Cell_Number)   
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
%         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
    end                      
end

clear Start_Time End_Time CellID c i 

FR(sum(FR~=0,2)==0,:) = 1;
FR(FR==0 | FR<0) = .0001;

touse = sum(binspike>0)>=cellcutoff; % has to have at least 5 cells spiking
binspike = binspike(:,touse);
times = times(touse);
TimeBins = length(times);

Number_Of_Bins=size(OutFR,2);
term2=zeros(Number_Of_Bins,TimeBins);

if Number_Of_Bins>TimeBins    
    for TBin=1:TimeBins        
        term2(:,TBin)=prod(FR.^binspike(:,TBin),1);
    end
else
    for PosBin=1:Number_Of_Bins           
        term2(PosBin,:)=prod((FR(:,PosBin)*ones(1,TimeBins)).^binspike,1);
    end
end


term3=exp(-EstBin*sum(FR,1)')*ones(1,TimeBins);

Matrix=term2.*term3;
% Index=round(B/(EstBin));

%normalize for each time bin
Est = Matrix./(ones(size(Matrix,1),1)*sum(Matrix,1));
% %make times when there is zero spiking equal to zeros/NaNs all the way down
% Matrix(:,sum(binspike>0)==0) = NaN;

[~,i] = histc(times,pos(:,1));
inx = find(i==0);
inx2 = inx<floor(length(times)/2);
i(inx(inx2)) = i(find(i~=0,1,'first')); i(inx(~inx2)) = i(find(i~=0,1,'last'));
binposition = binpos(i);
binvel = vel(i);
binarm = armpos(i);

exclude = false(size(binposition));
%excluding the very center of the maze. eh
% for iarm = 1:3
%     a = find(armposindex(:,iarm));
%     exclude(binposition>=a(1) & binposition<=a(3)) = true;
% end
    
[~,decoded_position] = max(Est); % change from max to mean and do seperately for each arm, then take the max

Mat2 = NaN(size(Est,2),3); Mat3 = Mat2;
for iarm = 1:3
    Mat = Est(armposindex(:,iarm),:);    
    I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
    Mat2(:,iarm) = I2+find(armposindex(:,iarm),1,'first')-1;
    Mat3(:,iarm) = sum(Est(armposindex(:,iarm),:))./sum(armposindex(:,iarm));
end
[~,mm] = max(Mat3,[],2);
decoded_position2 = decoded_position;
for iarm = 1:3
    decoded_position2(mm==iarm) = Mat2(mm==iarm,iarm);
end

t = times(~exclude)';
dp = decoded_position2(~exclude)';
bp = binposition(~exclude);
bv = binvel(~exclude);
ba = binarm(~exclude);
da = mm(~exclude);
diffarm = ba~=da;
toadd = [0 find(armposindex(:,1),1,'last') find(armposindex(:,2),1,'last')]';
err = abs(dp-bp);
err(diffarm) = (dp(diffarm)-toadd(da(diffarm)))+(bp(diffarm)-toadd(ba(diffarm))); 
%when it decodes a different arm, add the distance between them in real space instead of
%in linearalized space

% figure; histogram(err)

%low velocity, high velocity,
%local, nonlocal
%PFC cell 100ms afterwards, FR

other_cells = other_cells(RP_pSSDarm<.05);
SSDarm = RP_SSDarm(RP_pSSDarm<.05);
% SSDarm = RP_SSDarm;

alldat = NaN(length(SSDarm),3,2,2,3,2); 
alldatS = NaN(length(SSDarm),3,2,2,3,2,numshuff); 
%cells,arms,velocity (low, high), locality (local, nonlocal), measurment (mean, median, cmcutoff), raw or smoothed
Spike = spikedata(ismember(spikedata(:,2),other_cells),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mean
%data
if 1
errcutoff = mean(err); 
fr = NaN(length(other_cells),size(armposindex,1),2,2); 
frS = NaN(length(other_cells),size(armposindex,1),2,2,numshuff); 
occ = NaN(size(armposindex,1),2,2);
for ivel = 1:2
    if ivel==1; velind = bv<5; elseif ivel==2; velind = bv>5; end
    for ilocal = 1:2
        if ilocal ==1; locind = err<errcutoff; elseif ilocal ==2; locind = err>errcutoff; end
        t2 = t(velind & locind);        
        [~,dp2] = histc(dp(velind & locind),1:size(armposindex,1));
        spks = NaN(length(other_cells),length(t2)); 
        for itt = 1:length(t2)
            s = Spike(Spike(:,1)>=t2(itt) & Spike(:,1)<(t2(itt)+EstBin),2);
            spks(:,itt) = histc(s,other_cells);
        end
        
        for idp = 1:size(armposindex,1)
            fr(:,idp,ivel,ilocal) = sum(spks(:,dp2==idp),2)./(sum(dp2==idp)*EstBin); %EstBin);
        end
        
        dp2s = NaN(length(dp2),numshuff);
        for ishuff = 1:numshuff
            dp2s(:,ishuff) = dp2(randperm(length(dp2)));
        end
        spkss = repmat(spks,[1 1 numshuff]);
       for icell = 1:size(spkss,1)
           dat = squeeze(spkss(icell,:,:));
           for idp = 1:size(armposindex,1)
                frS(icell,idp,ivel,ilocal,:) = sum(dat(dp2s==idp),1)./(sum(dp2==idp)*.1); 
           end
       end         
        occ(:,ivel,ilocal) = hist(dp2,1:size(armposindex,1));
    end
end

Filter=fspecial('gaussian',[1 10],1); 
fr2 = NaN(size(fr)); fr3 = fr2;
fr3S = NaN(size(frS)); fr2S = fr3S;
for ivel = 1:2
    for ilocal = 1:2
        for icell = 1:size(fr,1)
            for iarm = 1:3
                fr(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                FR=fillmissing(fr(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                FR(FR<0) = 0;
    %             FR2 = [FR(1)*ones(20,1);FR;FR(end)*ones(20,1)];        
                FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                FRnew3 = filter2(Filter,FR2');        
                FRnew4 = mean([FRnew2' FRnew3'],2);
                fr2(icell,armposindex(:,iarm),ivel,ilocal)=FRnew4(21:end-20);
                fr3(icell,armposindex(:,iarm),ivel,ilocal) = fr2(icell,armposindex(:,iarm),ivel,ilocal);
                fr3(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                for ishuff = 1:numshuff
                   frS(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                    FR=fillmissing(frS(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                    FR(FR<0) = 0;   
                    FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                    FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                    FRnew3 = filter2(Filter,FR2');        
                    FRnew4 = mean([FRnew2' FRnew3'],2);
                    fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff)=FRnew4(21:end-20);
                    fr3S(icell,armposindex(:,iarm),ivel,ilocal,ishuff) = fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff);
                    fr3S(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                end
            end
        end
        for iarm = 1:3
            alldat(:,iarm,ivel,ilocal,1,1) = squeeze(nansum(fr(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldat(:,iarm,ivel,ilocal,1,2) = squeeze(nansum(fr3(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,1,1,:) = squeeze(nansum(frS(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,1,2,:) = squeeze(nansum(fr3S(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
        end
        
    end
end
end
%add to alldat

fr2 = fr3;
% Plot cells
if toplot
for icell = 1:length(other_cells)
    figure; hold on
    subplot(2,1,1); hold on
    yyaxis left
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')
    
    yyaxis right
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')
    xlabel('cm along track (Low Velocity)')        
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    title([thisdir(1) thisdir(3:end-4) ' PFC Cell ' num2str(other_cells(icell))])
    
    subplot(2,1,2); hold on
    yyaxis left
    for iarm = 1:3
        a = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')   
    yyaxis right
    for iarm = 1:3
        b = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')    
    xlabel('cm along track (High Velocity)')
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    set(gcf,'Position',[2050          33        1158         865])
    legend([a,b],'Local','Non-Local')
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\InternalFields\')
        mkdir('E:\XY_matdata\Figures\ForPaper\InternalFields\')
    end
    set(gcf,'renderer','Painters')
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Mean_' thisdir(1:end-4) '_Cell' num2str(other_cells(icell))])
end
end

%plot occupancy
if toplot
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')   
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Mean_' thisdir(1:end-4) '_Occupancy'])
end

%plot normalized occupancy
if toplot
    
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')
ylim([0 1])

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Normalized Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k') 
ylim([0 1])
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Mean_' thisdir(1:end-4) '_OccupancyNormalized'])
end

%plot error
if toplot
figure; histogram(err,'FaceColor','k'); hold on
yl = get(gca,'ylim');
xlabel('Error')
ylabel('Counts')
plot([mean(err) mean(err)],yl,'r-','LineWidth',3)
title([thisdir(1) thisdir(3:end-4) ' Mean Error (Cutoff) is ' num2str(mean(err))])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Mean_' thisdir(1:end-4) '_Cutoff'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Median
%data
if 1
errcutoff = median(err);
fr = NaN(length(other_cells),size(armposindex,1),2,2); 
frS = NaN(length(other_cells),size(armposindex,1),2,2,numshuff); 
occ = NaN(size(armposindex,1),2,2);
for ivel = 1:2
    if ivel==1; velind = bv<5; elseif ivel==2; velind = bv>5; end
    for ilocal = 1:2
        if ilocal ==1; locind = err<errcutoff; elseif ilocal ==2; locind = err>errcutoff; end
        t2 = t(velind & locind);        
        [~,dp2] = histc(dp(velind & locind),1:size(armposindex,1));
        spks = NaN(length(other_cells),length(t2)); 
        for itt = 1:length(t2)
            s = Spike(Spike(:,1)>=t2(itt) & Spike(:,1)<(t2(itt)+EstBin),2);
            spks(:,itt) = histc(s,other_cells);
        end
        
        for idp = 1:size(armposindex,1)
            fr(:,idp,ivel,ilocal) = sum(spks(:,dp2==idp),2)./(sum(dp2==idp)*EstBin); %EstBin);
        end
        
        dp2s = NaN(length(dp2),numshuff);
        for ishuff = 1:numshuff
            dp2s(:,ishuff) = dp2(randperm(length(dp2)));
        end
        spkss = repmat(spks,[1 1 numshuff]);
       for icell = 1:size(spkss,1)
           dat = squeeze(spkss(icell,:,:));
           for idp = 1:size(armposindex,1)
                frS(icell,idp,ivel,ilocal,:) = sum(dat(dp2s==idp),1)./(sum(dp2==idp)*.1); 
           end
       end         
        occ(:,ivel,ilocal) = hist(dp2,1:size(armposindex,1));
    end
end

Filter=fspecial('gaussian',[1 10],1); 
fr2 = NaN(size(fr)); fr3 = fr2;
fr3S = NaN(size(frS)); fr2S = fr3S;
for ivel = 1:2
    for ilocal = 1:2
        for icell = 1:size(fr,1)
            for iarm = 1:3
                fr(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                FR=fillmissing(fr(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                FR(FR<0) = 0;
    %             FR2 = [FR(1)*ones(20,1);FR;FR(end)*ones(20,1)];        
                FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                FRnew3 = filter2(Filter,FR2');        
                FRnew4 = mean([FRnew2' FRnew3'],2);
                fr2(icell,armposindex(:,iarm),ivel,ilocal)=FRnew4(21:end-20);
                fr3(icell,armposindex(:,iarm),ivel,ilocal) = fr2(icell,armposindex(:,iarm),ivel,ilocal);
                fr3(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                for ishuff = 1:numshuff
                   frS(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                    FR=fillmissing(frS(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                    FR(FR<0) = 0;   
                    FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                    FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                    FRnew3 = filter2(Filter,FR2');        
                    FRnew4 = mean([FRnew2' FRnew3'],2);
                    fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff)=FRnew4(21:end-20);
                    fr3S(icell,armposindex(:,iarm),ivel,ilocal,ishuff) = fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff);
                    fr3S(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                end
            end
        end
         for iarm = 1:3
            alldat(:,iarm,ivel,ilocal,2,1) = squeeze(nansum(fr(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldat(:,iarm,ivel,ilocal,2,2) = squeeze(nansum(fr3(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,2,1,:) = squeeze(nansum(frS(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,2,2,:) = squeeze(nansum(fr3S(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
        end
    end
end
end

%add to alldat

fr2 = fr3;
% plot cells
if toplot
for icell = 1:length(other_cells)
    figure; hold on
    subplot(2,1,1); hold on
    yyaxis left
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')
    
    yyaxis right
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')
    xlabel('cm along track (Low Velocity)')        
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    title([thisdir(1) thisdir(3:end-4) ' PFC Cell ' num2str(other_cells(icell))])
    
    subplot(2,1,2); hold on
    yyaxis left
    for iarm = 1:3
        a = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')   
    yyaxis right
    for iarm = 1:3
        b = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')    
    xlabel('cm along track (High Velocity)')
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    set(gcf,'Position',[2050          33        1158         865])
    legend([a,b],'Local','Non-Local')
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\InternalFields\')
        mkdir('E:\XY_matdata\Figures\ForPaper\InternalFields\')
    end
    set(gcf,'renderer','Painters')
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Median_' thisdir(1:end-4) '_Cell' num2str(other_cells(icell))])
end
end

%plot occupancy
if toplot
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')   
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Median_' thisdir(1:end-4) '_Occupancy'])
end

%plot normalized occupancy
if toplot
    
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')
ylim([0 1])

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Normalized Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k') 
ylim([0 1])
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Median_' thisdir(1:end-4) '_OccupanyNormalized'])
end

%plot error
if toplot
figure; histogram(err,'FaceColor','k'); hold on
yl = get(gca,'ylim');
xlabel('Error')
ylabel('Counts')
plot([median(err) median(err)],yl,'r-','LineWidth',3)
title([thisdir(1) thisdir(3:end-4) ' Median Error (Cutoff) is ' num2str(median(err))])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_Median_' thisdir(1:end-4) '_Cutoff'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cm cutoff

%data
if 1
errcutoff = cmcutoff;

fr = NaN(length(other_cells),size(armposindex,1),2,2); 
frS = NaN(length(other_cells),size(armposindex,1),2,2,numshuff); 
occ = NaN(size(armposindex,1),2,2);
for ivel = 1:2
    if ivel==1; velind = bv<5; elseif ivel==2; velind = bv>5; end
    for ilocal = 1:2
        if ilocal ==1; locind = err<errcutoff; elseif ilocal ==2; locind = err>errcutoff; end
        t2 = t(velind & locind);        
        [~,dp2] = histc(dp(velind & locind),1:size(armposindex,1));
        spks = NaN(length(other_cells),length(t2)); 
        for itt = 1:length(t2)
            s = Spike(Spike(:,1)>=t2(itt) & Spike(:,1)<(t2(itt)+EstBin),2);
            spks(:,itt) = histc(s,other_cells);
        end
        
        for idp = 1:size(armposindex,1)
            fr(:,idp,ivel,ilocal) = sum(spks(:,dp2==idp),2)./(sum(dp2==idp)*EstBin); %EstBin);
        end
        
        dp2s = NaN(length(dp2),numshuff);
        for ishuff = 1:numshuff
            dp2s(:,ishuff) = dp2(randperm(length(dp2)));
        end
        spkss = repmat(spks,[1 1 numshuff]);
       for icell = 1:size(spkss,1)
           dat = squeeze(spkss(icell,:,:));
           for idp = 1:size(armposindex,1)
                frS(icell,idp,ivel,ilocal,:) = sum(dat(dp2s==idp),1)./(sum(dp2==idp)*.1); 
           end
       end         
        occ(:,ivel,ilocal) = hist(dp2,1:size(armposindex,1));
    end
end

Filter=fspecial('gaussian',[1 10],1); 
fr2 = NaN(size(fr)); fr3 = fr2;
fr3S = NaN(size(frS)); fr2S = fr3S;
for ivel = 1:2
    for ilocal = 1:2
        for icell = 1:size(fr,1)
            for iarm = 1:3
                fr(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                FR=fillmissing(fr(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                FR(FR<0) = 0;
    %             FR2 = [FR(1)*ones(20,1);FR;FR(end)*ones(20,1)];        
                FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                FRnew3 = filter2(Filter,FR2');        
                FRnew4 = mean([FRnew2' FRnew3'],2);
                fr2(icell,armposindex(:,iarm),ivel,ilocal)=FRnew4(21:end-20);
                fr3(icell,armposindex(:,iarm),ivel,ilocal) = fr2(icell,armposindex(:,iarm),ivel,ilocal);
                fr3(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal) = NaN;
                for ishuff = 1:numshuff
                   frS(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                    FR=fillmissing(frS(icell,armposindex(:,iarm),ivel,ilocal),'linear',2)';
                    FR(FR<0) = 0;   
                    FR2 = [FR(20:-1:1);FR;FR(end:-1:end-19)];        
                    FRnew2 = filter2(Filter,FR2(end:-1:1)'); FRnew2 = FRnew2(end:-1:1);
                    FRnew3 = filter2(Filter,FR2');        
                    FRnew4 = mean([FRnew2' FRnew3'],2);
                    fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff)=FRnew4(21:end-20);
                    fr3S(icell,armposindex(:,iarm),ivel,ilocal,ishuff) = fr2S(icell,armposindex(:,iarm),ivel,ilocal,ishuff);
                    fr3S(icell,(occ(:,ivel,ilocal)<10),ivel,ilocal,ishuff) = NaN;
                end
            end
        end
         for iarm = 1:3
            alldat(:,iarm,ivel,ilocal,3,1) = squeeze(nansum(fr(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldat(:,iarm,ivel,ilocal,3,2) = squeeze(nansum(fr3(:,armposindex(:,iarm),ivel,ilocal),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,3,1,:) = squeeze(nansum(frS(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
            alldatS(:,iarm,ivel,ilocal,3,2,:) = squeeze(nansum(fr3S(:,armposindex(:,iarm),ivel,ilocal,:),2))./sum(armposindex(:,iarm) & occ(:,ivel,ilocal)>=10);
        end
    end
end
end
%add to alldat

fr2 = fr3;
%plot cells
if toplot
for icell = 1:length(other_cells)
    figure; hold on
    subplot(2,1,1); hold on
    yyaxis left
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')
    
    yyaxis right
    for iarm = 1:3
        plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')
    xlabel('cm along track (Low Velocity)')        
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    title([thisdir(1) thisdir(3:end-4) ' PFC Cell ' num2str(other_cells(icell))])
    
    subplot(2,1,2); hold on
    yyaxis left
    for iarm = 1:3
        a = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Local (Dashed)','Color','k')
    set(gca,'YColor','k')   
    yyaxis right
    for iarm = 1:3
        b = plot(find(armposindex(:,iarm)),fr2(icell,armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3);
    end
    ylabel('Non-Local (Solid)','Color','k')
    set(gca,'YColor','k')    
    xlabel('cm along track (High Velocity)')
    xlim([1 size(armposindex,1)])
    set(gca,'FontSize',18)
    set(gcf,'Position',[2050          33        1158         865])
    legend([a,b],'Local','Non-Local')
    if ~isfolder('E:\XY_matdata\Figures\ForPaper\InternalFields\')
        mkdir('E:\XY_matdata\Figures\ForPaper\InternalFields\')
    end
    set(gcf,'renderer','Painters')
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_' num2str(cmcutoff) 'cmcutoff_' thisdir(1:end-4) '_Cell' num2str(other_cells(icell))])
end
end

%plot occupancy
if toplot
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
plot(find(occ(:,1,1)<10),occ(occ(:,1,1)<10,1,1),'r*')
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
plot(find(occ(:,1,2)<10),occ(occ(:,1,2)<10,1,2),'r*')
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1),'--','color',newcol(iarm,:),'LineWidth',3)
end
plot(find(occ(:,2,1)<10),occ(occ(:,2,1)<10,2,1),'r*')
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')   
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2),'-','color',newcol(iarm,:),'LineWidth',3)
end
plot(find(occ(:,2,2)<10),occ(occ(:,2,2)<10,2,2),'r*')
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_' num2str(cmcutoff) 'cmcutoff_' thisdir(1:end-4) '_Occupancy'])
end

%plot normalized occupancy
if toplot
    
figure; hold on
subplot(2,1,1); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k')
ylim([0 1])

yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),1,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')
xlabel('cm along track (Low Velocity)')        
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
title([thisdir(1) thisdir(3:end-4) ' Normalized Occupancy'])

subplot(2,1,2); hold on
yyaxis left
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,1)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'--','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Local (Dashed)','Color','k')
set(gca,'YColor','k') 
ylim([0 1])
yyaxis right
for iarm = 1:3
    plot(find(armposindex(:,iarm)),occ(armposindex(:,iarm),2,2)./sum(sum(occ(armposindex(:,iarm),:,:),2),3),'-','color',newcol(iarm,:),'LineWidth',3)
end
ylabel('Non-Local (Solid)','Color','k')
set(gca,'YColor','k')    
xlabel('cm along track (High Velocity)')
xlim([1 size(armposindex,1)])
ylim([0 1])
set(gca,'FontSize',18)
set(gcf,'Position',[2050          33        1158         865])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_' num2str(cmcutoff) 'cmcutoff_' thisdir(1:end-4) '_OccupancyNormalized'])
end

%plot error
if toplot
figure; histogram(err,'FaceColor','k'); hold on
yl = get(gca,'ylim');
xlabel('Error')
ylabel('Counts')
plot([cmcutoff cmcutoff],yl,'r-','LineWidth',3)
title([thisdir(1) thisdir(3:end-4) ' Cutoff is ' num2str(cmcutoff)])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\EstBin' num2str(EstBin) '_CellCutoff' num2str(cellcutoff) '_' num2str(cmcutoff) 'cmcutoff_' thisdir(1:end-4) '_Cutoff'])

end