function make_PosteriorfromThetaCandEvents_manypluspfc(thisdir,cutoff)

load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','prospectiveTheta_armcombs',...
    'armsig_mx','prospectiveTheta_celldat','pos','binpos','prospectiveTheta_pcells')
thetaevents = times_armon_thetaof_headingarm_lap_thetahalf_all;

EstBin=0.02;
load(thisdir,'InFR','OutFR','hp_cells','hpinterneurons','spikedata','other_cells')

% load('Ripple_Events')
% Cand = [Ripple_Events(:,1)-.02 Ripple_Events(:,2)+.02];

celldat = squeeze(nanmean(prospectiveTheta_celldat,1));
celldatdiff = (abs(diff(celldat,[],1))./sum(celldat))';
celldatdiff(:,2) = prospectiveTheta_pcells;
celldatdiff(:,3) = armsig_mx(:,6);
incl = celldatdiff(:,2)<.2;
other_cells = other_cells(incl);


Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
spike = spikedata(ismember(spikedata(:,2),Cell_Number),:);
pfcspike = spikedata(ismember(spikedata(:,2),other_cells),:);
th = reshape(prospectiveTheta_armcombs,[size(prospectiveTheta_armcombs,1) 3 2]);

thmx = th(incl,:,:);
% celldatdiff = celldatdiff(incl,:);

FR = InFR+OutFR;
FR = FR(Cell_Number,:);
laps = unique(thetaevents(~isnan(thetaevents(:,6)),6));
for ilap = 1:length(laps)
    if sum(thetaevents(thetaevents(:,6)==laps(ilap),9)>.4)<2
        continue
    end
    Spike = [spike(:,2) spike(:,1)];
    Cand1 = thetaevents(thetaevents(:,6)==laps(ilap),1:2);
    Cand = [Cand1(1,1) Cand1(end,2)];


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
    FR(FR==0) = .0001;


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

    Mat2=term2.*term3;
    Index=round(B/(EstBin));
    Index2 = 4*Index(1)+1:4*Index(end)-3;
    timebins = Cand(1):.005:Cand(2)+.0025;
    [~,~,i]= histcounts(timebins(1:end-1),pos(:,1));
    timebins = Cand(1):.005:Cand(2)+.05;
    [~,~,ic]= histcounts(Cand1(:,1),timebins);
    [~,~,ic2]= histcounts(Cand1(:,2),timebins);
    armof = thetaevents(thetaevents(:,6)==laps(ilap),4);
    armof(thetaevents(thetaevents(:,6)==laps(ilap),9)<0) = NaN;
    as = unique(armof(~isnan(armof)));
    if length(as)<2
        continue
    end
    
    %normalize for each time bin
    Mat2 = Mat2./(ones(size(Mat2,1),1)*sum(Mat2,1));
%     Mat2 = (Mat2-nanmin(Mat2))./(nanmax(Mat2)-nanmin(Mat2));
    %make times when there is zero spiking equal to zeros/NaNs all the way down
    Mat2(:,sum(binspike>0)==0) = NaN;
    
    pfcspike2 = pfcspike(pfcspike(:,1)>=timebins(1) & pfcspike(:,1)<=timebins(end),:);
    [~,~,ipfc]= histcounts(pfcspike2(:,1),timebins);
    
    thisarm = unique(thetaevents(thetaevents(:,6)==laps(ilap),3));
    if size(thmx,1)>1
        [~,mx] = max(squeeze(thmx(:,thisarm,:)),[],2);
    else
        [~,mx] = max(squeeze(thmx(:,thisarm,:)),[],1);
    end
    armsn = setdiff(1:3,thisarm);
    plot1 = other_cells(mx==1);
    plot2 = other_cells(mx==2);
    thratio = thetaevents(thetaevents(:,6)==laps(ilap),[9]);
    if sum(armof==as(2) & thratio>cutoff)==0 || sum(armof==as(1) & thratio>cutoff)==0
        continue
    end
    
    figure; hold on
    
    
    aa = subplot(3,1,1:2);
    imagesc(Mat2); 
%     colormap jet; 
    colormap gray; 
    cm = get(gca,'colormap');
    cm2 = cm(end:-1:1,:);
    set(gca,'colormap',cm2)
    axis xy; hold on  
    
%     plot([ic(armof==as(2)) ic(armof==as(2))],[0 163],'k-','LineWidth',1);  
%     plot([ic(armof==as(1)) ic(armof==as(1))]',[0 163],'r-','LineWidth',1);
%     asp1 = plot([ic2(armof==as(1)) ic2(armof==as(1))]',[0 163],'r-','LineWidth',3);
%     asp2 = plot([ic2(armof==as(2)) ic2(armof==as(2))],[0 163],'k-','LineWidth',3);     
    
    plot([ic(armof==as(2) & thratio>cutoff) ic(armof==as(2) & thratio>cutoff)]'+1,[0 163],'k-','LineWidth',3);  
    plot([ic(armof==as(1) & thratio>cutoff) ic(armof==as(1) & thratio>cutoff)]'+1,[0 163],'r-','LineWidth',3);
    asp1 = plot([ic2(armof==as(1) & thratio>cutoff) ic2(armof==as(1) & thratio>cutoff)]'-1,[0 163],'r-','LineWidth',3);
    asp2 = plot([ic2(armof==as(2) & thratio>cutoff) ic2(armof==as(2) & thratio>cutoff)]'-1,[0 163],'k-','LineWidth',3);     
    plot(Index2,binpos(i),'r*')            
    legend([asp1(1) asp2(1)],num2str(as),'Location','northwest')    
    set(gca,'xtick',mean([ic(thratio>cutoff) ic2(thratio>cutoff)],2),'xticklabel',num2str(round(thratio(thratio>cutoff),2)));
    xlabel('Non-Local Proportion (>.34)')
    ylabel('Binned Position')
    set(gca,'FontSize',18)    
    set(gca,'clim',[0 .1])
    title(['Rat on ' num2str(thisarm) ' heading to Arm ' ...
        num2str(unique(thetaevents(thetaevents(:,6)==laps(ilap),5))) ', lap ' num2str(ilap)])
    
    
    bb = subplot(3,1,3);
    hold on
    yl = [.5 length(plot1)+length(plot2)+.5];
    
    icy2 = ic2(armof==as(1))+20; %was 40
    icy = ic2(armof==as(1) );    
    patch([icy icy2 icy2 icy]',(ones(length(icy),1)*[yl(1) yl(1) yl(2) yl(2)])','r','FaceAlpha',.1,'EdgeAlpha',0)
    icy2 = ic2(armof==as(2))+20;
    icy = ic2(armof==as(2));    
    patch([icy icy2 icy2 icy]',(ones(length(icy),1)*[yl(1) yl(1) yl(2) yl(2)])','k','FaceAlpha',.1,'EdgeAlpha',0)    
    icy2 = ic2(armof==as(1) & thratio>cutoff)+20; %was 40
    icy = ic2(armof==as(1)  & thratio>cutoff);    
    patch([icy icy2 icy2 icy]',(ones(length(icy),1)*[yl(1) yl(1) yl(2) yl(2)])','r','FaceAlpha',.2,'EdgeAlpha',0)
    icy2 = ic2(armof==as(2) & thratio>cutoff)+20;
    icy = ic2(armof==as(2) & thratio>cutoff);    
    patch([icy icy2 icy2 icy]',(ones(length(icy),1)*[yl(1) yl(1) yl(2) yl(2)])','k','FaceAlpha',.2,'EdgeAlpha',0)    
    for icell = 1:length(plot1)
        ind = pfcspike2(:,2)==plot1(icell);
        plot(ipfc(ind),icell*ones(sum(ind),1),'ro','MarkerFaceColor','r','MarkerSize',10)
    end
    tc = icell;
    for icell = tc:tc+length(plot2)-1
         ind = pfcspike2(:,2)==plot2(icell-tc+1);
        plot(ipfc(ind),icell*ones(sum(ind),1)+1,'ro','MarkerFaceColor','k','MarkerSize',10)
    end
    ylim(yl)    
    set(gca,'ytick',1:2:length(other_cells),'yticklabel',other_cells(1:2:end)-min(other_cells)+1) %armsn([ones(length(plot1),1);2*ones(length(plot2),1)]))
    ylabel('PFC Cell')
    xlabel('Time (Seconds)')
    set(gca,'FontSize',18)   
    linkaxes([aa bb],'x')
    set(gca,'xtick',1:50:length(timebins),'xticklabel',timebins(1:50:end)-timebins(1))
    set(gca,'xlim',[1 size(Mat2,2)])
    
    set(gcf,'Position',[ 1921         -79        1600        1083])   
    disp(ilap)
%     set(gcf,'renderer','Painters')
%     helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Theta\newExamplePFC_' num2str(thisdir(1:end-4)) '_lap' num2str(ilap)])
    close gcf
end