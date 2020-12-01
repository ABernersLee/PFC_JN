function Figure3_Internal2d(dirs,igroup,savefolder,toplotcells,tosave,EstBin,velcutoff,cellcutoff,spikecutoff,numshuff,errcutoff,binsize,PFCdelay)
%% get data out
cd(dirs.homedir)
d2 = dir('*.mat');
% binsize = 20; 
%20 = ~8x8cm % if figures dont have binsize on the end of the title, they are 20

% errcutoff use 1(8) 1.5 (12) and 2.5 (20)
% how big is this? I wrote that I used 5cm x 5cm
% bin of ~12.5 would be 5x5

% %testing
% for id = 1:11; load(d2(id).name,'params');armlen(1,id) = sqrt((params.arms(1,3)-params.arms(1,1)).^2+(params.arms(1,4)-params.arms(1,2)).^2)./161;
% armlen(2,id) = sqrt((params.arms(2,3)-params.arms(2,1)).^2+(params.arms(2,4)-params.arms(2,2)).^2)./81;
% armlen(3,id) = sqrt((params.arms(3,3)-params.arms(3,1)).^2+(params.arms(3,4)-params.arms(3,2)).^2)./81; end

%mean(mean(armlen)) = 2.5293 pixels/cm
    % 20 pixels = 7.9073 cm

% How did I think it was 4 pixels/cm before?

%testing another way:
% for id = 1:11; load(d2(id).name,'cm_conv'); cm(id,1) = 1./cm_conv; end
% mean(cm) =  2.5233 pixels/cm
    % 20 pixels = 7.9262 cm

%     PFCdelay = 0.06; %EstBin; %was EstBin (.04), checking with .06 on 9/18/2020 reviewer comment

pix2cm = 2.5;
% cm = num*binsize/pix2cm
SigmaField = 2;
m = errcutoff;
Two_D_Filter=fspecial('gaussian',[3 3],SigmaField);
vellab = {'Low Velocity';'High Velocity'};
loclab = {'Local';'Non-Local'};
if ~isfolder([savefolder '\Figure3\related\'])
    mkdir([savefolder '\Figure3\related\'])
end
figlab = ['_BinSize' num2str(EstBin) '_VelCutoff' num2str(velcutoff) '_cellcutoff' num2str(cellcutoff) '_spikecutoff' num2str(spikecutoff) '_errcutoff' num2str(m) '_numshuff' num2str(numshuff) '_binsize' num2str(binsize) '_PFCdelay' num2str(PFCdelay)];

repall = []; sigcells = []; othall = []; othall_lowvel = []; othallS = [];
toterr  = []; numbinsused = NaN(size(d2,1),2);
for id = 1:size(d2,1)
    load(d2(id).name,'other_cells_touse')
    if sum(other_cells_touse(:,igroup))==0
        continue
    end
    
    clear OpenFRinternal
    tocalc = true;
    if tocalc% get the internal fields        

        load(d2(id).name,'OpenFR','hp_cells','hpinterneurons','spikedata','pos','armpos','armposindex','vel','other_cells','other_cells_touse','RP_CandEventTimes')

        OGOpenFR = OpenFR;
        pos2 = pos(:,2:3);
        pos2(:,1) = ceil(pos2(:,1)/binsize);
        pos2(:,2) = ceil(pos2(:,2)/binsize);
        pos2 = pos2-min(pos2)+1;
        tm = diff(pos(:,1)); tm = [tm;tm(end)];

        pos3 = sub2ind(max(pos2),pos2(:,1),pos2(:,2));
        pos3(vel<velcutoff) = NaN;
        ind = 1:max(pos2(:,1))*max(pos2(:,2));
        occ = NaN(size(ind,2),1);
        f = zeros(size(ind,2),max(spikedata(:,2)));
        for j = 1:length(ind)
            occ(ind(j)) = sum(tm(pos3==ind(j)));
            [h,k] = histc(spikedata(pos3(spikedata(:,3)) == ind(j),2),1:max(spikedata(:,2)));
            if ~isempty(k)
                f(ind(j),k) = h(k);            
            end
            if rem(j,1000)==0
                disp(['j = ' num2str(j)])
            end
        end
        occ2 = reshape(occ,[max(pos2)]);        
        f2 = reshape(f,[max(pos2) size(f,2)]);
        OpenFR = f2./repmat(occ2,[1 1 size(f2,3)]);
        
        
        
        

        Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));

        for icell = 1:size(OpenFR,3)
            dat = OpenFR(:,:,icell);                
            datnan = isnan(dat);
            dat(datnan) = 0;
            dat = filter2(Two_D_Filter,dat);
%                 smx = max(max(dat));
            dat(datnan) = NaN;
            OpenFR(:,:,icell) = dat;
        end

        FR1 = reshape(OpenFR(:,:,Cell_Number),[size(OpenFR,1)*size(OpenFR,2) length(Cell_Number)])';


        nanFR = sum(~isnan(FR1))~=0;
        nanIND = find(nanFR);
        FR = FR1(:,nanFR);

        Cand = [min(spikedata(:,1)) max(spikedata(:,1))];
        CandSeq = Cand;
        N=round(diff(Cand')/EstBin)';
        t=find(mod(N,2));
        Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
        t=find(~mod(N,2));
        Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
        clear N t
        B=cumsum([0 diff(Cand')]);


        Spike = [spikedata(:,2) spikedata(:,1)];
        [~,I]=histc(Spike(:,2),sortrows(Cand(:)));    
        for i=1:2:max(I)
            Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
        end
        Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
        clear I i

        Start_Time=0;
        End_Time=B(end);
        TimeBins=round((End_Time-Start_Time)/EstBin)*4;
        times = Cand(1):EstBin/4:Cand(2);
        times = times+EstBin/2; times = times(1:end-1);

%             Bayesian Decoding - 5ms moving step, 20ms window estimate
        binspike=zeros(length(Cell_Number),TimeBins);
        for i=1:4
            for CellID=1:length(Cell_Number)   
                c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
                binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
        %         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
            end                      
        end            
        FR(sum(FR~=0,2)==0,:) = 1;
        FR(FR==0 | FR<0) = .0001;

        touse = sum(binspike>0)>=cellcutoff; % has to have at least 5 cells spiking
        touse = touse & sum(binspike)>=spikecutoff; %has to have at least x spikes 
        binspike = binspike(:,touse);
        times = times(touse);
        TimeBins = length(times);
        numbinsused(id,:) = [sum(touse) sum(~touse)];

        Number_Of_Bins=size(FR,2);
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
        clear term2 term3

        %normalize for each time bin
        Est = Matrix./(ones(size(Matrix,1),1)*sum(Matrix,1));      
        clear Matrix

        exclude = false(length(times),1);

        [~,i] = histc(times,pos(:,1));
        clear times
        inx = find(i==0);
        exclude(inx) = true;
        i(inx) = [];


        pos2 = pos(:,2:3);
        pos2(:,1) = ceil(pos2(:,1)/binsize);
        pos2(:,2) = ceil(pos2(:,2)/binsize);
        pos2 = pos2-min(pos2)+1;

        binposition = pos2(i,:);
        binvel = vel(i);
        binarm = armpos(i);

        [~,decoded_position1] = max(Est); 
        clear Est
        [X,Y] = ind2sub(max(pos2),nanIND(decoded_position1)');
        clear decoded_position1
        decoded_position2 = [X Y];
        clear X Y

        dp = decoded_position2(~exclude,:);
        clear decoded_position2
        bp = binposition(~exclude,:);
        clear binposition
        bv = binvel(~exclude);
        clear binvel

        %excluding the very center of the maze. eh
        % for iarm = 1:3
        %     a = find(armposindex(:,iarm));
        %     exclude(binposition>=a(1) & binposition<=a(3)) = true;
        % end

%             figure; hold on; plot(bp(ba==1,1),bp(ba==1,2),'r.');plot(bp(ba==2,1),bp(ba==2,2),'b.'); plot(bp(ba==3,1),bp(ba==3,2),'k.')          

        %arm rat truely is on
        baX = binarm(~exclude);
        ba1 = armpos;
        pind = max([dp; pos2]);
        x = sub2ind(pind,pos2(:,1),pos2(:,2)); %real linear indicies
        x2 = sub2ind(pind,bp(:,1),bp(:,2)); %real linear indicies excluding bins
        y = sub2ind(pind,dp(:,1),dp(:,2)); %decoded linear indicies   
        inds = cell(3,1);

        inds{1} = unique(x(ba1==1));
        inds{2} = unique(x(ba1==2));
        inds{3} = unique(x(ba1==3));

%             %arm rat is decoded to be on
        da = NaN(size(y));  
        ba = NaN(size(baX));
        for iarm = 1:3
            da(ismember(y,inds{iarm})) = iarm;
            ba(ismember(x2,inds{iarm})) = iarm;
        end
%             figure; hold on; plot(bp(:,1),bp(:,2),'k.');for iarm = 1:3; plot(bp(ismember(x2,inds{iarm}),1),bp(ismember(x2,inds{iarm}),2),'*'); end
%             figure; hold on; plot(dp(:,1),dp(:,2),'k.');for iarm = 1:3; plot(dp(ismember(y,inds{iarm}),1),dp(ismember(y,inds{iarm}),2),'*'); end            

        err = vecnorm(dp-bp,2,2);

        disp(['err = ' num2str(nanmedian(err)*binsize/pix2cm)])
        toterr = cat(1,toterr,[err(bv>velcutoff,1) ones(sum(bv>velcutoff),1)*id]);
    end

     TimeBins=round((End_Time-Start_Time)/EstBin)*4;
     Cell_Number = 1:max(Spike(:,1));
     Start_Time = Start_Time+PFCdelay; %EstBin;    %added delay *has been EstBin with EstBin as .04, testing with .06 %9/18/2020
     binspike2=zeros(length(Cell_Number),TimeBins);
        for i=1:4
            for CellID=1:length(Cell_Number)   
                c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);                    
                binspike2(CellID,i:4:TimeBins)=c(1:TimeBins/4);
        %         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
            end                      
        end
    clear binspike Spike
    bsp3 = binspike2(:,touse);
    clear binspike2 
    bsp3 = bsp3(:,~exclude);

    ind = 1:pind(1)*pind(2);
    pos4 = sub2ind(pind,dp(:,1),dp(:,2));
%         figure; histogram(err(ismember(pos4,ind)))
%         title(num2str(median(err(ismember(pos4,ind)))))
    OpenFRinternal = NaN(pind(1),pind(2),size(OpenFR,3),2,2);
    OpenFRinternalS = NaN(pind(1),pind(2),size(OpenFR,3),2,2,numshuff);
    for ivel = 1:2
        tic
        if ivel==1; velind = bv<velcutoff; elseif ivel==2; velind = bv>velcutoff; end
        
        pos3 = pos4;
        pos3(~velind,:) = NaN; %restrict to times with the correct velocity and locality

        occ = NaN(size(ind,2),2);
        occS = NaN(size(ind,2),2,numshuff);
        f = zeros(size(ind,2),size(OpenFR,3),2);
        fS = zeros(size(ind,2),size(OpenFR,3),2,numshuff);
        mind1 = err<m & ~isnan(pos3);
        mind2 = err>m & ~isnan(pos3);

        posNL = pos3(mind2);
        posShuff = NaN(size(pos3,1),numshuff);


        out = cell2mat(arrayfun(@(dummy) randperm(length(posNL)), 1:numshuff, 'UniformOutput', false)');
        out = out';
        out2 = repmat(posNL,[1 numshuff]);
        posShuffT = out2(out);
        while any(sum(posShuffT~=posNL)==0)
            out = cell2mat(arrayfun(@(dummy) randperm(length(posNL)), 1:numshuff, 'UniformOutput', false)');
            out = out';
            out2 = repmat(posNL,[1 numshuff]);
            posShuffT = out2(out);
        end            
        posShuff(repmat(mind2,[1 numshuff])) = posShuffT;
        clear out out2 posShuffT

%         disp(num2str([(size(bsp3,2)*10) (size(bsp3,2)*10)<3e6]))
        clear bsp4
        if (size(bsp3,1)*size(bsp3,2))>6e7 && binsize<20
            quickway = false;
            halfind = round(size(bsp3,1)/2);
        else
            bsp4 = repmat(bsp3,[1 1 10]);
            quickway = true;                    
        end
        [user,~] = memory;
        
        disp([num2str(round(user.MemUsedMATLAB,3,'significant'))])
        if user.MemUsedMATLAB>21e9 %was 12e9
           disp('out of memory')
        end

        for j = 1:length(ind)       
            if sum(pos3==ind(j))==0
                continue
            end
            for ilocal = 1:2      
                if ilocal == 1; mind = mind1; elseif ilocal ==2; mind = mind2; end
                occ(ind(j),ilocal) = sum(pos3==ind(j) & mind)*EstBin;
                spkb = bsp3(:,pos3==ind(j) & mind);
                f(ind(j),:,ilocal) = sum(spkb,2);
                clear spkb
            end

            %shuffle decoded positions of nonlocal    - broken up
            %because of memory
            ilocal = 2;                                
            occS(ind(j),ilocal,:) = sum(posShuff==ind(j) & mind2)*EstBin;  
            spkbP = NaN(size(bsp3,1),numshuff);
%                 fun = @(a,b) sum(a(:,posShuff(:,b)==ind(j) & mind2),2);
%                 spkbP = bsxfun(fun,bsp3,1:numshuff);
            if quickway
                portions = 0:10:numshuff;                        
                for ii = 1:length(portions)-1
                    spkbS = reshape(bsp4(:,posShuff(:,portions(ii)+1:portions(ii+1))==ind(j) & repmat(mind2,[1 10])),[size(bsp4,1) sum(posShuff(:,1)==ind(j) & mind2) 10]);
                    spkbP(:,portions(ii)+1:portions(ii+1)) = sum(spkbS,2);
                end
            else
%                     for ii = 1:numshuff                            
%                         spkbP(:,ii) = sum(bsp3(:,posShuff(:,ii)==ind(j) & mind2),2);
%                     end
%                  for icell= 1:size(bsp3,1)
%                     bsp5 = repmat(squeeze(bsp3(icell,:))',[1 numshuff]);
%                     spkbP(icell,:) = sum(bsp5,1);
%                  end
                for iii = 1:2
                  if iii==1;  bsp4 = repmat(bsp3(1:halfind,:),[1 1 10]); elseif iii==2;  bsp4 = repmat(bsp3(halfind+1:end,:),[1 1 10]); end
                  portions = 0:10:numshuff;                        
                    for ii = 1:length(portions)-1
                        spkbS = reshape(bsp4(:,posShuff(:,portions(ii)+1:portions(ii+1))==ind(j) & repmat(mind2,[1 10])),[size(bsp4,1) sum(posShuff(:,1)==ind(j) & mind2) 10]);
                        if iii ==1
                            spkbP(1:halfind,portions(ii)+1:portions(ii+1)) = sum(spkbS,2);
                        elseif iii==2
                            spkbP(halfind+1:end,portions(ii)+1:portions(ii+1)) = sum(spkbS,2);
                        end
                    end
                end
                clear bsp4
            end

            fS(ind(j),:,ilocal,:) = spkbP;
            clear spkbP spkbS
        end

        occ2 = reshape(occ,[pind 2]);
        occS2 = reshape(occS,[pind 2 numshuff]);

        if toplotcells
            figure; hold on
            subplot(3,1,1)
            imagesc(occ2(:,:,1)); axis xy; colorbar
            title([vellab{ivel} ' ' loclab{1}])
            subplot(3,1,2)
            imagesc(occ2(:,:,2)); axis xy; colorbar
            title([vellab{ivel} ' ' loclab{2}])
            subplot(3,1,3)
            imagesc(occ2(:,:,2)-occ2(:,:,1)); axis xy; colorbar
            title('Difference in Occupancy')
            set(gcf,'Position',[2562           6         494         879])
            helper_saveandclosefig([savefolder '\Figure3\related\2DInternal_' vellab{ivel} '_' d2(id).name(1:end-4) '_Occ']) 
        end

        f2 = reshape(f,[pind size(f,2) 2]);
        fS2 = reshape(fS,[pind size(f,2) 2 numshuff]);

        OpenFRinternal_temp1 = f2(:,:,:,1)./repmat(occ2(:,:,1),[1 1 size(f2,3)]);
        OpenFRinternal_temp2 = f2(:,:,:,2)./repmat(occ2(:,:,2),[1 1 size(f2,3)]);

        OpenFRinternal_tempS1 = squeeze(fS2(:,:,:,1,:))./repmat(occS2(:,:,1),[1 1 size(fS2,3) size(fS2,5)]);
        OpenFRinternal_tempS2 = squeeze(fS2(:,:,:,2,:))./repmat(occS2(:,:,2),[1 1 size(fS2,3) size(fS2,5)]);
        OpenFRinternal_tempS1 = reshape(OpenFRinternal_tempS1,[size(OpenFRinternal_tempS1,1) size(OpenFRinternal_tempS1,2) size(OpenFRinternal_tempS1,3) 1 1 size(OpenFRinternal_tempS1,4)]);
        OpenFRinternal_tempS2 = reshape(OpenFRinternal_tempS2,[size(OpenFRinternal_tempS2,1) size(OpenFRinternal_tempS2,2) size(OpenFRinternal_tempS2,3) 1 1 size(OpenFRinternal_tempS2,4)]);

        OpenFRinternal(:,:,:,ivel,:) = cat(5,OpenFRinternal_temp1,OpenFRinternal_temp2);

        OpenFRinternalS(:,:,:,ivel,:,:) = cat(5,OpenFRinternal_tempS1,OpenFRinternal_tempS2);                
                
        clear spikedata OpenFRinternal_tempS1 OpenFRinternal_tempS2 fS2 f2 occS2
        ttt = toc;
        disp(['Done with day ' num2str(id) ' ivel ' num2str(ivel) ' in ' num2str(ttt/60) ' min'])
        tic
    end
    
    
    if tosave
%             save(d2(id).name,'OpenFRinternal','-append')        
    else
%         load(d2(id).name,'OpenFRinternal')
    end
    
    
    load(d2(id).name,'other_cells','RP_moduarm','other_cells_touse','pos','armpos','RP_pSSDarm','RP_SSDarm')
    replay = RP_moduarm';

    snan = repmat(sum(sum(~isnan(OpenFRinternal(:,:,:,2,:)),5),3)==0,[1 1 sum(other_cells_touse(:,igroup))]);

    for icell = 1:size(OpenFR,3)
        dat = OpenFR(:,:,icell);                
        datnan = isnan(dat);
        dat(datnan) = 0;
        dat = filter2(Two_D_Filter,dat);
%                 smx = max(max(dat));
        dat(datnan) = NaN;
        OpenFR(:,:,icell) = dat;
    end
    OpenFRPFC = OpenFR(:,:,other_cells(other_cells_touse(:,igroup)));
    OpenFRPFC(snan) = NaN;



    pos2 = pos(:,2:3); clear pos
    pos2(:,1) = ceil(pos2(:,1)/binsize);
    pos2(:,2) = ceil(pos2(:,2)/binsize);
    pos2 = pos2-min(pos2)+1;
    posu = unique(pos2,'rows');
    x = sub2ind(max(pos2),pos2(:,1),pos2(:,2)); %real linear indicies
    x2 = sub2ind(max(pos2),posu(:,1),posu(:,2)); %real linear indicies
    for iarm = 1:3
        inds{iarm} = unique(x(armpos==iarm));
    end
%             %arm rat is on in matrix
    aa = NaN(size(x2));
    for iarm = 1:3
        aa(ismember(x2,inds{iarm})) = iarm;
    end
%             figure; hold on; for iarm = 1:3; plot(posu(aa==iarm,2),posu(aa==iarm,1),'.'); end; legend

    fr = NaN(size(OpenFRPFC,3),3);
    for icell = 1:size(fr,1)
        for iarm = 1:3
            dat = OpenFRPFC(:,:,icell);
            dat = dat(x2);
            fr(icell,iarm) = nanmean(dat(aa==iarm));
        end
    end

    ArmBarPlot_LocalNonLocal_MeanSEM = NaN(size(OpenFRinternal,3),3,2,2);
    ArmBarPlot_LocalNonLocal_MeanSEMS = NaN(size(OpenFRinternal,3),3,2,2,numshuff);
    ArmBarPlot_LocalNonLocal_MeanSEM_lowvel = NaN(size(OpenFRinternal,3),3,2,2);

    for icell = 1:size(OpenFRinternal,3)
        dat = OpenFRinternal(:,:,icell,2,1);
        dat = dat(x2);            
        ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,1,1) = [nanmean(dat(aa==1)) nanmean(dat(aa==2)) nanmean(dat(aa==3))];
        ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,1,2) = [nanstd(dat(aa==1))./sqrt(sum(~isnan(dat(aa==1)))) ...
            nanstd(dat(aa==2))./sqrt(sum(~isnan(dat(aa==2)))) nanstd(dat(aa==3))./sqrt(sum(~isnan(dat(aa==3))))];
        dat = OpenFRinternal(:,:,icell,2,2);
        dat = dat(x2);            
        ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,2,1) = [nanmean(dat(aa==1)) nanmean(dat(aa==2)) nanmean(dat(aa==3))];
        ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,2,2) = [nanstd(dat(aa==1))./sqrt(sum(~isnan(dat(aa==1)))) ...
            nanstd(dat(aa==2))./sqrt(sum(~isnan(dat(aa==2)))) nanstd(dat(aa==3))./sqrt(sum(~isnan(dat(aa==3))))];


        dat = squeeze(OpenFRinternalS(:,:,icell,2,1,:));
        dat = reshape(dat,[size(dat,1)*size(dat,2) numshuff]);
        dat = dat(x2,:);            
        ArmBarPlot_LocalNonLocal_MeanSEMS(icell,:,1,1,:) = [nanmean(dat(aa==1,:))' nanmean(dat(aa==2,:))' nanmean(dat(aa==3,:))']';
        ArmBarPlot_LocalNonLocal_MeanSEMS(icell,:,1,2,:) = [(nanstd(dat(aa==1,:))./sqrt(sum(~isnan(dat(aa==1,:)))))' ...
            (nanstd(dat(aa==2,:))./sqrt(sum(~isnan(dat(aa==2,:)))))' (nanstd(dat(aa==3,:))./sqrt(sum(~isnan(dat(aa==3,:)))))']';
        dat = OpenFRinternalS(:,:,icell,2,2,:);
        dat = reshape(dat,[size(dat,1)*size(dat,2) numshuff]);
        dat = dat(x2,:);            
        ArmBarPlot_LocalNonLocal_MeanSEMS(icell,:,2,1,:) = [nanmean(dat(aa==1,:))' nanmean(dat(aa==2,:))' nanmean(dat(aa==3,:))']';
        ArmBarPlot_LocalNonLocal_MeanSEMS(icell,:,2,2,:) = [(nanstd(dat(aa==1,:))./sqrt(sum(~isnan(dat(aa==1,:)))))' ...
            (nanstd(dat(aa==2,:))./sqrt(sum(~isnan(dat(aa==2,:)))))' (nanstd(dat(aa==3,:))./sqrt(sum(~isnan(dat(aa==3,:)))))']';


        dat = OpenFRinternal(:,:,icell,1,1);
        dat = dat(x2);            
        ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(icell,:,1,1) = [nanmean(dat(aa==1)) nanmean(dat(aa==2)) nanmean(dat(aa==3))];
        ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(icell,:,1,2) = [nanstd(dat(aa==1))./sqrt(sum(~isnan(dat(aa==1)))) ...
            nanstd(dat(aa==2))./sqrt(sum(~isnan(dat(aa==2)))) nanstd(dat(aa==3))./sqrt(sum(~isnan(dat(aa==3))))];
        dat = OpenFRinternal(:,:,icell,1,2);
        dat = dat(x2);            
        ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(icell,:,2,1) = [nanmean(dat(aa==1)) nanmean(dat(aa==2)) nanmean(dat(aa==3))];
        ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(icell,:,2,2) = [nanstd(dat(aa==1))./sqrt(sum(~isnan(dat(aa==1)))) ...
            nanstd(dat(aa==2))./sqrt(sum(~isnan(dat(aa==2)))) nanstd(dat(aa==3))./sqrt(sum(~isnan(dat(aa==3))))];
    end
    
    
    if  tosave
%         save(d2(id).name,'ArmBarPlot_LocalNonLocal_MeanSEM','-append')
    else
%         load(d2(id).name,'ArmBarPlot_LocalNonLocal_MeanSEM')
    end
    
    
    repall = cat(1,repall,replay(other_cells_touse(:,igroup),:));
    
    local = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,1,1);
    nonlocal = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,2,1);    
    othall = cat(1,othall,cat(3,fr(other_cells_touse(:,igroup),:),local(other_cells_touse(:,igroup),:),nonlocal(other_cells_touse(:,igroup),:)));
    
    localS = ArmBarPlot_LocalNonLocal_MeanSEMS(other_cells,:,1,1,:);
    nonlocalS = ArmBarPlot_LocalNonLocal_MeanSEMS(other_cells,:,2,1,:);    
    othallS = cat(1,othallS,cat(3,NaN(sum(other_cells_touse(:,igroup)),size(fr,2),1,size(localS,5)),localS(other_cells_touse(:,igroup),:,:,:),nonlocalS(other_cells_touse(:,igroup),:,:,:)));
    
    local_L = ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(other_cells,:,1,1);
    nonlocal_L = ArmBarPlot_LocalNonLocal_MeanSEM_lowvel(other_cells,:,2,1);
    othall_lowvel = cat(1,othall_lowvel,cat(3,fr(other_cells_touse(:,igroup),:),local_L(other_cells_touse(:,igroup),:),nonlocal_L(other_cells_touse(:,igroup),:)));
    
    sigcells = cat(1,sigcells,[RP_pSSDarm(other_cells_touse(:,igroup)) RP_SSDarm(other_cells_touse(:,igroup)) ]);
    
    
    if toplotcells && igroup ==1
        load(d2(id).name,'spikedata')
        figure; hold on
        dat = OpenFRinternal(:,:,1,1,1);
        imagesc(dat); 
    %         colormap autumn
        colormap gray
        cm = get(gca,'colormap');
    %         cm = cm(end:-1:1,:);
        close gcf
    %     cm2 = [cm; [1 1 1]];
        cm2 = [cm(end-1:-1:1,:); [1 1 1]];
        
        
        figure; hold on        
        datnan = isnan(dat);
        imagesc(datnan); 
        set(gca,'colormap',cm)
        axis xy
        helper_saveandclosefig([savefolder '\Figure3\related\Internal_' d2(id).name(1:end-4)])

        other_cells = other_cells(other_cells_touse(:,igroup));
        RP_pSSDarm = RP_pSSDarm(other_cells_touse(:,igroup));
        other_cells_sig = other_cells(RP_pSSDarm<.05);

        for icell = 1:max(spikedata(:,2))
            if ~ismember(icell,other_cells_sig)
                continue
            end
            if ismember(icell,other_cells_sig); cellid = 'PFCsig'; elseif ismember(icell,other_cells); cellid = 'PFC'; elseif ismember(icell,hpinterneurons); cellid = 'HP INT'; elseif ismember(icell,hp_cells); cellid = 'HP'; end
            
            figure; hold on
            subplot(2,2,1); hold on
            bar(replay(icell==other_cells,:),'FaceColor','w')
            subplot(2,2,2); hold on
            bar(nonlocal_L(icell==other_cells,:),'FaceColor','w')
            subplot(2,2,3); hold on
            bar(local(icell==other_cells,:),'FaceColor','w')
            subplot(2,2,4); hold on
            bar(nonlocal(icell==other_cells,:),'FaceColor','w')
            set(gcf,'renderer','Painters')
            suptitle([d2(id).name(1) ' ' d2(id).name(3:end-4) ', Cell' num2str(other_cells(icell==other_cells)) ' Raw'])
            helper_saveandclosefig([savefolder '\Figure3\related\BarPlots_ByArm_2DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(other_cells(icell==other_cells)) figlab])                     
            
            figure; hold on
            smdat = NaN(size(OpenFRinternal,1),size(OpenFRinternal,2),2,2);
            for ivel = 1:2
                for ilocal = 1:2            
                    dat = OpenFRinternal(:,:,icell,ivel,ilocal);
                    bb(ivel,ilocal) = subplot(2,2,(ivel-1)*2+ilocal);
                    datnan = isnan(dat);
                    dat(datnan) = 0;
    %                 dat = fillmissing(dat,'nearest');
                    dat = filter2(Two_D_Filter,dat);
                    dat2 = dat; dat2(datnan) = NaN;
                    smdat(:,:,ivel,ilocal) = dat2;
                    smx = max(max(dat));
                    dat(datnan) = smx+3*range(range(dat))/size(cm2,1);
    %                 dat(datnan) = NaN;
                    imagesc(dat); 
    %                 smdat(:,:,ivel,ilocal) = dat;
                    set(gca,'colormap',cm2)
    %                 colormap jet
                    axis xy
                    cmi = get(gca,'clim');
    %                 set(gca,'clim',[0 cmi(2)])
                    colorbar('Ticks',[min(min(dat)) max(max(dat))],'TickLabels',[round(min(min(dat)),2,'significant'),round(max(max(dat)),2,'significant')])
    %                 title([d2(id).name(1:end-4) ', Cell ' num2str(icell) ' ' cellid ', ' vellab{ivel} ' ' loclab{ilocal}])
                    title([vellab{ivel} ' ' loclab{ilocal}])

                if ivel==2
                dat1 = OpenFRinternal(:,:,icell,1,2);
                dat2 = OpenFRinternal(:,:,icell,2,ilocal);
                [r,p] = corr(dat1(:),dat2(:),'rows','complete');


                %corr to low velocity non-local - smoothed
                dat1 = smdat(:,:,1,2);
                dat2 = smdat(:,:,2,ilocal); %high velocity local
        %         dat3 = smdat(:,:,2,2); %high velocity non-local
                [r2,p2] = corr(dat1(:),dat2(:),'rows','complete');        
    %         r2 = corr(dat1(:),dat3(:),'rows','complete');        
                xlabel(['r = ' num2str(round(r,2,'significant')) ',p = ' num2str(round(p,2,'significant')) ', '  ...
                    'r2 = ' num2str(round(r2,2,'significant')) ',p2 = ' num2str(round(p2,2,'significant')) ])
                end     
                end
            end        


            suptitle([d2(id).name(1:end-4) ', Cell ' num2str(icell) ' ' cellid])
            set(gcf,'Position',[2134          -6        1123         780])
            helper_saveandclosefig([savefolder '\Figure3\related\2DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(icell) figlab])                     
        end
       
    end
    
        ttt = toc;
    disp(['Done with day ' num2str(id) ' in ' num2str(ttt/60) ' min'])
end    
%%
figure; hold on;
subplot(2,2,[1 2])
bar(numbinsused*EstBin/60)
xlabel('Session'); ylabel('Minutes')
legend('Used','Excluded')
subplot(2,2,3)
bar(mean(numbinsused*EstBin/60))
yl = get(gca,'ylim');
text(.5,yl(2)/2,[num2str(100*mean(numbinsused(:,1))./mean(numbinsused(:,2))) '%'])
ylabel('Average Minutes')
set(gca,'xticklabel',{'Used';'Excluded'})
subplot(2,2,4); hold on
bar(1,min(numbinsused(:,1)*EstBin/60))
text(.6,min(numbinsused(:,1)*EstBin/60)*2,num2str(min(numbinsused(:,1)*EstBin/60)))
bar(2,max(numbinsused(:,1)*EstBin/60))
text(1.6,max(numbinsused(:,1)*EstBin/60)/2,num2str(max(numbinsused(:,1)*EstBin/60)))
bar(3,mean(numbinsused(:,1)*EstBin/60))
text(2.6,mean(numbinsused(:,1)*EstBin/60)/2,num2str(mean(numbinsused(:,1)*EstBin/60)))
ylabel('Minutes')
set(gca,'xtick',1:3,'xticklabel',{'Minimum','Maximum','Average'})
helper_saveandclosefig([savefolder '\Figure3\UsedData' figlab])


figure; hold on; 
h = histogram(toterr(:,1)*binsize/pix2cm,'FaceColor','k','LineWidth',1,'Normalization','probability'); 
yl = get(gca,'ylim');
hold on; plot([.5*binsize/pix2cm 1*binsize/pix2cm],[0 yl(2)*1.1],'-r','LineWidth',.5); text(1*binsize/pix2cm,yl(2)*1.1,[num2str(1*binsize/pix2cm) 'cm'],'Color','r','FontSize',14)
hold on; plot([1.5*binsize/pix2cm 1.5*binsize/pix2cm],[0 yl(2)*1.05],'-g','LineWidth',.5);text(1.5*binsize/pix2cm,yl(2)*1.05,[num2str(1.5*binsize/pix2cm) 'cm'],'Color','g','FontSize',14)
hold on; plot([2.5*binsize/pix2cm 2.5*binsize/pix2cm],[0 yl(2)],'-b','LineWidth',.5);text(2.5*binsize/pix2cm,yl(2)*1,[num2str(2.5*binsize/pix2cm) 'cm'],'Color','b','FontSize',14)
hold on; plot([4.5*binsize/pix2cm 4.5*binsize/pix2cm],[0 yl(2)*.95],'-k','LineWidth',.5);text(4.5*binsize/pix2cm,yl(2)*.95,[num2str(4.5*binsize/pix2cm) 'cm'],'Color','k','FontSize',14)
ylim([yl(1) yl(2)*1.14])
xlabel(['Error (cm)'])
ylabel(['Probability of decoded error'])
set(gca,'FontSize',18)
helper_saveandclosefig([savefolder '\Figure3\ErrorAllSessions_HighVelOnly_binsize' num2str(binsize)])
% load([savefolder '\igroup'  num2str(igroup) figlab '.mat'])
%% difference plots

zlab = {'Across Cells';'Median Subtracted';'Zscored';'Norm'};
siglab = {'All Cells';'Sig Cells'};
vellab = {'High Vel';'Low Vel'};
ind = eye(3)==1;
for ivel = 1:2
    if ivel==1
        othall2 = othall; othallS2 = othallS;
    elseif ivel == 2
        othall2 = othall_lowvel;
    end
    for iz = 4:-1:3
        if iz == 1
            r = repall;
            o = othall2;
        elseif iz == 2
            r = repall-median(repall,2);
            o = othall2-median(othall2,2);
        elseif iz == 3
            r = zscore(repall,[],2);
            o = zscore(othall2,[],2);
            if ivel==1
                oS = zscore(othallS2,[],2);
            end
        elseif iz == 4
            r = (repall-min(repall,[],2))./range(repall,2);
            o = (othall2-min(othall2,[],2))./range(othall2,2);
            if ivel==1
                oS = (othallS2-min(othallS2,[],2))./range(othallS2,2);
            end
        end

        for isig = 2:-1:1
           if isig == 1
               indd = sum(repall<0,2)~=3;
               indd = true(size(indd));               
           elseif isig == 2
%                indd = sigcells(:,1)<.2 & sum(repall<0,2)~=3;
               indd = sigcells(:,1)<.05;               
           end
           rr = r(indd,:); oo = o(indd,:,:);
           S = sigcells(indd,:);               
%                S = max(repall(indd,:),2);
           if ivel==1; oos = oS(indd,:,:,:);end
                
            nonl = oo(:,:,3); l = oo(:,:,2);             
            nonl2 = repmat(nonl,[1 1 3]); l2 = repmat(l,[1 1 3]);
            rip = reshape(rr,[size(rr,1) 1 size(rr,2)]);
            rip2 = repmat(rip,[1 3 1]);
            dl = abs(l2-rip2);
            dnl = abs(nonl2-rip2);
            NL = sum(dnl(:,ind),2)./sum(sum(dnl,2),3);
            L = sum(dl(:,ind),2)./sum(sum(dl,2),3);                        
            dif = L-NL;
            [rho,p] = corr(S(dif~=0,1),dif(dif~=0),'rows','complete');
            
            
            figure; hold on
% %             b1 = bar(1,nanmean(L),'k','EdgeColor','k');
% %             b2 = bar(2,nanmean(NL),'k','EdgeColor','k');
% %             b2.FaceColor = 'w'; b1.FaceColor = 'w';
            violinplot1([L NL],1:2);
            errorbar(1,nanmean(L),nanstd(L)./sqrt(sum(~isnan(L))),'k','LineWidth',3)
            plot(1,nanmean(L),'k.','MarkerSize',30)
            errorbar(2,nanmean(NL),nanstd(NL)./sqrt(sum(~isnan(NL))),'k','LineWidth',3)
            plot(2,nanmean(NL),'k.','MarkerSize',30)
            
%             boxplot2([L NL],[1 2])
%             boxplot2(L,1)
%             boxplot2(NL,2)
            
            
            set(gca,'xtick',1:2,'xticklabel',{'Local';'NonLocal'})
            [p2,~,stats] = signrank(L,NL,'tail','right');
            if isig ==1
                if isfield('stats','zval')
                    title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ', rho = ' num2str(rho) ' p = ' num2str(p) ' z = ' num2str(stats.zval) ' n = ' num2str(length(L))])    
                else
                    title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ', rho = ' num2str(rho) ' p = ' num2str(p) ' SR = ' num2str(stats.signedrank) ' n = ' num2str(length(L))])    
                end
            else
                if isfield('stats','zval')
                    title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ' z = ' num2str(stats.zval) ' n = ' num2str(length(L))])    
                else
                    title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ' SR = ' num2str(stats.signedrank) ' n = ' num2str(length(L))])    
                end
            end
            set(gcf,'renderer','Painters')
            helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_NL_' ...
                   siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab '_nobar'])
               
            frr = oo(:,:,1);
            frr2 = repmat(frr,[1 1 3]);                        
            dfr = abs(frr2-rip2);
            B = sum(dfr(:,ind),2)./sum(sum(dfr,2),3);                       
            
            
            figure; hold on
            bar(1,nanmean(L),'w','EdgeColor','k');
            bar(2,nanmean(NL),'w','EdgeColor','k');
            bar(3,nanmean(B),'w','EdgeColor','k');
%             b2.FaceColor = 'w'; b1.FaceColor = 'w'; b1.FaceColor = 'w';
            errorbar(1,nanmean(L),nanstd(L)./sqrt(sum(~isnan(L))),'k','LineWidth',3)
            errorbar(2,nanmean(NL),nanstd(NL)./sqrt(sum(~isnan(NL))),'k','LineWidth',3)
            errorbar(3,nanmean(B),nanstd(B)./sqrt(sum(~isnan(B))),'k','LineWidth',3)
            set(gca,'xtick',1:3,'xticklabel',{'Local';'NonLocal';'Behavior'})
            [p2,~,stats] = signrank(B,NL,'tail','right');
            p3 = kruskalwallis([L NL B],[],'off');            
            if isfield('stats','zval')
                title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2)  ', pall = ' num2str(p3) ' z = ' num2str(stats.zval) ' n = ' num2str(length(L))])    
            else
                title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2)  ', pall = ' num2str(p3) ' SR = ' num2str(stats.signedrank) ' n = ' num2str(length(L))])    
            end
            
            helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_NL_BehFR_' ...
                   siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])
               
            if ivel==1
                  nonlS = oos(:,:,3,:); 
%                   lS = oos(:,:,2,:);             
                  nonlS2 = repmat(nonlS,[1 1 3 1]); 
%                   lS2 = repmat(lS,[1 1 3 1]);                  
                ripS2 = repmat(rip2,[1 1 1 size(nonlS,4)]);
%                 dlS1 = abs(lS2-ripS2);
%                 dlS = reshape(dlS1,[size(lS2,1) size(lS2,2)*size(lS2,3) size(lS2,4)]);
                dnlS1 = abs(nonlS2-ripS2);
                dnlS = reshape(dnlS1,[size(nonlS2,1) size(nonlS2,2)*size(nonlS2,3) size(nonlS2,4)]);
                NLS = squeeze(sum(dnlS(:,ind(:),:),2))./squeeze(sum(sum(dnlS1,2),3));
%                 LS = squeeze(sum(dlS(:,ind(:),:),2))./squeeze(sum(sum(dlS1,2),3));                         
%                 difS = squeeze(LS-NLS);
                difS = squeeze(L-NLS);
                
                figure; hold on;
                histogram(nanmean(difS),20,'FaceColor','w')
                yl = get(gca,'ylim');
                plot([nanmean(dif) nanmean(dif)],yl,'k','LineWidth',3)
                p= (sum(nanmean(difS)>=nanmean(dif))+1)./(size(difS,2)+1);
                title(['Shuffle ' zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p)])    
                
            set(gcf,'renderer','Painters')
                helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_NL_Shuffle_' ...
                   siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])
               
               
               diffcellsp = (sum(difS>=repmat(dif,[1 size(difS,2)]),2)+1)./(1+size(difS,2));
               pb = 1-binocdf(sum(diffcellsp<.05)-1,sum(~isnan(diffcellsp)),.05);
               srp = signrank(dif,mean(difS,2),'tail','right');
               
               figure; hold on
               subplot(2,1,1); hold on
               histogram(diffcellsp,0:.05:1,'FaceColor','w','LineWidth',2)
               yl = get(gca,'ylim');
               plot([.05 .05],[yl(1) yl(2)+.5],'r--','LineWidth',3)
               title(['Binomial test, ' num2str(sum(diffcellsp<.05)) ' of ' num2str(sum(~isnan(diffcellsp))) ', p = ' num2str(round(pb,2,'significant'))])
               ylabel('Sig PFC Neurons')
               xlabel('Permutation test p-value')               
               set(gca,'FontSize',18)
               
               subplot(2,1,2); hold on
               plot([mean(difS,2) dif]','color',[.5 .5 .5])
               errorbar(1:2,mean([mean(difS,2) dif]),std([mean(difS,2) dif])./sqrt(length([mean(difS,2) dif])),'k','LineWidth',3)
               title(['One-sided signed rank test, p = ' num2str(round(srp,2,'significant'))])
               set(gca,'xtick',1:2,'xticklabel',{'Shuffle';'Real'})
               ylabel('Difference score')
               xlim([.7 2.3])
               set(gca,'FontSize',18)
               set(gcf,'Position',[  2118          22         554         837])
                helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_NL_Shuffle_CellByCell_' ...
                   siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])
            end
            
        end    
    end
end


%% group data
  RS = []; SSD = []; SIG = []; FR = []; FRs = [];
for id = 1:size(d2,1)
    clear OpenFRinternal
%     load(d2(id).name,'spikedata','RP_SSDarm','hp_cells','other_cells','hpinterneurons','pos','OpenFRinternal')  
    load(d2(id).name,'RP_pSSDarm','RP_SSDarm','OpenFRinternal','other_cells','other_cells_touse')

    other_cells = other_cells(other_cells_touse(:,igroup));
    RP_pSSDarm = RP_pSSDarm(other_cells_touse(:,igroup));
    
    sig = RP_pSSDarm<.05;
%     sig = RP_pSSDarm<.1;
    
    rs = NaN(length(other_cells),2,2); % cells, local non-local, raw smoothed
    fr1 = NaN(length(other_cells),size(OpenFRinternal,1)*size(OpenFRinternal,2),2,2);
    fr2 = fr1;
    for icell = 1:length(other_cells)
        
        
        smdat = NaN(size(OpenFRinternal,1),size(OpenFRinternal,2),2,2);
        for ilocal = 1:2
            for ivel = 1:2            
                dat = OpenFRinternal(:,:,other_cells(icell),ivel,ilocal);
                fr1(icell,:,ivel,ilocal) = dat(:);
                datnan = isnan(dat);
%                 dat(datnan) = 0;
                dat = fillmissing(dat,'nearest');
                dat = filter2(Two_D_Filter,dat);                
                dat(datnan) = NaN;
                smdat(:,:,ivel,ilocal) = dat;
                fr2(icell,:,ivel,ilocal) = dat(:);
            end     
        end                
        
        %corr to low velocity non-local
        dat1 = OpenFRinternal(:,:,other_cells(icell),1,2);
        dat2 = OpenFRinternal(:,:,other_cells(icell),2,1); %high velocity local
        dat3 = OpenFRinternal(:,:,other_cells(icell),2,2); %high velocity non-local
        rs(icell,1,1) = corr(dat1(:),dat2(:),'rows','complete','type','Spearman');        
        rs(icell,2,1) = corr(dat1(:),dat3(:),'rows','complete','type','Spearman');        
        
        %corr to low velocity non-local - smoothed
        dat1 = smdat(:,:,1,2);
        dat2 = smdat(:,:,2,1); %high velocity local
        dat3 = smdat(:,:,2,2); %high velocity non-local
        rs(icell,1,2) = corr(dat1(:),dat2(:),'rows','complete','type','Spearman');        
        rs(icell,2,2) = corr(dat1(:),dat3(:),'rows','complete','type','Spearman');        
        
        
    end

    SIG = cat(1,SIG,sig);
    SSD = cat(1,SSD,RP_SSDarm);
    RS = cat(1,RS,rs);
%     FR = cat(1,FR,fr1);
%     FRs = cat(1,FRs,fr2);
end
    
%% difference between high velocity local and non-local with low velocity non-local

%raw 

figure; hold on
subplot(2,1,1); hold on
histogram(RS(:,2,1)-RS(:,1,1),10,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(:,2,1)-RS(:,1,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.9 .9])
[p,~,stats] = signrank(RS(:,2,1),RS(:,1,1),'tail','right');
if isfield('stats','zval')
    title(['All Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(size(RS,1))])
else
    title(['All Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(size(RS,1))])
end

subplot(2,1,2); hold on
histogram(RS(SIG==1,2,1)-RS(SIG==1,1,1),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,1)-RS(SIG==1,1,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.9 .9])
[p,~,stats] = signrank(RS(SIG==1,2,1),RS(SIG==1,1,1),'tail','right');
if isfield('stats','zval')
    title(['Sig Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(sum(SIG==1))])
else
    title(['Sig Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(sum(SIG==1))])
end
set(gcf,'Position',[ 2160         289         361         487])
set(gcf,'renderer','Painters')
if ~isfolder([savefolder '\Figure3\related\'])
    mkdir([savefolder '\Figure3\related\'])
end
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Raw' figlab])
%%
% smoothed
figure; hold on

subplot(2,1,1); hold on
histogram(RS(:,2,2)-RS(:,1,2),10,'facecolor','k')
xlim([-1 1])
yl = get(gca,'ylim');
md = nanmean(RS(:,2,2)-RS(:,1,2));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
set(gca,'ylim',yl)
% xlim([-.3 .3])
[p,~,stats] = signrank(RS(:,2,2),RS(:,1,2),'tail','right');
if isfield('stats','zval')
title(['All Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(size(RS,1))])
else
title(['All Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(size(RS,1))])
end
subplot(2,1,2); hold on
histogram(RS(SIG==1,2,2)-RS(SIG==1,1,2),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,2)-RS(SIG==1,1,2));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-1 1])
[p,~,stats] = signrank(RS(SIG==1,2,2),RS(SIG==1,1,2),'tail','right');
if isfield('stats','zval')
    title(['Sig Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(sum(SIG==1))])
else
    title(['Sig Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(sum(SIG==1))])
end
set(gcf,'Position',[ 2160         289         361         487])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Smoothed' figlab])

%% difference between 0 and high velocity non-local with low velocity non-local

%raw 
figure; hold on
subplot(2,1,1); hold on
histogram(RS(:,2,1),10,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(:,2,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.65 .65])
[p,~,stats] = signrank(RS(:,2,1),0,'tail','right');
if isfield('stats','zval')
title(['All Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(size(RS,1))])
else
    title(['All Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(size(RS,1))])
end

subplot(2,1,2); hold on
histogram(RS(SIG==1,2,1),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.65 .65])
[p,~,stats] = signrank(RS(SIG==1,2,1),0,'tail','right');
if isfield('stats','zval')
title(['Sig Cells: ' num2str(p) ', z = ' num2str(stats.zval) ', n = ' num2str(sum(SIG==1))])
else
    title(['All Cells: ' num2str(p) ', SR = ' num2str(stats.signedrank) ', n = ' num2str(size(RS,1))])
end
set(gcf,'Position',[ 2160         289         361         487])
set(gcf,'renderer','Painters')

helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Raw_diff0' figlab])
%%
% smoothed

figure; hold on
subplot(2,1,1); hold on
histogram(RS(:,2,2),10,'facecolor','k')
yl = get(gca,'ylim');
md = nanmedian(RS(:,2,2));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-1 1])
set(gca,'ylim',yl)
p = signrank(RS(:,2,2),0,'tail','right');
title(['All Cells: ' num2str(p)])
subplot(2,1,2); hold on
histogram(RS(SIG==1,2,2),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmedian(RS(SIG==1,2,2));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-1 1])
p = signrank(RS(SIG==1,2,2),0,'tail','right');
title(['Sig Cells: ' num2str(p)])
set(gcf,'Position',[ 2160         289         361         487])
set(gcf,'renderer','Painters')
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Smoothed_diff0' figlab])
