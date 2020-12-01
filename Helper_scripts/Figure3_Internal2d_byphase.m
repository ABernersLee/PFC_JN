function Figure3_Internal2d_byphase(dirs,igroup,savefolder,toplotcells,tosave,velcutoff,cellcutoff,spikecutoff,numshuff)

% newcol = [75 0 130;255 130 0;34 139 34]/255;
cd(dirs.homedir)
d2 = dir('*.mat');
binsize = 20;
SigmaField = 2;
noreplacement = false;
Two_D_Filter=fspecial('gaussian',[3 3],SigmaField);
EstBin = .06; 
vellab = {'Low Velocity';'High Velocity'};
loclab = {'Local';'Non-Local'};
if ~isfolder([savefolder '\Figure3\related\'])
    mkdir([savefolder '\Figure3\related\'])
end
figlab = ['_Theta_Phase_BinSize' num2str(EstBin) '_VelCutoff' num2str(velcutoff) '_cellcutoff' num2str(cellcutoff) '_spikecutoff' num2str(spikecutoff) '_numshuff' num2str(numshuff)];
dayran = [];
repall = []; sigcells = []; othall = []; othall_lowvel = []; othallS = []; errbyphase = [];
% tocalc = true;
for id = 1:size(d2,1)
    load(d2(id).name,'other_cells_touse')
    if sum(other_cells_touse(:,igroup))==0
        continue
    end
    clear OpenFRinternal
    tocalc = true;
    if tocalc % get the internal fields        
        
        skipped = 0; notskipped = 0;
        
        load(d2(id).name,'OpenFR','hp_cells','hpinterneurons','spikedata','pos','armpos','armposindex','vel','other_cells','other_cells_touse','RP_CandEventTimes','times_armon_thetaof_headingarm_lap_thetahalf_all');


        Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
        clear times_armon_thetaof_headingarm_lap_thetahalf_all
        Event = Th(:,2); % triggered off the end of theta


        Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));

        for icell = 1:size(OpenFR,3)
            dat = OpenFR(:,:,icell);                
            datnan = isnan(dat);
            dat(datnan) = 0;
            dat = filter2(Two_D_Filter,dat);            
            dat(datnan) = NaN;
            OpenFR(:,:,icell) = dat;
        end

        FR1 = reshape(OpenFR(:,:,Cell_Number),[size(OpenFR,1)*size(OpenFR,2) length(Cell_Number)])';        
        nanFR = sum(~isnan(FR1))~=0;
        nanIND = find(nanFR);
        FR = FR1(:,nanFR);

        Cand = sortrows([Event-.12 Event-.6; Event-.6 Event]); %first half, second half            
        B=cumsum([0 diff(Cand')]);


        Spike = [spikedata(:,2) spikedata(:,1)];
        [~,I]=histc(Spike(:,2),sortrows(Cand(:)));    
        for i=1:2:max(I)
            Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
        end
        Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
        % clear Cand I i

        Start_Time=0;
        End_Time=size(Cand,1)*EstBin;
        TimeBins=round((End_Time-Start_Time)/EstBin);
        times = mean(Cand,2);

        % Bayesian Decoding - EstBin ms window estimate
        binspike=zeros(length(Cell_Number),TimeBins);
        for CellID=1:length(Cell_Number)   
            c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time:EstBin:End_Time);
            binspike(CellID,1:TimeBins)=c(1:TimeBins);                            
        end
        firsthalf = repmat([1 0],[1 size(Cand,1)/2])';
        FR(sum(FR~=0,2)==0,:) = 1;
        FR(FR==0 | FR<0) = .0001;

        excl = true(length(times),1);
        ripdat = sortrows(RP_CandEventTimes(:,1:2))';
        [~,i] = histc(times,ripdat(:));
        excl(~mod(i,2),:) = false;
        clear i


        touse = sum(binspike>0)>=cellcutoff; % has to have at least 5 cells spiking
        touse = touse & sum(binspike)>=spikecutoff; %has to have at least x spikes 
        touse = touse & ~excl'; %and can't be in a ripple event
        binspike = binspike(:,touse);
        times = times(touse);
        TimeBins = length(times);
        firsthalf = firsthalf(touse);

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
        % Index=round(B/(EstBin));

        %normalize for each time bin
        Est = Matrix./(ones(size(Matrix,1),1)*sum(Matrix,1));
        % %make times when there is zero spiking equal to zeros/NaNs all the way down
        % Matrix(:,sum(binspike>0)==0) = NaN;

        exclude = false(length(times),1);

        [~,i] = histc(times,pos(:,1));
        inx = find(i==0);
        exclude(inx) = true;
        i(inx) = [];
%             inx2 = inx<floor(length(times)/2);
%             i(inx(inx2)) = i(find(i~=0,1,'first')); i(inx(~inx2)) = i(find(i~=0,1,'last'));


        pos2 = pos(:,2:3);
        pos2(:,1) = ceil(pos2(:,1)/binsize);
        pos2(:,2) = ceil(pos2(:,2)/binsize);
        pos2 = pos2-min(pos2)+1;

        binposition = pos2(i,:);
        binvel = vel(i);
        binarm = armpos(i);

         [~,decoded_position1] = max(Est); % change from max to mean and do seperately for each arm, then take the max                       
        [X,Y] = ind2sub(max(pos2),nanIND(decoded_position1)');
        decoded_position2 = [X Y];

%             t = times(~exclude)';
        dp = decoded_position2(~exclude,:);
        bp = binposition(~exclude,:);
        bv = binvel(~exclude);

        %arm rat truely is on
        firsthalf1 = firsthalf(~exclude);        
        err = vecnorm(dp-bp,2,2);
        errbyphase = cat(1,errbyphase,[err firsthalf1 id*ones(size(err,1),1)]);
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

         TimeBins=round((End_Time-Start_Time)/EstBin);
         Cell_Number = 1:max(Spike(:,1));
         Start_Time = Start_Time+EstBin;    %added delay of 60ms
         binspike2=zeros(length(Cell_Number),TimeBins);                             
        for CellID=1:length(Cell_Number)   
            c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time:EstBin:End_Time+EstBin);
            binspike2(CellID,1:TimeBins)=c(1:TimeBins);                            
        end
            bsp3 = binspike2(:,touse);
            bsp3 = bsp3(:,~exclude);
            
        ind = 1:pind(1)*pind(2);
        pos4 = sub2ind(pind,dp(:,1),dp(:,2));
%         figure; histogram(err(ismember(pos4,ind)))
%         title(num2str(median(err(ismember(pos4,ind)))))
        OpenFRinternal = NaN(pind(1),pind(2),size(OpenFR,3),2,2);
        OpenFRinternalS = NaN(pind(1),pind(2),size(OpenFR,3),2,2,numshuff);
        ivel = 2;
        if ivel==1; velind = bv<velcutoff; elseif ivel==2; velind = bv>velcutoff; end
        pos3 = pos4;
        pos3(~velind,:) = NaN; %restrict to times with the correct velocity and locality

        occ = NaN(size(ind,2),2);
        occS = NaN(size(ind,2),2,numshuff);
        f = zeros(size(ind,2),size(OpenFR,3),2);
        fS = zeros(size(ind,2),size(OpenFR,3),2,numshuff);
        for j = 1:length(ind)
%                   
            skipit = false;
            num1 = sum(pos3==ind(j) & firsthalf1==1);
            num2 = sum(pos3==ind(j) & firsthalf1==0);
            h = [num1 num2];
            numperm = nchoosek(sum(h),min(h));
            if noreplacement                
                
                if num1==0 || num2==0
%                     skipped = skipped+1;
    %                 continue
                    skipit = true;
                end

                if numperm<numshuff     
%                     skipped = skipped+1;
    %                 continue
                    skipit = true;
                end
                
            else
                if (num1<1 || num2<1) || numperm<50 
                    skipit = true;
                end
            end
            
            if skipit
                skipped = skipped+1;
                continue
            end
%             disp(num2str([num1 num2 numperm]))
            notskipped = notskipped+1;

            for ilocal = 1:2      
                if ilocal == 1; mind =firsthalf1==1; elseif ilocal ==2; mind = firsthalf1==0; end
                occ(ind(j),ilocal) = sum(pos3==ind(j) & mind)*EstBin;
                spkb = bsp3(:,pos3==ind(j) & mind);
                f(ind(j),:,ilocal) = sum(spkb,2);

                if rem(j,1000)==0
                    disp(['j = ' num2str(j)])
                end
            end

             errpos = firsthalf1(pos3==ind(j));
            for ishuff = 1:numshuff
                errs = firsthalf1;
                errs(pos3==ind(j)) = errpos(randperm(length(errpos)));       
                
                if noreplacement
                    if sum((errs==1)~=(firsthalf1==1))==0
                        while sum((errs==1)~=(firsthalf1==1))==0
                            errs(pos3==ind(j)) = errpos(randperm(length(errpos)));
    %                         disp(['Same ' num2str(ishuff)])
                        end
                    end

                    if ishuff ==1
                        errsshuffsave = errs;
                    else
                        while any(sum((errsshuffsave==1)~=(errs==1))==0)
                            errs(pos3==ind(j)) = errpos(randperm(length(errpos)));
    %                         disp(['Same ' num2str(ishuff)])
                        end
                          errsshuffsave = cat(2,errsshuffsave,errs);
                    end
                end

                for ilocal = 1:2
                    if ilocal == 1; mindS = errs==1; elseif ilocal ==2; mindS = errs==0; end
                    occS(ind(j),ilocal,ishuff) = sum(pos3==ind(j) & mindS)*EstBin;
                    spkbS = bsp3(:,pos3==ind(j) & mindS);
                    fS(ind(j),:,ilocal,ishuff) = sum(spkbS,2);
                end
            end
            clear errsshuffsave

        end
        occ2 = reshape(occ,[pind 2]);     
        occS2 = reshape(occS,[pind 2 numshuff]);     
        if toplotcells
            figure; hold on
            subplot(3,1,1)
            imagesc(occ2(:,:,1)/EstBin); axis xy; colorbar
            title([vellab{ivel} ' ' loclab{1}])
            subplot(3,1,2)
            imagesc(occ2(:,:,2)/EstBin); axis xy; colorbar
            title([vellab{ivel} ' ' loclab{2}])
            subplot(3,1,3)
            imagesc(occ2(:,:,2)/EstBin-occ2(:,:,1)/EstBin); axis xy; colorbar
            title('Difference in Occupancy')
            set(gcf,'Position',[2562           6         494         879])
            helper_saveandclosefig([savefolder '\Figure3\related\2DInternal_ThetaPhase_' vellab{ivel} '_' d2(id).name(1:end-4) '_Occ']) 
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

    end
    
    if tosave
%             save(d2(id).name,'OpenFRinternal','-append')        
    else
% %         load(d2(id).name,'OpenFRinternal')
    end
    
    load(d2(id).name,'other_cells','other_cells_touse')        
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
        
        
    load(d2(id).name,'pos','armpos')
    pos2 = pos(:,2:3);
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
    
    
    fr = NaN(size(OpenFRPFC,3),3);
    for icell = 1:size(fr,1)
        for iarm = 1:3
            dat = OpenFRPFC(:,:,icell);
            dat = dat(x2);
            fr(icell,iarm) = nanmean(dat(aa==iarm));
        end
    end
            
%             figure; hold on; for iarm = 1:3; plot(posu(aa==iarm,2),posu(aa==iarm,1),'.'); end; legend

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
    
    
    load(d2(id).name,'other_cells','RP_moduarm','RP_pSSDarm','RP_SSDarm','other_cells_touse','hp_cells','hpinterneurons')
    replay = RP_moduarm';
    repall = cat(1,repall,replay(other_cells_touse(:,igroup),:));
    
    local = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,1,1);
    nonlocal = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,2,1);    
    othall = cat(1,othall,cat(3,fr(other_cells_touse(:,igroup),:),local(other_cells_touse(:,igroup),:),nonlocal(other_cells_touse(:,igroup),:)));
    
    localS = ArmBarPlot_LocalNonLocal_MeanSEMS(other_cells,:,1,1,:);
    nonlocalS = ArmBarPlot_LocalNonLocal_MeanSEMS(other_cells,:,2,1,:);    
    othallS = cat(1,othallS,cat(3,NaN(sum(other_cells_touse(:,igroup)),size(fr,2),1,size(localS,5)),localS(other_cells_touse(:,igroup),:,:,:),nonlocalS(other_cells_touse(:,igroup),:,:,:)));
    
    
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
%         helper_saveandclosefig([savefolder '\Figure3\related\Internal_' d2(id).name(1:end-4)])

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
%             subplot(2,2,2); hold on
%             bar(nonlocal_L(icell==other_cells,:),'FaceColor','w')
            subplot(2,2,3); hold on
            bar(local(icell==other_cells,:),'FaceColor','w')
            subplot(2,2,4); hold on
            bar(nonlocal(icell==other_cells,:),'FaceColor','w')
            set(gcf,'renderer','Painters')
            suptitle([d2(id).name(1) ' ' d2(id).name(3:end-4) ', Cell' num2str(other_cells(icell==other_cells)) ' Raw'])
            helper_saveandclosefig([savefolder '\Figure3\related\BarPlots_ByArm_2DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(other_cells(icell==other_cells)) figlab])                     
            
            figure; hold on
            smdat = NaN(size(OpenFRinternal,1),size(OpenFRinternal,2),2,2);
            ivel = 2;
            for ilocal = 1:2            
                dat = OpenFRinternal(:,:,icell,ivel,ilocal);
                bb(ivel,ilocal) = subplot(2,1,ilocal);
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
            end                     

            suptitle([d2(id).name(1:end-4) ', Cell ' num2str(icell) ' ' cellid])
            set(gcf,'Position',[  2220          38         430         707])
            helper_saveandclosefig([savefolder '\Figure3\related\2DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(icell) figlab])                     
        end
       
    end
    
    disp(['Done with day ' num2str(id) ' skipped ' num2str(skipped) ', ran ' num2str(notskipped)])
    dayran = [dayran; id notskipped];
end    

disp(dayran)

% load([savefolder '\igroup'  num2str(igroup) figlab '.mat'])
%%
p1 = ranksum(errbyphase(errbyphase(:,2)==1,1),errbyphase(errbyphase(:,2)==0,1),'tail','left');
ebp = NaN(max(errbyphase(:,3)),2);
for id = 1:max(errbyphase(:,3))
    ebp(id,:) = [mean(errbyphase(errbyphase(:,2)==1 & errbyphase(:,3)==id,1)) mean(errbyphase(errbyphase(:,2)==0 & errbyphase(:,3)==id,1))];
end
[p2,~,stats] = signrank(ebp(:,1),ebp(:,2),'tail','left','method','approximate');
[p2] = signrank(ebp(:,1),ebp(:,2),'tail','left');

figure; hold on
subplot(1,2,1); hold on
histogram(errbyphase(errbyphase(:,2)==1,1),15) %,'FaceColor','k')
histogram(errbyphase(errbyphase(:,2)==0,1),15) %,'FaceColor','r')
title(['p = ' num2str(p1)])
ylabel('Count')
xlabel('Error')

subplot(1,2,2); hold on
plot(ebp','.-k','MarkerSize',10)
for ib = 1:2
errorbar(ib,mean(ebp(:,ib)),std(ebp(:,ib))./sqrt(size(ebp,1)),'k','LineWidth',3)
end
xlim([.8 2.2])
ylabel('Decoding Error')
xlabel('Theta Phase')
set(gca,'xtick',1:2,'xticklabel',{'First Half';'Second Half'})
title(['p = ' num2str(p2) ', z = ' num2str(stats.zval) ', = ' num2str(size(ebp,1))])
set(gcf,'Position',[2185         495         758         249])
helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_Error' figlab])
%%
% othlabel = {'Behavioral FR';'Local';'Non-Local'};
zlab = {'Across Cells';'Median Subtracted';'Zscored';'Norm'};
siglab = {'All Cells';'Sig Cells'};
vellab = {'High Vel';'Low Vel'};
ind = eye(3)==1;
ivel =1;
    
othall2 = othall; othallS2 = othallS;
    

for iz = 3:4
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

    for isig = 1:2
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
        [rho,p] = corr(S(dif~=0,2),dif(dif~=0),'rows','complete');


        figure; hold on
%         b1 = bar(1,nanmedian(L),'k','EdgeColor','k');
%         b2 = bar(2,nanmedian(NL),'k','EdgeColor','k');
%         b2.FaceColor = 'w'; b1.FaceColor = 'w';
         violinplot1([L NL],1:2);
            errorbar(1,nanmean(L),nanstd(L)./sqrt(sum(~isnan(L))),'k','LineWidth',3)
            plot(1,nanmean(L),'k.','MarkerSize',30)
            errorbar(2,nanmean(NL),nanstd(NL)./sqrt(sum(~isnan(NL))),'k','LineWidth',3)
            plot(2,nanmean(NL),'k.','MarkerSize',30)
          
%         errorbar(1,nanmedian(L),nanstd(L)./sqrt(sum(~isnan(L))),'k','LineWidth',3)
%         errorbar(2,nanmedian(NL),nanstd(NL)./sqrt(sum(~isnan(NL))),'k','LineWidth',3)
        set(gca,'xtick',1:2,'xticklabel',{'Local';'NonLocal'})
        [p2,~,stats] = signrank(L,NL,'tail','right');
        if isig ==1
            title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ', rho = ' num2str(rho) ' p = ' num2str(p) ' z = ' num2str(stats.zval) ' n = ' num2str(length(L))])     
        else
            title([zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2) ' z = ' num2str(stats.zval) ' n = ' num2str(length(L))])    
        end
        
            set(gcf,'renderer','Painters')
        helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_' ...
               siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])

        frr = oo(:,:,1);
        frr2 = repmat(frr,[1 1 3]);                        
        dfr = abs(frr2-rip2);
        B = sum(dfr(:,ind),2)./sum(sum(dfr,2),3);                       


        figure; hold on
        bar(1,nanmedian(L),'w','EdgeColor','k');
        bar(2,nanmedian(NL),'w','EdgeColor','k');
        bar(3,nanmedian(B),'w','EdgeColor','k');
%             b2.FaceColor = 'w'; b1.FaceColor = 'w'; b1.FaceColor = 'w';
        errorbar(1,nanmedian(L),nanstd(L)./sqrt(sum(~isnan(L))),'k','LineWidth',3)
        errorbar(2,nanmedian(NL),nanstd(NL)./sqrt(sum(~isnan(NL))),'k','LineWidth',3)
        errorbar(3,nanmedian(B),nanstd(B)./sqrt(sum(~isnan(B))),'k','LineWidth',3)
        set(gca,'xtick',1:3,'xticklabel',{'Local';'NonLocal';'Behavior'})
        p2 = signrank(B,NL,'tail','right');
        p3 = kruskalwallis([L NL B],[],'off');            
        title(['Theta Phase - ' zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p2)  ', pall = ' num2str(p3)])    
        set(gcf,'renderer','Painters')
        helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_BehFR_' ...
               siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])

        if ivel==1
              nonlS = oos(:,:,3,:); lS = oos(:,:,2,:);             
              nonlS2 = repmat(nonlS,[1 1 3 1]); lS2 = repmat(lS,[1 1 3 1]);                  
            ripS2 = repmat(rip2,[1 1 1 size(nonlS,4)]);
            dlS1 = abs(lS2-ripS2);
            dlS = reshape(dlS1,[size(lS2,1) size(lS2,2)*size(lS2,3) size(lS2,4)]);
            dnlS1 = abs(nonlS2-ripS2);
            dnlS = reshape(dnlS1,[size(nonlS2,1) size(nonlS2,2)*size(nonlS2,3) size(nonlS2,4)]);
            NLS = squeeze(sum(dnlS(:,ind(:),:),2))./squeeze(sum(sum(dnlS1,2),3));
            LS = squeeze(sum(dlS(:,ind(:),:),2))./squeeze(sum(sum(dlS1,2),3));                        
            difS = squeeze(LS-NLS);

            figure; hold on;
            histogram(nanmean(difS),20,'FaceColor','w')
            yl = get(gca,'ylim');
            plot([nanmean(dif) nanmean(dif)],yl,'k','LineWidth',3)
            p= (sum(nanmean(difS)>=nanmean(dif))+1)./(size(difS,2)+1);
            title(['Theta Phase - Shuffle ' zlab{iz} ' ' siglab{isig} ' - ' vellab{ivel} ' p = ' num2str(p)])    
            set(gcf,'renderer','Painters')
            helper_saveandclosefig([savefolder '\Figure3\BarPlots_ReplayVSLocalNonLocal_Shuffle_' ...
               siglab{isig} '_' zlab{iz} '_' vellab{ivel} figlab])
        end

    end    
end



%% group data
if 0
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
p = signrank(RS(:,2,1),RS(:,1,1),'tail','right');
title(['All Cells: ' num2str(p)])

subplot(2,1,2); hold on
histogram(RS(SIG==1,2,1)-RS(SIG==1,1,1),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,1)-RS(SIG==1,1,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.9 .9])
p = signrank(RS(SIG==1,2,1),RS(SIG==1,1,1),'tail','right');
title(['Sig Cells: ' num2str(p)])
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
p = signrank(RS(:,2,2),RS(:,1,2),'tail','right');
title(['All Cells: ' num2str(p)])
subplot(2,1,2); hold on
histogram(RS(SIG==1,2,2)-RS(SIG==1,1,2),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,2)-RS(SIG==1,1,2));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-1 1])
p = signrank(RS(SIG==1,2,2),RS(SIG==1,1,2),'tail','right');
title(['Sig Cells: ' num2str(p)])
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
p = signrank(RS(:,2,1),0,'tail','right');
title(['All Cells: ' num2str(p)])

subplot(2,1,2); hold on
histogram(RS(SIG==1,2,1),8,'facecolor','k')
yl = get(gca,'ylim');
md = nanmean(RS(SIG==1,2,1));
plot([md md],yl,'r-','LineWidth',2)
plot([0 0],yl,'k--')
xlim([-.65 .65])
p = signrank(RS(SIG==1,2,1),0,'tail','right');
title(['Sig Cells: ' num2str(p)])
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
end