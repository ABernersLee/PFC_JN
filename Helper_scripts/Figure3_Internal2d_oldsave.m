function Figure3_Internal2d_oldsave(dirs,igroup,savefolder,toplot,toplotcells,tosave)

% newcol = [75 0 130;255 130 0;34 139 34]/255;
cd(dirs.homedir)
d2 = dir('*.mat');
binsize = 20;
SigmaField = 2;
Two_D_Filter=fspecial('gaussian',[3 3],SigmaField);
velcutoff = 5;
EstBin = .06; 
cellcutoff = 3; %0; %3; %3; (3,5, and with err m as the median of each session (1.4-3.1)) not downsampled. m=2 is similar.
spikecutoff = 5; %0; %5; %5;
vellab = {'Low Velocity';'High Velocity'};
loclab = {'Local';'Non-Local'};
if ~isfolder([savefolder '\Figure3\related\'])
    mkdir([savefolder '\Figure3\related\'])
end


repall = []; sigcells = []; othall = [];
% tocalc = true;
for id = 1:size(d2,1)
    clear OpenFRinternal
%     load(d2(id).name,'spikedata','hp_cells','other_cells','hpinterneurons','pos','OpenFRinternal','RP_pSSDarm')                          
%     if exist('OpenFRinternal','var'); tocalc = false; else; tocalc = true; end
    tocalc = true;
    if tocalc && igroup == 1% get the internal fields        
        
        if 1 % get the error indicies
        
            load(d2(id).name,'RP_pSSDarm','OpenFR','hp_cells','hpinterneurons','spikedata','pos','armpos','armposindex','vel','other_cells','other_cells_touse')

            Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
            
            for icell = 1:size(OpenFR,3)
                dat = OpenFR(:,:,icell);                
                datnan = isnan(dat);
                dat(datnan) = 0;
                dat = filter2(Two_D_Filter,dat);
                smx = max(max(dat));
                dat(datnan) = NaN;
                OpenFR(:,:,icell) = dat;
            end

            FR1 = reshape(OpenFR(:,:,Cell_Number),[size(OpenFR,1)*size(OpenFR,2) length(Cell_Number)])';
            
            % FR2 = reshape(FR',[size(OpenFR,1) size(OpenFR,2) length(Cell_Number)]);
            nanFR = sum(~isnan(FR1))~=0;
            nanFR2 = reshape(nanFR,[size(OpenFR,1) size(OpenFR,2)]);
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
            % clear Cand I i

            Start_Time=0;
            End_Time=B(end)-EstBin;
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

            

            FR(sum(FR~=0,2)==0,:) = 1;
            FR(FR==0 | FR<0) = .0001;

            touse = sum(binspike>0)>=cellcutoff; % has to have at least 5 cells spiking
            touse = touse & sum(binspike)>=spikecutoff; %has to have at least x spikes
            binspike = binspike(:,touse);
            times = times(touse);
            TimeBins = length(times);

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

            [~,i] = histc(times,pos(:,1));
            inx = find(i==0);
            inx2 = inx<floor(length(times)/2);
            i(inx(inx2)) = i(find(i~=0,1,'first')); i(inx(~inx2)) = i(find(i~=0,1,'last'));

            
            pos2 = pos(:,2:3);
            pos2(:,1) = ceil(pos2(:,1)/binsize);
            pos2(:,2) = ceil(pos2(:,2)/binsize);
            pos2 = pos2-min(pos2)+1;

            binposition = pos2(i,:);
            binvel = vel(i);
            binarm = armpos(i);

            
            %excluding the very center of the maze. eh
            % for iarm = 1:3
            %     a = find(armposindex(:,iarm));
            %     exclude(binposition>=a(1) & binposition<=a(3)) = true;
            % end

            [mval,decoded_position1] = max(Est); % change from max to mean and do seperately for each arm, then take the max

            exclude = false(size(binposition,1),1);
%             exclude(mval<.2) = true;
            
            [X,Y] = ind2sub(max(pos2),nanIND(decoded_position1)');
            decoded_position2 = [X Y];

            t = times(~exclude)';
            dp = decoded_position2(~exclude,:);
            bp = binposition(~exclude,:);
            bv = binvel(~exclude);
                        
%             da = mm(~exclude);
%             diffarm = ba~=da;
            % toadd = [0 find(armposindex(:,1),1,'last') find(armposindex(:,2),1,'last')]';
            % err = abs(dp-bp);
%             figure; hold on; plot(bp(ba==1,1),bp(ba==1,2),'r.');plot(bp(ba==2,1),bp(ba==2,2),'b.'); plot(bp(ba==3,1),bp(ba==3,2),'k.')          

            %arm rat truely is on
            ba = binarm(~exclude);
            pind = max([dp; bp]);
            x = sub2ind(pind,bp(:,1),bp(:,2)); %real linear indicies
            y = sub2ind(pind,dp(:,1),dp(:,2)); %decoded linear indicies                                    
            ind1 = setdiff(unique(x(ba==1)),unique([x(ba==2);x(ba==3)]));
            ind2 = setdiff(unique(x(ba==2)),unique([x(ba==1);x(ba==3)]));
            ind3 = setdiff(unique(x(ba==3)),unique([x(ba==1);x(ba==2)]));
            
%             %arm rat is decoded to be on
%             da = NaN(size(y));    
%             da(ismember(y,ind1)) = 1;
%             da(ismember(y,ind2)) = 2;
%             da(ismember(y,ind3)) = 3;
            
            err = vecnorm(dp-bp,2,2);
%             err = diffarm;
        end
        
        
%         Spike = [spikedata(:,2) spikedata(:,1)];
%         [~,I]=histc(Spike(:,2),sortrows(Cand(:)+.1));    
%         for i=1:2:max(I)
%             Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
%         end
%         Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events
%         % clear Cand I i

%         Spike = [spikedata(:,2) spikedata(:,1)];
%         [~,I]=histc(Spike(:,2),sortrows(Cand(:)+EstBin));    
%         for i=1:2:max(I)
%             Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); %the spikes inside candidate events are at the times of those events
%         end
%         Spike(~mod(I,2),:)=[]; % only takes spikes within the candidate events, shifted by EstBin
%         % clear Cand I i
%             
% %         Start_Time = Start_Time+.1;
%          TimeBins=round((End_Time-Start_Time)/EstBin)*4;
%          Cell_Number = 1:size(OpenFR,3);
%          Start_Time = Start_Time+EstBin;
%          binspike2=zeros(length(Cell_Number),TimeBins);
% %          bintimes2 = NaN(TimeBins,1);
%             for i=1:4
%                 for CellID=1:length(Cell_Number)   
%                     c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
%                     binspike2(CellID,i:4:TimeBins)=c(1:TimeBins/4);
% %                     jnk = Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins;
% %                     if length(jnk)>length(i:4:TimeBins)
% %                         bintimes2(i:4:TimeBins) = jnk(1:end-1)+(EstBin/2);
% %                     else
% %                         bintimes2(i:4:TimeBins) = jnk+(EstBin/2);
% %                     end
%             %         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
%                 end                      
%             end
%             bsp3 = binspike2(:,touse);
%             bsp3 = bsp3(:,~exclude);
%         [~,~,tind] = histcounts(bintimes2+Cand(1),t);
        
%         bsp4 = bsp3(:,bv>velcutoff);
%         da4 = da(bv>velcutoff);
%         ba4 = ba(bv>velcutoff);
%         barplotarms = NaN(size(bsp4,1),3,2);
        ind = 1:pind(1)*pind(2);
        pos4 = sub2ind(pind,dp(:,1),dp(:,2));
%         figure; histogram(err(ismember(pos4,ind)))
%         title(num2str(median(err(ismember(pos4,ind)))))
        OpenFRinternal = NaN(pind(1),pind(2),size(OpenFR,3),2,2);
        for ivel = 1:2
            if ivel==1; velind = bv<velcutoff; elseif ivel==2; velind = bv>velcutoff; end
%             for ilocal = 1:2
%                 if ilocal ==1; locind = err<errcutoff; elseif ilocal ==2; locind = err>errcutoff; end
%                 pos4 = dp;
%                 pos4(~velind | ~locind,:) = NaN; %restrict to times with the correct velocity and locality
                pos3 = pos4;
                pos3(~velind,:) = NaN; %restrict to times with the correct velocity and locality
                     
                occ = NaN(size(ind,2),2);
                f = zeros(size(ind,2),size(OpenFR,3),2);
                for j = 1:length(ind)
                    m = nanmedian(err(pos3==ind(j)));                    
                    if m==0
                        m = .0001;
                    end           
                    
                    m = median(err(ismember(pos4,ind)));
%                     m = 2;
                    
                    if ~isnan(m)
%                         disp(m)
                        if sum(err>m)>sum(err<m)                            
                            mind1 = err<=m;
                            mind2 = err>m; 
                        elseif sum(err>m)<sum(err<m)                            
                            mind1 = err<m;
                            mind2 = err>=m;
                        elseif sum(err>m)==sum(err<m)                            
                            mind1 = err<m;
                            mind2 = err>m;
                        end
                            
                        
                        %downsample to equal 
%                         num1 = sum(pos3==ind(j) & mind1);
%                         num2 = sum(pos3==ind(j) & mind2);
%                         if (num2-num1)<0
%                             jnk = find(pos3==ind(j) & mind1);
%                             mind1(jnk(randperm(num1-num2))) = false;
%                         elseif (num2-num1)>0
%                             jnk = find(pos3==ind(j) & mind2);
%                             mind2(jnk(randperm(num2-num1))) = false;                            
%                         end
%                         
%                         num1 = sum(pos3==ind(j) & mind1);
%                         num2 = sum(pos3==ind(j) & mind2);
%                         if (num1-num2)~=0
%                             disp(num2str(num1-num2))
%                         end
% %                         
                        for ilocal = 1:2      
                            if ilocal == 1; mind = mind1; elseif ilocal ==2; mind = mind2; end
%                             if m==0
%                                 disp(j)
%                             end
                            occ(ind(j),ilocal) = sum(pos3==ind(j) & mind)*EstBin;

                            t2 = t(pos3==ind(j) & mind);   
                            spkb = bsp3(:,pos3==ind(j) & mind);
%                             occ(ind(j),ilocal) = sum(ismember(tind,find(pos3==ind(j) & mind)))*EstBin;
%                             spks = zeros(max(size(OpenFR,3)),length(t2));
%                             for itt = 1:length(t2)
%                                 s = spikedata(spikedata(:,1)>=t2(itt) & spikedata(:,1)<(t2(itt)+EstBin),2);
%                                 spks(:,itt) = histc(s,1:size(OpenFR,3));
%                             end
                            f(ind(j),:,ilocal) = sum(spkb,2);

                            if rem(j,1000)==0
                                disp(['j = ' num2str(j)])
                            end
                        end
                    end
                end
                occ2 = reshape(occ,[pind 2]);     
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
                
                OpenFRinternal_temp1 = f2(:,:,:,1)./repmat(occ2(:,:,1),[1 1 size(f2,3)]);
                OpenFRinternal_temp2 = f2(:,:,:,2)./repmat(occ2(:,:,2),[1 1 size(f2,3)]);
                
                OpenFRinternal(:,:,:,ivel,:) = cat(5,OpenFRinternal_temp1,OpenFRinternal_temp2);
                
%             end
        end
        if tosave
        save(d2(id).name,'OpenFRinternal','-append')
        end

    else
        load(d2(id).name,'OpenFRinternal')
    end
    
    if tosave && igroup ==1
        
        ArmBarPlot_LocalNonLocal_MeanSEM = NaN(size(OpenFRinternal,3),3,2,2);
        for icell = 1:size(OpenFRinternal,3)
            dat = OpenFRinternal(:,:,icell,2,1);
            dat = dat(:);
            ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,1,1) = [nanmean(dat(ind1)) nanmean(dat(ind2)) nanmean(dat(ind3))];
            ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,1,2) = [nanstd(dat(ind1))./sqrt(sum(~isnan(dat(ind1)))) ...
                nanstd(dat(ind2))./sqrt(sum(~isnan(dat(ind2)))) nanstd(dat(ind3))./sqrt(sum(~isnan(dat(ind3))))];
            dat = OpenFRinternal(:,:,icell,2,2);
            dat = dat(:);
            ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,2,1) = [nanmean(dat(ind1)) nanmean(dat(ind2)) nanmean(dat(ind3))];
            ArmBarPlot_LocalNonLocal_MeanSEM(icell,:,2,2) = [nanstd(dat(ind1))./sqrt(sum(~isnan(dat(ind1)))) ...
                nanstd(dat(ind2))./sqrt(sum(~isnan(dat(ind2)))) nanstd(dat(ind3))./sqrt(sum(~isnan(dat(ind3))))];
        end
    
        save(d2(id).name,'ArmBarPlot_LocalNonLocal_MeanSEM','-append')
    else
        load(d2(id).name,'ArmBarPlot_LocalNonLocal_MeanSEM')
    end
    
    
    load(d2(id).name,'other_cells','InFR','OutFR','armposindex','RP_moduarm','RP_pSSDarm','other_cells_touse')
    replay = RP_moduarm';
    FR = InFR(other_cells,:)+OutFR(other_cells,:);
    fr = NaN(size(FR,1),3);
    for iarm = 1:3
        fr(:,iarm) = nanmean(FR(:,armposindex(:,iarm)),2);
    end
    local = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,1,1);
    nonlocal = ArmBarPlot_LocalNonLocal_MeanSEM(other_cells,:,2,1);
    repall = cat(1,repall,replay(other_cells_touse(:,igroup),:));
    othall = cat(1,othall,cat(3,fr(other_cells_touse(:,igroup),:),local(other_cells_touse(:,igroup),:),nonlocal(other_cells_touse(:,igroup),:)));
    sigcells = cat(1,sigcells,RP_pSSDarm(other_cells_touse(:,igroup))<.05);
    
    
    if toplotcells && igroup ==1
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
            helper_saveandclosefig([savefolder '\Figure3\related\2DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(icell) ])                     
        end
       
    end
    disp(['Done with day ' num2str(id)])
end    



othlabel = {'Behavioral FR';'Local';'Non-Local'};
zlab = {'Across Cells';'Median Subtracted';'Zscored';'Norm'};
siglab = {'All Cells';'Sig Cells'};
for iz = 1:4
    if iz == 1
        r = repall;
        o = othall;
    elseif iz == 2
        r = repall-median(repall,2);
        o = othall-median(othall,2);
    elseif iz == 3
        r = zscore(repall,[],2);
        o = zscore(othall,[],2);
    elseif iz == 4
        r = (repall-min(repall,[],2))./range(repall,2);
        o = (othall-min(othall,[],2))./range(othall,2);
    end
    
    for isig = 2
       if isig == 1
           rr = r; oo = o;
       elseif isig == 2
           rr = r(sigcells==1,:); oo = o(sigcells==1,:,:);
       end
       dat1 = rr(:); 
       for ioth = 1:3
           dat2 = oo(:,:,ioth); dat2 = dat2(:);
           figure; hold on
           plot(dat1,dat2,'ok','MarkerSize',10,'LineWidth',3)
           [rho,p] = corr(dat1,dat2,'rows','complete','type','Spearman'); %,'type','Kendall');
           xlabel('Modulation by Replay')
           ylabel(['Modulation by ' othlabel{ioth}])
           set(gca,'FontSize',18)
           xl = get(gca,'xlim'); yl = get(gca,'ylim');
           t1 = text(xl(2)*.4,yl(2)*.7,['r = ' num2str(round(rho,2,'significant')) ' p = ' num2str(round(p,2,'significant'))],'FontSize',18);
           if p<.05
                t1.Color = 'r';
           end
           title([zlab{iz} ' ' siglab{isig} ' - Replay vs ' othlabel{ioth}])
           set(gcf,'Position',[680   431   928   547])
           set(gcf,'renderer','Painters')
           helper_saveandclosefig([savefolder '\Figure3\BarPlots_Scatter_ReplayVS' ...
               othlabel{ioth} ' ' siglab{isig} ' ' zlab{iz}])
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
if toplot
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
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Raw'])
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
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Smoothed'])

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

helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Raw_diff0'])
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
helper_saveandclosefig([savefolder '\Figure3\2DInternal_GroupData_Smoothed_diff0'])
end