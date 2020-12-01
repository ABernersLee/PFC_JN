%Figure3_Internal1d

% newcol = [75 0 130;255 130 0;34 139 34]/255;

d2 = dir('*.mat');
binsize = 20;
SigmaField = 1;
One_D_Filter=fspecial('gaussian',[3 1],SigmaField);
velcutoff = 5;
EstBin = .1; 
cellcutoff = 3;
spikecutoff = 5;
vellab = {'Low Velocity';'High Velocity'};
loclab = {'Local';'Non-Local'};
newcol = [75 0 130;34 139 34;255 130 0]/255;
toplot = true;
% tocalc = true;
for id = 1:size(d2,1)
    clear OpenFRinternal
    load(d2(id).name,'spikedata','hp_cells','other_cells','hpinterneurons','pos','FRinternal','RP_pSSDarm')    
                      thisdir = d2(id).name;
%     if exist('FRinternal','var'); tocalc = false; else; tocalc = true; end
    tocalc = true;
    if tocalc % get the internal fields        
        
        if 1 % get the error indicies
        
            load(d2(id).name,'RP_pSSDarm','RP_SSDarm','InFR','OutFR','hp_cells','hpinterneurons','spikedata','pos','binpos','armpos','armposindex','vel','other_cells')


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

            

            FR(sum(FR~=0,2)==0,:) = 1;
            FR(FR==0 | FR<0) = .0001;

            touse = sum(binspike>0)>=cellcutoff; % has to have at least 5 cells spiking
            touse = touse & sum(binspike)>=spikecutoff;
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
            % %make times when fewer than 5 cells are spiking equal to zeros/NaNs all the way down
            

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
%                 [~,I2] = max(Mat);
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
%             err(diffarm) = (dp(diffarm)-toadd(da(diffarm)))+(bp(diffarm)-toadd(ba(diffarm)));             
            err = diffarm;
            [~,dp2] = histc(dp,0:size(armposindex,1));
        end
    
        Cell_Number = 1:size(InFR,1);
            TimeBins=round((End_Time-Start_Time)/EstBin)*4;
         binspike2=zeros(length(Cell_Number),TimeBins);
%          bintimes2 = NaN(TimeBins,1);
            for i=1:4
                for CellID=1:length(Cell_Number)   
                    c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
                    binspike2(CellID,i:4:TimeBins)=c(1:TimeBins/4);
%                     jnk = Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins;
%                     if length(jnk)>length(i:4:TimeBins)
%                         bintimes2(i:4:TimeBins) = jnk(1:end-1)+(EstBin/2);
%                     else
%                         bintimes2(i:4:TimeBins) = jnk+(EstBin/2);
%                     end
            %         binspike(CellID,i:4:TimeBins)=i*10*ones(length(1:TimeBins/4),1);         %testing
                end                      
            end
            bsp3 = binspike2(:,touse);

        FRinternal = NaN(max(dp2(:,1)),size(InFR,1),2,2);
        for ivel = 1:2
            if ivel==1; velind = bv<velcutoff; elseif ivel==2; velind = bv>velcutoff; end
%             for ilocal = 1:2
%                 if ilocal ==1; locind = err<errcutoff; elseif ilocal ==2; locind = err>errcutoff; end
                pos4 = dp2;
%                 pos4(~velind | ~locind,:) = NaN; %restrict to times with the correct velocity and locality
                pos4(~velind,:) = NaN; %restrict to times with the correct velocity and locality                
                pos3 = pos4; %sub2ind(max(dp2),pos4(:,1));
                               
                
                ind = 1:max(dp2(:,1));
                occ = NaN(size(ind,2),2);
                f = zeros(size(ind,2),size(InFR,1),2);
                for j = 1:length(ind)
                    m = median(err(pos3==ind(j)));
%                     m = 2;
                    if m==0
                        m = .0001;
                    end           
                    m= 0;
                    
                    if ~isnan(m)
                        mind1 = err<=m;
                        mind2 = err>m;
                        
                        %downsample to make the same number in each bin
                        num1 = sum(pos3==ind(j) & mind1);
                        num2 = sum(pos3==ind(j) & mind2);
                        if (num2-num1)<0
                            jnk = find(pos3==ind(j) & mind1);
                            mind1(jnk(randperm(num1-num2))) = false;
                        elseif (num2-num1)>0
                            jnk = find(pos3==ind(j) & mind2);
                            mind2(jnk(randperm(num2-num1))) = false;                            
                        end
                        
                        num1 = sum(pos3==ind(j) & mind1);
                        num2 = sum(pos3==ind(j) & mind2);
                        if (num1-num2)~=0
                            disp(num2str(num1-num2))
                        end
                        
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
                   
                figure; hold on
                subplot(3,1,1)
                plot(occ(:,1),'k'); 
                title([vellab{ivel} ' ' loclab{1}])
                subplot(3,1,2)
                plot(occ(:,2),'k'); 
                title([vellab{ivel} ' ' loclab{2}])
                subplot(3,1,3)
                plot(occ(:,2)-occ(:,1),'k.');
                title('Difference in Occupancy')
                set(gcf,'Position',[2562           6         494         879])
%                 helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Figure3\related\1DInternal_' vellab{ivel} '_' d2(id).name(1:end-4) '_Occ']) 
                        close gcf
                
                FRinternal_temp1 = f(:,:,1)./repmat(occ(:,1),[1 size(f,2)]);
                FRinternal_temp2 = f(:,:,2)./repmat(occ(:,2),[1 size(f,2)]);
                
                FRinternal(:,:,ivel,:) = cat(4,FRinternal_temp1,FRinternal_temp2);
                
%             end
        end
%         save(d2(id).name,'FRinternal','-append')
    else
        load(d2(id).name,'FRinternal')
    end

    other_cells_sig = other_cells(RP_pSSDarm<.05);
    
    
   
    fr2 = FRinternal;
    % Plot cells
    if toplot
    for icell = 1:length(other_cells_sig)
        figure; hold on
        subplot(1,3,1); hold on
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),1,2);
            errorbar(iarm,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'Color',newcol(iarm,:),'LineWidth',3)
        end
        xlim([.5 3.5])
        subplot(1,3,2); hold on
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),2,1);
            errorbar(iarm,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'Color',newcol(iarm,:),'LineWidth',3)
        end
        xlim([.5 3.5])
        subplot(1,3,3); hold on
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),2,2);
            errorbar(iarm,nanmean(dat),nanstd(dat)./sqrt(length(dat)),'Color',newcol(iarm,:),'LineWidth',3)
        end
        xlim([.5 3.5])
        suptitle([thisdir(1) thisdir(3:end-4) ' PFC Cell ' num2str(other_cells_sig(icell))])
        helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\1DInternal_Errorbars_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(other_cells_sig(icell))]) 
        
        figure; hold on
        subplot(2,1,1); hold on
        yyaxis left
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),1,1);
%             dat = filter(One_D_Filter,1,dat);
            plot(find(armposindex(:,iarm)),dat,'--','color',newcol(iarm,:),'LineWidth',3)
        end
        ylabel('Local (Dashed)','Color','k')
        set(gca,'YColor','k')

        yyaxis right
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),1,2);
%             dat = filter(One_D_Filter,1,dat);
            plot(find(armposindex(:,iarm)),dat,'-','color',newcol(iarm,:),'LineWidth',3)
        end
        ylabel('Non-Local (Solid)','Color','k')
        set(gca,'YColor','k')
        xlabel('cm along track (Low Velocity)')        
        xlim([1 size(armposindex,1)])
        set(gca,'FontSize',18)
        title([thisdir(1) thisdir(3:end-4) ' PFC Cell ' num2str(other_cells_sig(icell))])

        subplot(2,1,2); hold on
        yyaxis left
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),2,1);
%             dat = filter(One_D_Filter,1,dat);
            a = plot(find(armposindex(:,iarm)),dat,'--','color',newcol(iarm,:),'LineWidth',3);
        end
        ylabel('Local (Dashed)','Color','k')
        set(gca,'YColor','k')   
        yyaxis right
        for iarm = 1:3
            dat  = fr2(armposindex(:,iarm),other_cells_sig(icell),2,2);
%             dat = filter(One_D_Filter,1,dat);
            b = plot(find(armposindex(:,iarm)),dat,'-','color',newcol(iarm,:),'LineWidth',3);
        end
        ylabel('Non-Local (Solid)','Color','k')
        set(gca,'YColor','k')    
        xlabel('cm along track (High Velocity)')
        xlim([1 size(armposindex,1)])
        set(gca,'FontSize',18)
        set(gcf,'Position',[2050          33        1158         865])
        legend([a,b],'Local','Non-Local')

        set(gcf,'renderer','Painters')
        helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\InternalFields\1DInternal_' cellid '_' d2(id).name(1:end-4) '_Cell' num2str(other_cells_sig(icell))]) 
    end
    end
    dat1 = fr2(:,other_cells_sig(icell),1,2);
dat2 = fr2(:,other_cells_sig(icell),2,1);
dat3 = fr2(:,other_cells_sig(icell),2,2);
[r,p] = corr(dat1,dat2,'rows','complete');
[r,p] = corr(dat1,dat3,'rows','complete');
[r,p] = corr((dat1-min(dat1))./(max(dat1)-min(dat1)),(dat2-min(dat2))./(max(dat2)-min(dat2)),'rows','complete');
[r,p] = corr((dat1-min(dat1))./(max(dat1)-min(dat1)),(dat3-min(dat3))./(max(dat3)-min(dat3)),'rows','complete');

end
    
%add quantification 
%and shuffle?



% below hasn't been edited

%% group data

  RS = []; SSD = []; SIG = []; FR = []; FRs = [];
for id = 1:size(d2,1)
    clear OpenFRinternal
%     load(d2(id).name,'spikedata','RP_SSDarm','hp_cells','other_cells','hpinterneurons','pos','OpenFRinternal')  
    load(d2(id).name,'RP_pSSDarm','RP_SSDarm','FRinternal','other_cells')

    
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
        rs(icell,1,1) = corr(dat1(:),dat2(:),'rows','complete');        
        rs(icell,2,1) = corr(dat1(:),dat3(:),'rows','complete');        
        
        %corr to low velocity non-local - smoothed
        dat1 = smdat(:,:,1,2);
        dat2 = smdat(:,:,2,1); %high velocity local
        dat3 = smdat(:,:,2,2); %high velocity non-local
        rs(icell,1,2) = corr(dat1(:),dat2(:),'rows','complete');        
        rs(icell,2,2) = corr(dat1(:),dat3(:),'rows','complete');        
        
        
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
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Figure3\related\2DInternal_GroupData_Raw'])
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
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Figure3\related\2DInternal_GroupData_Smoothed'])

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

helper_savefig(['E:\XY_matdata\Figures\ForPaper\Figure3\related\2DInternal_GroupData_Raw_diff0'])
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
helper_savefig(['E:\XY_matdata\Figures\ForPaper\Figure3\related\2DInternal_GroupData_Smoothed_diff0'])
