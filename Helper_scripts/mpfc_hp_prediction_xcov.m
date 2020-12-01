function xcovlaps = mpfc_hp_prediction_xcov(thisdir,binsize,wind)
%%
load(thisdir,'spikedata','hp_cells','hpinterneurons','other_cells','laps_singlepass','armpos','pos'); %'headingarm','armpos');


hpcells = hp_cells(~ismember(hp_cells,hpinterneurons));
% indinds = indinds-(theta_seqeuence_zero/360)*.125;
% indinds2 = indinds2-(theta_seqeuence_zero/360)*.125;
% s = spikedata(ismember(spikedata(:,2),other_cells),:);
% binsize = .04; 
lapspikes_pfc = cell(max(laps_singlepass),1); lapspikes_hp = lapspikes_pfc;
lapspikesmean_pfc = NaN(max(laps_singlepass),length(other_cells));
lapspikesmean_hp = NaN(max(laps_singlepass),length(hpcells));
laptimes = []; fmto = NaN(max(laps_singlepass),2);
for ilap = 1:max(laps_singlepass)
    
    for icell = 1:length(other_cells)
        slap = spikedata(spikedata(:,2)==other_cells(icell) & laps_singlepass(spikedata(:,3))==ilap,:);        
        bins = pos(find(laps_singlepass==ilap,1,'first'),1):binsize:pos(find(laps_singlepass==ilap,1,'last'),1);
        [h,c] = histcounts(slap(:,1),bins);
        if icell == 1
%             laptimes = cat(1,laptimes,[mean([c(1:end-1); c(2:end)])' ilap*ones(size(h,2),1)]);
            laptimes = cat(1,laptimes,[c' ilap*ones(size(c,2),1)]);
            lapspikes_pfc{ilap,1} = NaN(size(h,2),1);
        end
        lapspikes_pfc{ilap,1}(:,icell) = h;
        lapspikesmean_pfc(ilap,icell) = mean(h);
    end
    for icell = 1:length(hpcells)
        slap = spikedata(spikedata(:,2)==hpcells(icell) & laps_singlepass(spikedata(:,3))==ilap,:);        
        bins = pos(find(laps_singlepass==ilap,1,'first'),1):binsize:pos(find(laps_singlepass==ilap,1,'last'),1);
        h = histcounts(slap(:,1),bins);
        if icell == 1            
            lapspikes_hp{ilap,1} = NaN(size(h,2),1);
        end
        lapspikes_hp{ilap,1}(:,icell) = h;
        lapspikesmean_hp(ilap,icell) = mean(h);
    end
    fmto(ilap,:) = [armpos(find(laps_singlepass==ilap,1,'first')) armpos(find(laps_singlepass==ilap,1,'last'))];
end

%%

xcovlaps = NaN(max(laps_singlepass),round(wind/binsize)*2+1);
for iarm = 1:3    
    laplaps = find(fmto(:,1)==iarm);        
    otharms = setdiff(1:3,iarm);
    for ii = 1:length(laplaps)   
        trainlaps = laplaps(~ismember(laplaps,laplaps(ii)));
        traindat1 = lapspikesmean_pfc(trainlaps,:);
%         traindat = zscore(traindat1); 
        traindat = traindat1;
        trainlabel = fmto(trainlaps,2);
        testdat = lapspikes_pfc{laplaps(ii),1};
        if length(trainlaps)==1 || length(unique(trainlabel))==1
            continue
        end
       excl = std(traindat,[],1)<10^-10 | std(traindat(trainlabel==otharms(1),:),1)<10^-10  | std(traindat(trainlabel==otharms(2),:),1)<10^-10;
       traindat(:,excl) = []; testdat(:,excl) = [];
       [~,~,P] = classify(testdat,traindat,trainlabel,'diaglinear');  
       
        traindat1 = lapspikesmean_hp(trainlaps,:);
%         traindat = zscore(traindat1,[],1);        
        traindat = traindat1;
        testdat = lapspikes_hp{laplaps(ii),1};
        excl = std(traindat,[],1)<10^-10 | std(traindat(trainlabel==otharms(1),:),1)<10^-10  | std(traindat(trainlabel==otharms(2),:),1)<10^-10;
       traindat(:,excl) = []; testdat(:,excl) = [];
       [~,~,H] = classify(testdat,traindat,trainlabel,'diaglinear');  
%         PP = (P(:,1)-P(:,2))./sum(P,2);
%         HH = (H(:,1)-H(:,2))./sum(H,2);
        PP = (P(:,1)-P(:,2));
        HH = (H(:,1)-H(:,2));
        [c1,lags] = xcov(HH,PP,round(wind/binsize),'coeff');
        xcovlaps(laplaps(ii),:) = c1;        
    end
end
figure; plot(-wind:binsize:wind,nanmean(xcovlaps))
