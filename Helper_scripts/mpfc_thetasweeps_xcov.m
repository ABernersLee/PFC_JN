function [xcovlaps,lags] = mpfc_thetasweeps_xcov(thisdir,binsize,wind,cutoff)
load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','spikedata','other_cells','laps_singlepass','armpos','pos'); %'headingarm','armpos');
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
clear times_armon_thetaof_headingarm_lap_thetahalf_all

b = [];
for ilap = 1:max(laps_singlepass)        
        bins = pos(find(laps_singlepass==ilap,1,'first'),1):binsize:pos(find(laps_singlepass==ilap,1,'last'),1);
        b = [b;bins(1:end-1)'];
end
lapspikes_pfc = NaN(length(b),length(other_cells),2); lapspikes_hp = lapspikes_pfc;
laptimes = []; fmto = NaN(max(laps_singlepass),2);
for icell = 1:length(other_cells)
    lp = [];
    for ilap = 1:max(laps_singlepass)
        slap = spikedata(spikedata(:,2)==other_cells(icell) & laps_singlepass(spikedata(:,3))==ilap,:);        
        bins = pos(find(laps_singlepass==ilap,1,'first'),1):binsize:pos(find(laps_singlepass==ilap,1,'last'),1);
        [h,c] = histcounts(slap(:,1),bins);
        if icell == 1            
           laptimes = cat(1,laptimes,[c' ilap*ones(size(c,2),1)]);                        
           fmto(ilap,:) = [armpos(find(laps_singlepass==ilap,1,'first')) armpos(find(laps_singlepass==ilap,1,'last'))];
        end
        lp = cat(1,lp,[h' ilap*ones(size(h,2),1)]);                
    end
    lapspikes_pfc(:,icell,:) = lp;    
end
lapspikes_pfc2 = lapspikes_pfc;
lapspikes_pfc2(:,1) = zscore(lapspikes_pfc2(:,1));
lapspikesmean_pfc =NaN(max(laps_singlepass),length(other_cells));
for ilap = 1:max(laps_singlepass)
    for icell = 1:size(lapspikes_pfc2,2)
        lapspikesmean_pfc(ilap,icell) = nanmean(lapspikes_pfc2(lapspikes_pfc2(:,icell,2)==ilap,icell,1));
    end
end

%%
% wind = .25;
Th(Th(:,9)<cutoff,:) = [];
xcovlaps = NaN(max(laps_singlepass),round(wind/binsize)*2+1,3,3);
for iarm = 1:3    
    laplaps = find(fmto(:,1)==iarm);        
    
    for ii = 1:length(laplaps)   
        trainlaps = laplaps(~ismember(laplaps,laplaps(ii)));
        traindat = lapspikesmean_pfc(trainlaps,:);        
        trainlabel = fmto(trainlaps,2);
        testdat = lapspikes_pfc2(lapspikes_pfc2(:,1,2)==laplaps(ii),:,1);
        if length(trainlaps)==1 || length(unique(trainlabel))==1
            continue
        end
%        traindat = dat(:,train);       
       excl = std(traindat,[],1)<10^-10;
       traindat(:,excl) = []; testdat(:,excl) = [];
%         guess = classify(testdat,traindat,trainlabel,'diaglinear');  
        [~,~,P] = classify(testdat,traindat,trainlabel,'diaglinear');  
        Thlap = Th(Th(:,6)==laplaps(ii),[1 2 4]);
        Thlap(:,4) = Thlap(:,3);
        Thlap(:,3) = mean([Thlap(:,1) Thlap(:,2)],2);
        harm = fmto(laplaps(ii),2);
        harmind = find(setdiff(1:3,iarm)==harm);        
        for stnd = 1:3
            %theta sweep and mpfc both going towards heading direction
            thetasweep1 = Thlap(Thlap(:,4)==harm,1:3);
            thislaptimes = laptimes(laptimes(:,2)==laplaps(ii),1);
            h = histcounts(thetasweep1(:,stnd),thislaptimes);
%             [c1,lags] = xcov(h,guess==harm,round(wind/binsize),'coeff');
            [c1,lags] = xcov(h,P(:,harmind),round(wind/binsize),'coeff');

            %theta sweep and mpfc both going away from heading direction
            thetasweep1 = Thlap(Thlap(:,4)~=harm,1:3);
            thislaptimes = laptimes(laptimes(:,2)==laplaps(ii),1);
            h = histcounts(thetasweep1(:,stnd),thislaptimes);
%             cc = xcov(h,guess~=harm,round(wind/binsize),'coeff');
            cc = xcov(h,P(:,setdiff(1:3,iarm)~=harm),round(wind/binsize),'coeff');
            c2 = nanmean([c1;cc]);
%             c2 = mean([c1;cc]);
%             c2 = cc;
            
            %mpfc is going in the heading direction
            thetasweep1 = Thlap(:,1:3);
            thislaptimes = laptimes(laptimes(:,2)==laplaps(ii),1);
            h = histcounts(thetasweep1(:,stnd),thislaptimes);
%             c3 = xcov(h,guess==harm,round(wind/binsize),'coeff');
            c3 = xcov(h,P(:,harmind),round(wind/binsize),'coeff');
            
            xcovlaps(laplaps(ii),:,:,stnd) = [c1;c2;c3]';
        end

    end
end
