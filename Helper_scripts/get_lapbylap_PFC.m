function get_lapbylap_PFC(thisdir,label,laptype,velcutoff)

%across laps get for each PFC cell:
    %spatial information 
    
    %FR
    %theta locking mrv
    %theta phase preference
    %theta preceission correlation x dir
    %theta preceission slope x dir
        %and across each arm and heading arm
        
    %replay event modulation
    %arm-replay modulation pSSD
    %arm-replay modulation SSD
        %modu across each arm
       
    
% lapdata = Cell Num X lapnum X 10
% laparmdata = Cell Num X lapnum X 8 X 3

laplab = {'laps_coverspace','laps_twoarms','laps_singlepass'};

load(thisdir,'HPTheta_spikephase','spikedata','pos','other_cells','linposcat','armpos', ...
    [label '_PFCcandspikes_list'], ...
[label '_PFCreplayspikes_binned'], [label  '_replay_stnd'],...
[label  '_replay_replayarm'],[label  '_replay_jointreplayarm'], ...
[label  '_replay_singlebothjoint'],[label '_CandEventAssign'],[label '_CandEventTimes'],...
    [label  '_replay_jointjumptime'],[label  '_replay_shuffle_p'],'laps_coverspace','laps_singlepass','laps_twoarms','binpos','vel','dirdat')

eval(['replayarm = ' label '_replay_replayarm;'])
eval(['CTimes = ' label '_CandEventTimes(' label '_CandEventAssign==1);'])
eval(['CandList1 = ' label '_PFCcandspikes_list;'])
eval(['ReplayBinned1 = ' label '_PFCreplayspikes_binned;'])
eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
clear([label  '_replay_replayarm'],[label '_CandEventAssign'],[label '_CandEventTimes'],[label  '_replay_jointreplayarm'],[label '_PFCcandspikes_list'], [label '_PFCreplayspikes_binned'])


eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
% ReplayBinned1(:,:,singlebothjoint==3 | shuffle_p>.05) = [];

eval(['CTimes2 = ' label '_replay_stnd(:,1);'])
CTimes2(singlebothjoint==3 | shuffle_p>.05,:) = [];

eval(['laps = ' laplab{laptype} ';'])
clear('laps_coverspace','laps_singlepass','laps_twoarms')
[~,~,bin] = histcounts(CTimes,pos(:,1));
iind = find(bin==0); bin(iind(iind>(size(bin,1)/2))) = bin(find(bin~=0,1,'last')); bin(iind(iind<(size(bin,1)/2))) = bin(1);
CandEvLap = laps(bin);

[~,~,bin2] = histcounts(CTimes2,pos(:,1));
iind = find(bin2==0); bin2(iind(iind>(size(bin2,1)/2))) = bin2(find(bin2~=0,1,'last')); bin2(iind(iind<(size(bin2,1)/2))) = bin2(1);
ReplayLap = laps(bin2);

binsize = 10;
modind = [0 .3];
baseind = [-.5 -.1];

lapdata = NaN(length(other_cells),max(laps),12);
laparmdata = NaN(length(other_cells),max(laps),8,3);
s = spikedata(ismember(spikedata(:,2),other_cells),:);
% cellfr = hist(s(vel(s(:,3))>velcutoff,2),other_cells)./(sum(vel>velcutoff)/30);
Th = HPTheta_spikephase(ismember(HPTheta_spikephase(:,2),other_cells),:);
for ilap = 1:max(laps)
    CandList = CandList1(ismember(CandList1(:,3),find(CandEvLap==ilap)),:);
    ReplayBinned = ReplayBinned1(:,:,ReplayLap==ilap);
    
    %Firing Rate
    lapdata(:,ilap,1) = hist(s(laps(s(:,3))==ilap & vel(s(:,3))>velcutoff,2),other_cells)./(sum(laps==ilap & vel>velcutoff)/30);
    for iarm = 1:3
        laparmdata(:,ilap,1,iarm) = hist(s(laps(s(:,3))==ilap & armpos(s(:,3))==iarm & vel(s(:,3))>velcutoff,2),other_cells)./(sum(armpos==iarm &laps==ilap & vel>velcutoff)/30);
    end
    
    %Spatial Information
    ppos = hist(binpos(laps==ilap & vel>velcutoff),min(binpos):max(binpos))./(sum(laps==ilap & vel>velcutoff));        
    for icell = 1:length(other_cells)
        pfrC = hist(binpos(s(laps(s(:,3))==ilap & vel(s(:,3))>velcutoff & s(:,2)==other_cells(icell),3)),min(binpos):max(binpos))./...
            (hist(binpos(laps==ilap & vel>velcutoff),min(binpos):max(binpos))./30);        
        lapdata(icell,ilap,2) = (nansum([ppos'.*pfrC']'.*log2(pfrC./nanmean(pfrC))))./nanmean(pfrC);
    end
    clear ppos
    
    %Theta
    for icell = 1:length(other_cells)
        %Theta Locking
        
        ph = Th(Th(:,2)==other_cells(icell) & laps(Th(:,3))==ilap & vel(Th(:,3))>velcutoff,6);
        if length(ph)>2
        lapdata(icell,ilap,3) = circ_r(deg2rad(ph));
        %Theta phase pref
        lapdata(icell,ilap,4) = abl_circ_mean_ci(deg2rad(ph),[],[],[],[],1);
        end
        
        %Theta locking and phase pref by arm
        for iarm = 1:3
            ph = Th(Th(:,2)==other_cells(icell) & armpos(Th(:,3))==iarm & laps(Th(:,3))==ilap & vel(Th(:,3))>velcutoff,6);
            laparmdata(icell,ilap,2,iarm) = circ_r(deg2rad(ph));     
            if ~isempty(ph)
                laparmdata(icell,ilap,3,iarm) = abl_circ_mean_ci(deg2rad(ph),[],[],[],[],1);
            end
        end
        
        %Theta precession and slope
        for inout = 1:2
            ph = Th(Th(:,2)==other_cells(icell) & laps(Th(:,3))==ilap & dirdat(Th(:,3))==(inout-1)& vel(Th(:,3))>velcutoff,6);
            phpos = linposcat(Th(Th(:,2)==other_cells(icell) & dirdat(Th(:,3))==(inout-1) & laps(Th(:,3))==ilap & vel(Th(:,3))>velcutoff,3));
            if length(ph)>2
                [r,sl,~,~] = get_pp_slope(phpos,deg2rad(ph),deg2rad(binsize));
                lapdata(icell,ilap,5+inout-1) = r;
                lapdata(icell,ilap,7+inout-1) = sl;
            end
            %by arm
            for iarm = 1:3
                ph = Th(Th(:,2)==other_cells(icell) & laps(Th(:,3))==ilap & armpos(Th(:,3))==iarm & dirdat(Th(:,3))==(inout-1)& vel(Th(:,3))>velcutoff,6);
                 phpos = linposcat(Th(Th(:,2)==other_cells(icell) & armpos(Th(:,3))==iarm & dirdat(Th(:,3))==(inout-1) & laps(Th(:,3))==ilap & vel(Th(:,3))>velcutoff,3));
                 if length(ph)>2
                    [r,sl,~,~] = get_pp_slope(phpos,deg2rad(ph),deg2rad(binsize));
                    laparmdata(icell,ilap,4+inout-1,iarm) = r;
                    laparmdata(icell,ilap,6+inout-1,iarm) = sl;
                 end
            end
        end                             
    end
    
    %ripple/sd event or replay of arm triggered
    for icell = 1:length(other_cells)
        [moduE,~] = forcell_get_sig_modu_event(other_cells(icell),modind,baseind,CandList);
        [moduA,p_SSD,armSSD1,armSSD2] = forcell_get_sig_modu_byarm(icell,modind,baseind,replayarm(ReplayLap==ilap),ReplayBinned);
        lapdata(icell,ilap,9) =moduE;
        lapdata(icell,ilap,10) =p_SSD;
        lapdata(icell,ilap,11) =armSSD1;
        lapdata(icell,ilap,12) =armSSD2;
        laparmdata(icell,ilap,8,:) = moduA;
    end
%     disp(['Done with lap ' num2str(ilap) ' in get_lapbylap_PFC ' label])    
end

eval([label 'laptype' num2str(laptype) '_lapdata = lapdata;'])
eval([label 'laptype' num2str(laptype) '_laparmdata = laparmdata;'])
save(thisdir,[label 'laptype' num2str(laptype) '_lapdata'],[label 'laptype' num2str(laptype) '_laparmdata'],'-append')

disp(['Done with get_lapbylap_PFC ' label])