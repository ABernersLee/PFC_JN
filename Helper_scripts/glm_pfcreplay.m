% discrim_pfcfwdreplays
function glm_pfcreplay(dirs)

label = 'RP';


cd(dirs.homedir)
d2 = dir('*.mat');

ts = [];
ps = [];
for id = 1:size(d2,1)
        
    shuff_acc = [];
    real_acc = [];
    thisdir = d2(id).name;
    load(thisdir,[label '_replay_singlebothjoint'],[label '_Cand_sig_modu_include'], ...
        [label '_PFCreplayspikes_binned'],[label  '_replay_shuffle_p']...
        ,[label  '_replay_replayarm'],[label '_replay_stnd'],...
        [label  '_replay_corr'],[label  '_replay_dirbias'])
%     load(thisdir,[label '_HPreplayspikes_binned'])
  load(d2(id).name, [label '_moduarm'],[label '_pSSDarm'],[label '_SSDarm'],[label '_Cand_sig_modu_include'])
    eval(['armsig = [' label '_pSSDarm ' label '_SSDarm ' label '_Cand_sig_modu_include];'])       

    eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
%     eval(['HPreplayspikes_binned = ' label '_HPreplayspikes_binned;'])
%     eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
    eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
    eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
    eval(['replayarm = ' label '_replay_replayarm;'])
    eval(['Event = ' label '_replay_stnd(:,1);'])
    eval(['dirbias = ' label '_replay_dirbias;'])
    eval(['corrs = ' label '_replay_corr;'])

    eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
%     clear([label '_PFCreplayspikes_list'])
    clear([label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])
%     load(thisdir,'RP_CandEventTimes','HP_Ripple','RP_CandEventAssign')
%     RP_CandEventAssign(end,:) = [];
    load(thisdir,'HP_Ripple')
    
       
    % 1 is pure inbound, -1 is pure outbound
    % dirbias = dirbias(touse);

    % outbound is positive, inbound is negative
    % corrs = corrs(touse);

%     fwd = (dirbias>0 & corrs<0) | (dirbias<0 & corrs>0);
%     fwd = (dirbias<0 & corrs>0);
    % rev = (dirbias>0 & corrs>0) | (dirbias<0 & corrs<0);
    
    PFCreplayspikes_binned = PFCreplayspikes_binned(Cand_sig_modu_include(:,1)==1,:,:);
    
    %down-sampling the HP cells to match the PFC cell count
%     numcells = 1:size(HPreplayspikes_binned,1); numcells = numcells(randperm(length(numcells)));
%     HPreplayspikes_binned = HPreplayspikes_binned(numcells(1:size(PFCreplayspikes_binned,1)),:,:);
    

    
    
    load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos','linposcat')
    [~,armp] = max(armposindex,[],2); clear armposindex
    lappro = NaN(size(laps_singlepass,1),3);
    lapout = NaN(size(laps_singlepass,1),2);
    % lapdir = NaN(size(laps_singlepass,1),1);
    lapi = NaN(max(laps_singlepass),2);
    for ilap = 1:max(laps_singlepass)
       leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
       lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle      
       lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform

       nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
       thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
       lapi(ilap,:) = [thisarm nextarm];
    end

    lapi(:,3) = [NaN;lapi(1:end-1,1)];
    Event(singlebothjoint==3 | shuffle_p>=.05) = [];
    replayarm(singlebothjoint==3 | shuffle_p>=.05) = [];
    dirbias(singlebothjoint==3 | shuffle_p>=.05) = [];
    corrs(singlebothjoint==3 | shuffle_p>=.05) = [];
    
    [~,~,i] = histcounts(Event,HP_Ripple(:,1));
    RipPow = HP_Ripple(i,2);
    
    [~,~,i] = histcounts(Event,pos(:,1));
    % Event(i==0) = [];
    replayarm(i==0) = [];
    armon = armpos(i(i~=0));
    lapon = laps_singlepass(i(i~=0));
%     ulapon = unique(lapon);
    
    headingarm = NaN(size(lapon,1),1);
    armon2 = headingarm; lastarm = armon2;
    for i = find(~isnan(lapon))
        headingarm(i) = lapi(lapon(i),2);
        armon2(i) = lapi(lapon(i),1);
        lastarm(i) = lapi(lapon(i),3);
    end
%     armon2 = lapi(lapon,1);
    touse = true(size(headingarm,1),1);
%     touse(isnan(lapon)) = false;
    touse(armon~=armon2 | isnan(lapon) | lapon==1) = false;   
    
    
%     Event = Event(touse,1);
    %These are the things that will go into the LM
    replayarm = replayarm(touse);
    armon = armon(touse);
    headingarm = headingarm(touse);
    dirbias = dirbias(touse);
    corrs = corrs(touse);
    lastarm = lastarm(touse);
    RipPow = RipPow(touse);
    
    PFCreplayspikes_binned = PFCreplayspikes_binned(:,:,touse);
%     HPreplayspikes_binned = HPreplayspikes_binned(:,:,touse);
    
    
    modind1 = [0 .2];
    baseind1 = [-.5 -.1];
    binsize = .02;
    window = [-2 2];
    ind = [window(1):binsize:window(2)];
    modind = ind>=modind1(1) & ind<=(modind1(2));
    baseind = ind>=baseind1(1) & ind<=baseind1(2);

         
    pfcfr2 = squeeze((nanmean(PFCreplayspikes_binned(:,modind,:),2)./range(modind1))-(nanmean(PFCreplayspikes_binned(:,baseind,:),2)./range(baseind1)))';
    pfcfr = zscore(pfcfr2);
%     pfcfr = zscore(squeeze((nanmean(HPreplayspikes_binned(:,modind,:),2)./range(modind1))-(nanmean(HPreplayspikes_binned(:,baseind,:),2)./range(baseind1)))');        
   
%     X = [replayarm armon lastarm headingarm]; % with 2, (10)
    X = [replayarm armon]; % lastarm headingarm dirbias>0 corrs>0]; % with 1
    D = x2fx(X,'interaction');
%     X = [replayarm armon dirbias>0 corrs>0]; % with 1
    
%     pday = NaN(size(pfcfr,2),size(X,2));
    pday = NaN(size(pfcfr,2),size(D,2)-1);
    
    tbl2 = pday;
    for icell = 1:size(pfcfr,2)
    %     [b,dev,stats] = glmfit(X,pfcfr(:,icell));
%         [pday(icell,:),tbl,stats] = anovan(pfcfr(:,icell),X,'model',1,'display','off'); %,'model',[1 0 0 0 0;0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 1 1 0 0 0]);
%         tbl2(icell,:)  = ridge(pfcfr2(:,icell),X,0);
        tbl2(icell,:)  = ridge(pfcfr2(:,icell),D(:,2:end),0);
%         for i = 2:7
%             tbl2(icell,i-1) = tbl{i,6};
%         end
    end
    ts = [ts;tbl2];
    ps = [ps;pday];
end
% ts(isnan(ps)) = NaN;


nanmean(abs(ts))
sum(ps<.05)./sum(~isnan(ps))