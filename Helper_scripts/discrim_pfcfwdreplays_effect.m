function discrim_pfcfwdreplays_effect(dirs,numshuff) %,ispfc)
% discrim_pfcfwdreplays % excact script that produced the fact that PFC and
% not HP can decode the animal's choice from the response to replays
%had exactly matched number of cells
%have changed since then

label = 'RP';


cd(dirs.homedir)
d2 = dir('*.mat');

% realall = []; shuffall = [];
% numshuff= 500;
realacc = [];
shuffacc = [];
realacc2 = [];
shuffacc2 = [];
daydat = [];
pid = NaN(size(d2,1),1);
numdaycells = NaN(size(d2,1),1);
for id = 1:size(d2,1)
        
    shuff_acc = [];
    real_acc = [];
    thisdir = d2(id).name;
    load(thisdir,[label '_replay_singlebothjoint'],[label '_Cand_sig_modu_include'], ...
        [label '_PFCreplayspikes_binned'],[label  '_replay_shuffle_p']...
        ,[label  '_replay_replayarm'],[label '_replay_stnd'],...
        [label  '_replay_corr'],[label  '_replay_dirbias'])
    
    
%     if ispfc
        eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
%         savelab = 'allPFCcells';
%     else
%         savelab = 'allHPcells';
        load(thisdir,[label '_HPreplayspikes_binned'])
        eval(['HPreplayspikes_binned = ' label '_HPreplayspikes_binned;'])
%     end
    savelab = 'downsampledHPcells';
    
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

    % 1 is pure inbound, -1 is pure outbound
    % dirbias = dirbias(touse);

    % outbound is positive, inbound is negative
    % corrs = corrs(touse);

%     fwd = (dirbias>0 & corrs<0) | (dirbias<0 & corrs>0);
%     fwd = (dirbias<0 & corrs>0);
%     fwd = corrs>0;
    % rev = (dirbias>0 & corrs>0) | (dirbias<0 & corrs<0);
    fwd = true(length(dirbias),1);


    fwd(singlebothjoint==3 | shuffle_p>=.05) = [];
    Event(singlebothjoint==3 | shuffle_p>=.05) = [];
    replayarm(singlebothjoint==3 | shuffle_p>=.05) = [];

    touse = fwd;
    
    PFCreplayspikes_binned = PFCreplayspikes_binned(:,:,touse);
    
%     HPreplayspikes_binned = HPreplayspikes_binned(:,:,touse);
    
%     PFCreplayspikes_binned =     PFCreplayspikes_binned(Cand_sig_modu_include(:,1)==1,:,:); %commented     1/9/18
    
    %down-sampling the HP cells to match the PFC cell count
    numcells = 1:size(HPreplayspikes_binned,1); numcells = numcells(randperm(length(numcells)));
    HPreplayspikes_binned = HPreplayspikes_binned(numcells(1:size(PFCreplayspikes_binned,1)),:,:);
    
    Event = Event(touse,1);
    replayarm = replayarm(touse);
    
    
    load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos')
%     [~,armp] = max(armposindex,[],2); clear armposindex
%     lappro = NaN(size(laps_singlepass,1),3);
%     lapout = NaN(size(laps_singlepass,1),2);
%     lapdir = NaN(size(laps_singlepass,1),1);
    lapi = NaN(max(laps_singlepass),2);
    for ilap = 1:max(laps_singlepass)
%        leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
%        lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle      
%        lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform

       nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
       thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
       lapi(ilap,:) = [thisarm nextarm];
    end
    lapi(:,3) = [NaN;lapi(1:end-1,1)];
    
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
    
    modind1 = [0 .2];
    baseind1 = [-.5 -.1];
    binsize = .02;
    window = [-2 2];
    ind = [window(1):binsize:window(2)];
    modind = ind>=modind1(1) & ind<=(modind1(2));
    baseind = ind>=baseind1(1) & ind<=baseind1(2);

%     pfcfr = zscore(squeeze((nanmean(PFCreplayspikes_binned(:,modind,:),2)./range(modind1))-(nanmean(PFCreplayspikes_binned(:,baseind,:),2)./range(baseind1)))');        
    pfcfr = zscore(squeeze((nanmean(HPreplayspikes_binned(:,modind,:),2)./range(modind1))-(nanmean(HPreplayspikes_binned(:,baseind,:),2)./range(baseind1)))');        
    for iarm = 1:3     
        for iarm2 = 1:3
        repind = armon2==iarm & replayarm==iarm2 & armon==armon2; %ssetdiff(1:3,[headingarm armon]);
%         repind = armon2==iarm & replayarm==headingarm & armon==armon2; %ssetdiff(1:3,[headingarm armon]);
        
%         repind(replayarm==armon | isnan(lapon) | lapon==1) = false;    
        repind(isnan(lapon) | lapon==1) = false;    
        pfc = pfcfr(repind,:);
%         pfc = replayarm(repind);
%         Labels = lapi(lapon(repind),2);
        Labels = headingarm(repind);
%         Labels = lastarm(repind);
%         pfc = pfc(Labels~=iarm,:);
%         Labels(Labels==iarm) = [];
        if all(histc(Labels,unique(Labels))>2) && length(unique(Labels))>1
%             folds = 10;
            folds = length(Labels);
            jnk = ones(ceil(length(Labels)/folds),1)*[1:folds];
            inx = jnk(:); inx = inx(randperm(length(inx))); inx = inx(1:length(Labels));
            
            
            %OG has this commented
            nums = NaN(folds,1);
            for ifold = 1:folds
                trainlabels = Labels(inx~=ifold);
                nums(ifold) = length(unique(trainlabels));
            end
            
            if sum(nums==2)~=length(nums)
                disp(['Not enough labels on arm ' num2str(iarm) ' arm2 ' num2str(arm2)])
                continue
            end
            
            for ifold = 1:folds
                train = inx~=ifold; test = inx==ifold;
                traindata = pfc(train,:); testdata = pfc(test,:);
                trainlabels = Labels(train);
                Mdl = fitcdiscr(traindata,trainlabels,'discrimType','pseudoLinear');           
                
                ns = histc(trainlabels,unique(trainlabels));
                if sum(ns)<171                    
                   numu =  factorial(sum(ns))./(factorial(ns(1))*factorial(ns(2)));
                else numu = numshuff+1;
                    
                end
                
                if numu<numshuff
                    disp([numu ifold iarm iarm2])
                    continue                
                end
                
                guess = predict(Mdl,testdata);
                corrReal = guess==Labels(test);
                real_acc = [real_acc; corrReal]; 
                corrShuff = NaN(length(corrReal),numshuff);
                shufflabels = NaN(length(trainlabels),numshuff);
                
                for ishuff = 1:numshuff
                    ind = randperm(length(trainlabels));
%                     disp('start perm')
                    %OG has this commented, second part of this is new
                    while sum(trainlabels(ind)~=trainlabels)==0 || size(unique([shufflabels(:,1:ishuff-1) trainlabels(ind)]','rows'),1)<ishuff
                        ind = randperm(length(trainlabels));
                    end
%                     disp('end perm')
                    shufflabels(:,ishuff) = trainlabels(ind);           
                    
                    Mdl = fitcdiscr(traindata,shufflabels(:,ishuff),'discrimType','pseudoLinear');                        
                    guess = predict(Mdl,testdata);
                    corrShuff(:,ishuff) = guess==Labels(test);                                
                end
                shuff_acc = [shuff_acc; corrShuff];

            end
        else
%             disp(['skip arm ' num2str(iarm) ' day ' num2str(id)])
%             disp(['skip arm ' num2str(iarm) ' arm ' num2str(iarm2) ' day ' num2str(id)])
        end
            disp(['done arm ' num2str(iarm) ' arm ' num2str(iarm2) ' day ' num2str(id)])
        end
    end
    if ~isempty(real_acc)
    pid(id) = (sum(nanmean(shuff_acc)>=nanmean(real_acc))+1)/(size(shuff_acc,2)+1);
    realacc = [realacc;nanmean(real_acc)];
    realacc2 = [realacc2;real_acc];
    daydat = [daydat;id*ones(size(real_acc))];
    shuffacc = [shuffacc;nanmean(shuff_acc)];
    shuffacc2 = [shuffacc2;shuff_acc];
    
    end
    numdaycells(id) = size(pfc,2);
    disp(['Done with day ' num2str(id) ' p ' num2str(pid(id)) ' # cells: ' num2str(size(pfc,2))])
    
    
        
%     figure; histogram(nanmean(shuff_acc),'FaceColor','k'); 
%     hold on; 
%     yl = get(gca,'ylim');
%     plot([nanmean(real_acc) nanmean(real_acc)],yl,'r-','LineWidth',3)
%     text(nanmean(real_acc),yl(2)*.8,['      p = ' num2str(round(pid(id),2,'significant'))])
%     title(['Day ' num2str(id)])
end

p = (sum(nanmean(shuffacc)>=nanmean(realacc))+1)/(size(shuffacc,2)+1);
figure; histogram(nanmean(shuffacc),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realacc) nanmean(realacc)],yl,'r-','LineWidth',3)
text(nanmean(realacc),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title(['Session Averages ' savelab])
% helper_saveandclosefig('E:\XY_matdata\AllDays\matfiles\RP_predictbehavior_HP_AllReplays')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\fixedbug\discrim\allreplays_armon_and_replayarm_given\discrim_pfcfwdreplays_effect_PFCsessionmeans_new1_leaveoneout_' num2str(numshuff) '_' savelab])

p = (sum(nanmean(shuffacc2)>=nanmean(realacc2))+1)/(size(shuffacc2,2)+1);
figure; histogram(nanmean(shuffacc2),'FaceColor','k'); 
hold on;
yl = get(gca,'ylim');
plot([nanmean(realacc2) nanmean(realacc2)],yl,'r-','LineWidth',3)
text(nanmean(realacc2),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title(['All Replays ' savelab])
% helper_saveandclosefig('E:\XY_matdata\AllDays\matfiles\RP_predictbehavior_PFC_AllReplays')
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\fixedbug\discrim\allreplays_armon_and_replayarm_given\discrim_pfcfwdreplays_effect_PFCallreplays_new1_leaveoneout_' num2str(numshuff) '_' savelab])

save(['E:\XY_matdata\AllDays\matfiles\RP_predictbehavior_PFC_leaveoneout_' num2str(numshuff) '_' savelab '.mat'],'shuffacc','realacc','shuffacc2','realacc2','daydat','numdaycells')


touse = numdaycells>3;
shuffacc = shuffacc(touse,:);
realacc = realacc(touse);
p = (sum(nanmean(shuffacc)>=nanmean(realacc))+1)/(size(shuffacc,2)+1);
figure; histogram(nanmean(shuffacc),'FaceColor','k');
hold on;
yl = get(gca,'ylim');
plot([nanmean(realacc) nanmean(realacc)],yl,'r-','LineWidth',3)
text(nanmean(realacc),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title(['Session Averages, All days with 4+ Cells ' savelab])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\fixedbug\discrim\allreplays_armon_and_replayarm_given\discrim_pfcfwdreplays_effect_PFCsessionmeans_new1_leaveoneout_' num2str(numshuff) '_4pluscells_' savelab])

touse2 = ismember(daydat,find(touse));
shuffacc2 = shuffacc2(touse2,:);
realacc2 = realacc2(touse2);

p = (sum(nanmean(shuffacc2)>=nanmean(realacc2))+1)/(size(shuffacc2,2)+1);
figure; histogram(nanmean(shuffacc2),'FaceColor','k');
hold on;
yl = get(gca,'ylim');
plot([nanmean(realacc2) nanmean(realacc2)],yl,'r-','LineWidth',3)
text(nanmean(realacc2),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title(['All Replays, All days with 4+ Cells ' savelab])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ArmCombinations\fixedbug\discrim\allreplays_armon_and_replayarm_given\discrim_pfcfwdreplays_effect_PFCallreplays_new1_leaveoneout_' num2str(numshuff) '_4pluscells_' savelab])