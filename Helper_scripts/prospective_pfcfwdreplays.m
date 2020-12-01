function prospective_pfcfwdreplays(thisdir,label,ilab)
disp(['Start prospective_pfcfwdreplays ' label])
numshuff = 1000;
load(thisdir,[label '_replay_singlebothjoint'],[label '_Cand_sig_modu_include'], ...
    [label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],...    
    [label  '_replay_shuffle_p']...
    ,[label  '_replay_replayarm'],[label '_replay_stnd'],'pos','armpos',...
    [label  '_replay_corr'],[label  '_replay_dirbias'],[label  '_replay_corr'],[label  '_replay_maxjump'],[ label  '_replay_armcov'])


eval(['PFCreplayspikes_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFCreplayspikes_list = ' label '_PFCreplayspikes_list;'])
% eval(['PFCreplayspikes_binned = ' label '_PFCcandspikes_binned;'])
% eval(['PFCreplayspikes_list = ' label '_PFCcandspikes_list;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
eval(['dirbias = ' label '_replay_dirbias;'])
eval(['corrs = ' label '_replay_corr;'])
eval(['mj = ' label  '_replay_maxjump;'])
eval(['wc = ' label  '_replay_corr;'])

eval(['armcov1 = ' label  '_replay_armcov;'])

eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
clear([label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_Cand_sig_modu_include'])

% pfc = other_cells; clear other_cells


% 1 is pure inbound, -1 is pure outbound
% dirbias = dirbias(touse);

% outbound is positive, inbound is negative
% corrs = corrs(touse);

% fwd = (dirbias>0 & corrs<0) | (dirbias<0 & corrs>0);
% fwd = (dirbias<0 & corrs>0);
% fwd = corrs>0;
% rev = (dirbias>0 & corrs>0) | (dirbias<0 & corrs<0);
% fwd = true(length(dirbias),1);
% 
% if ilab==1    
%     touse = ~(singlebothjoint==3);
%     Event(singlebothjoint==3) = [];
%     replayarm(singlebothjoint==3) = [];
% elseif ilab ==2
%     touse = ~(singlebothjoint==3 | singlebothjoint==0);
%     Event(singlebothjoint==3 | singlebothjoint==0) = [];
%     replayarm(singlebothjoint==3 | singlebothjoint==0) = [];   
% elseif ilab ==3
%     touse = ~(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0);
%     Event(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0) = [];
%     replayarm(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0) = [];    
% elseif ilab == 4
%     touse = ~(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0 | abs(corrs)<.5);
%     Event(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0 | abs(corrs)<.5) = [];
%     replayarm(singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0 | abs(corrs)<.5) = [];    
% end

% replaylab = {'AllArmEvents';'All.3Events';'AllSigEvents';'All.3SigEvents';'All.5SigEvents';'All.5Events'}; 
% % replaylab = {'AllArmEvents';'All.3Events';'All.3SigEvents';'All.3SigEvents_JD.6';'All.5SigEvents';'All.6SigEvents_JD.4'};
% replaylab = {'AllArmEvents';'Wc3ArmCov3';'WC4ArmCov5mj7';'WC5ArmCov5mj4'}; %;'All.5SigEventsArmcov5'}; 
if ilab==1    
    indd = singlebothjoint==3;
elseif ilab ==2
    indd = singlebothjoint==3 | abs(wc)<.3 | armcov1<.3;
%     indd = singlebothjoint==3 | abs(wc)<.3 | armcov1<.3; % singlebothjoint==0;
elseif ilab == 3
    indd = singlebothjoint==3 | abs(wc)<.4 | armcov1<.5 | mj>.7;
%     indd = singlebothjoint==3 | shuffle_p>=.05 | armcov1<.3;
elseif ilab == 4
    indd = singlebothjoint==3 | abs(wc)<.5 | armcov1<.5 | mj>.4;
%     indd = singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0;
elseif ilab == 5
    indd = singlebothjoint==3 | shuffle_p>=.05 | singlebothjoint==0 | abs(wc)<.5; 
elseif ilab == 6
    indd = singlebothjoint==3 | singlebothjoint==0 | abs(wc)<.5; 
end
 
% indd = singlebothjoint==3 | singlebothjoint==0;
indd = singlebothjoint==3;

Event(indd) = [];
replayarm(indd) = [];

% fwd(singlebothjoint==3 | shuffle_p>=.05) = [];
% Event(singlebothjoint==3 | shuffle_p>=.05) = [];
% replayarm(singlebothjoint==3 | shuffle_p>=.05) = [];
% fwd(singlebothjoint==3) = []; %changed 2/4/19
% Event(singlebothjoint==3) = []; %changed 2/4/19
% replayarm(singlebothjoint==3) = []; %changed 2/4/19

% touse = fwd;
PFCreplayspikes_binned = PFCreplayspikes_binned(:,:,~indd);
% Event = Event(touse,1);
% replayarm = replayarm(touse);

[~,~,i] = histcounts(Event,pos(:,1));
% Event(i==0) = [];
replayarm(i==0) = [];
armon = armpos(i(i~=0));

modind1 = [0 .2];
baseind1 = [-.5 -.1];
binsize = .02;
window = [-2 2];
ind = [window(1):binsize:window(2)];
modind = ind>=modind1(1) & ind<=(modind1(2));
baseind = ind>=baseind1(1) & ind<=baseind1(2);

% cellstouse = Cand_sig_modu_include(:,3)==1 & Cand_sig_modu_include(:,1)==1 & Cand_sig_modu_include(:,2)>0;
% % cellstouse = Cand_sig_modu_include(:,3)==1 & Cand_sig_modu_include(:,1)==1 & Cand_sig_modu_include(:,2)<0;
% % cellstouse = Cand_sig_modu_include(:,2)<0;
% % cellstouse = Cand_sig_modu_include(:,2)>0;
% % cellstouse = true(size(Cand_sig_modu_include,1),1);
% % cellstouse = Cand_sig_modu_include(:,3)==1 & Cand_sig_modu_include(:,1)==1;

pfcfwdreplaymodu_armXpos = NaN(size(PFCreplayspikes_binned,1),3,2);
pfcfwdreplaymodu_armXpos_shuff = NaN(size(PFCreplayspikes_binned,1),3,2,numshuff);
for iarm = 1:3
    otharms = setdiff(1:3,iarm);
    IND = armon ==iarm;
    if sum(IND)>0
    for ioth = 1:2
        repind = IND & replayarm==otharms(ioth);
        if sum(repind)>1
            pfcfwdreplaymodu_armXpos(:,iarm,ioth) = nanmean(squeeze(nanmean(PFCreplayspikes_binned(:,modind,repind),2)./range(modind1))...
                -squeeze(nanmean(PFCreplayspikes_binned(:,baseind,repind),2)./range(baseind1)),2);    
        end
    end    
    replayarm2 = replayarm;
    replayarm3 = replayarm(IND);
    allperms = NaN(length(replayarm2),numshuff);
    
    for ii = 1:numshuff        
        replayarm2(IND) = replayarm3(randperm(sum(IND)));
        for ioth = 1:2
            repind = IND & replayarm2==otharms(ioth);
            if sum(repind)>1
                while sum(replayarm~=replayarm2)==0 || sum(sum(allperms(:,1:ii-1)~=replayarm2)==0)>0
                    replayarm2(IND) = replayarm3(randperm(sum(IND)));
                    repind = IND & replayarm2==otharms(ioth);
                    disp(['ReplayListSame ii = ' num2str(ii)])
                end
                allperms(:,ii) = replayarm2;
               if sum(repind)>0
                pfcfwdreplaymodu_armXpos_shuff(:,iarm,ioth,ii) = nanmean(squeeze(nanmean(PFCreplayspikes_binned(:,modind,repind),2)./range(modind1))...
                    -squeeze(nanmean(PFCreplayspikes_binned(:,baseind,repind),2)./range(baseind1)),2);    
               end 
            end
        end
    end
    end
end
pfcfwdreplaymodu_armXpos = reshape(pfcfwdreplaymodu_armXpos,[size(pfcfwdreplaymodu_armXpos,1) 6]); 
pfcfwdreplaymodu_armXpos_shuff = reshape(pfcfwdreplaymodu_armXpos_shuff,[size(pfcfwdreplaymodu_armXpos_shuff,1) 6 numshuff]); 
eval([label '_pfcfwdreplaymodu_armXpos = pfcfwdreplaymodu_armXpos;'])
eval([label '_pfcfwdreplaymodu_armXpos_shuff = pfcfwdreplaymodu_armXpos_shuff;'])

save(thisdir,[label '_pfcfwdreplaymodu_armXpos'],[label '_pfcfwdreplaymodu_armXpos_shuff'],'-append')
disp(['Done with prospective_pfcfwdreplays ' label])