function run_all_linear_classifiers(thisdir,label,types,id,runcontrols,rundownsample)

load(thisdir,[label  '_replay_replayarm'],[label  '_replay_singlebothjoint'],[label '_replay_dirbias'], ...
[label '_replay_corr'],[label '_replay_stnd'],'armpos','pos',[label  '_replay_shuffle_p']...
    ,[label  '_replay_corr'],[label  '_replay_maxjump'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['dirbias = ' label '_replay_dirbias;'])
eval(['corr = ' label '_replay_corr;'])
eval(['st = ' label '_replay_stnd(:,1);'])
eval(['shuffle_p = ' label  '_replay_shuffle_p;'])
eval(['mj = ' label  '_replay_maxjump;'])
eval(['wc = ' label  '_replay_corr;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])

[~,~,i] = histcounts(st,pos(:,1));
if sum(i==0)>0
    inx = find(i==0);
    i(inx(inx<(size(pos,1)/2))) = find(i~=0,1,'last');
    i(inx(inx>(size(pos,1)/2))) = find(i~=0,1,'first');
end
replaypos = armpos(i);


for itype = types
%     if 0
    disp(['Start run_all_linear_classifiers Arm '])
    linear_classifier_of_replay_withPFC(thisdir,label,itype,replayarm,'Arm',rundownsample)
    disp(['Done run_all_linear_classifiers Arm '])
%     end
    if runcontrols
        disp(['Start run_all_linear_classifiers Arm controls'])
        linear_classifier_of_replay_withPFC_controls(thisdir,label,replayarm,replaypos,'ArmControls',id)
        disp(['Done run_all_linear_classifiers Arm controls'])
    end
%     Labels= (dirbias(singlebothjoint~=3 & shuffle_p<.05)>0)+1;
%     linear_classifier_of_replay_withPFC(thisdir,label,itype,Labels,'Dir')
%     disp(['Done run_all_linear_classifiers Dir ' num2str(itype) ' ' label])
%     
%     Labels= replaypos(singlebothjoint~=3 & shuffle_p<.05);
%     linear_classifier_of_replay_withPFC(thisdir,label,itype,Labels,'RatPos')
%     disp(['Done run_all_linear_classifiers RatPos ' num2str(itype) ' ' label])
end