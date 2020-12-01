function MakeFigures(dirs,cutoff)

cd(dirs.homedir)
d2 = dir('*.mat');


for ii = 0:3
    plot_armcombs_relationship(dirs,'RP',ii,['_' num2str(cutoff) 'cutoff_new'])
end

armcombs_acrosspvals(dirs,cutoff)

% plot_lapbylap_PFC_all('RP',3)       
 

%(change from pie to bars)
proportion_replay_modulated_all('RP')


linear_classifier_of_replay_all('RP','Arm',1)
linear_classifier_of_replay_all_controls('RP','ArmControls_time')
linear_classifier_of_replay_all_controls('RP','ArmControls_armpos')

% linear_classifier_of_replay_all(label,'Dir',1)
% linear_classifier_of_replay_all(label,'RatPos',1)

plot3PFC_arm_modulation('RP')


% plot_Velocity_Vs_Value(label)



NonLocalSweeps_ThetaCycle_SummaryFigs
ThetaCycleHalf_pfc_SummaryFigs
NonLocalThetaDifference_SummaryFigs(dirs,cutoff)
Theta_Sweeps_Predict_Future_Arm(dirs,cutoff)
plot_MRVs(dirs)
% make_PosteriorfromThetaCandEvents_manypluspfc(thisdir)
fwdvsrev_pfc_plot(dirs)

%% Fig 1
%To Do: 
% -get ripples out and plot
% -get hp place cells out and plot spikes
% - do decoding over this period
Fig1_Examples


%% table

alldays_alldat_table

%%

%need to add plots to this
% Jadhav_Frank_Correlations_Ripple_PF(dirs)

%%
smoothsize = 30; withmodpatch = 1; 
for id =  1:size(d2,1)
    thisdir = d2(id).name;        
%     plot_PFC_Event_triggered(thisdir,'RP',withmodpatch,smoothsize)    
    plot_PFC_ArmReplay_triggered(thisdir,'RP',withmodpatch,smoothsize)
end
%     plot_PFC_Event_triggered(thisdir,'SD',withmodpatch,smoothsize)
%     if id>2
%         plot_prospective_FR(thisdir,5,8)        
%         plot_prospective_pfcfwdreplays(thisdir,'RP',smoothsize)
%         plot_prospective_pfcfwdreplays(thisdir,'SD',smoothsize)
%     end
%     

    
    
%     plot_PFC_ArmReplay_triggered(thisdir,'SD',withmodpatch,smoothsize)
%     disp(id)
%     for timetrig = 1:2
%         plot_PFC_JointReplay_triggered(thisdir,'RP',withmodpatch,smoothsize,timetrig)
%         plot_PFC_JointReplay_triggered(thisdir,'SD',withmodpatch,smoothsize,timetrig)
%     end        
% end
