function RunAnalysis(dirs)

%% Directories
cd(dirs.homedir)
d2 = dir('*.mat');

%% Defining some aspects of the analysis
velcutoff = 5;
num_shuffles = 5000;
exampleplots = false; otherplots =true;
ilab = 1; %used 1 for JN paper
replaylab = {'AllArmEvents';'Wc3ArmCov3';'WC4ArmCov5mj7';'WC5ArmCov5mj4'};
grouplabels = {['AllCells_' replaylab{ilab}];['BestTTDays_' replaylab{ilab}];['BestDays_' replaylab{ilab}]};
smoothsize = 30; withmodpatch = 1;
std_cutoff = 2;
cutoff = 0;
igroup = 1;  %used both 1 and 3 for JN paper
if ~isfolder(['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}])    
    mkdir(['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}]) 
end
savefolder = ['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}];
%% figures and data within day
tic        
for id = 1:size(d2,1)
    try
    thisdir = d2(id).name;     
        
    
    %make directional fields
    make_directionalfields(thisdir,velcutoff,true,2);    
    
    %make Theta and Ripple mat files 
    make_HP_csc(dirs,thisdir)        
    
    %make Ripple Candidate events
    make_RPCandEvents_ripple(thisdir,std_cutoff)      
    
    %get theta phase for each spike
    make_HPtheta_spikephase2(thisdir)
    
    %decode posterior of whole session
    make_PosteriorfromCandEvents(thisdir,'RP')              
    
    %get replays from the posterior (using a modified version of XYW way)
    make_singleandjointreplay_fromPosterior(thisdir,'RP',num_shuffles)   
    
    %make two matricies that are pfc cell firing triggered off of the
    %replays identified in an easy to plot format
    make_PFC_replayeventtriggeredmat(thisdir,'RP')   
    
    %get significance of pfc cells being modulated by ripples
    make_PFC_candeventtriggeredmat(thisdir,'RP')
            
    %get the significance of differential modulation for each pfc cell    
    get_PFC_armtriggered_modusig(thisdir,'RP',ilab);  
    
    %plot the pfc cells triggered to different replays
    plot_PFC_ArmReplay_triggered(thisdir,'RP',withmodpatch,smoothsize,savefolder,id)      
    
    %plot the expected pfc cells triggered to different replays
    plot_PFC_ArmReplay_triggered_expected(thisdir,withmodpatch,smoothsize,savefolder,id,true);        
            
    %classify which arm is being replayed in hp with only pfc cells (and
    %controls to account for time within session, position of rat, and ripple power)
    run_all_linear_classifiers(thisdir,'RP',1,id,true)        
    
    %decodes theta sequences, adjusts theta to center around theta
    %sequence, then identifies any non-local theta sweeps
    Theta_Shifts_Sequences_Xcorr(thisdir,velcutoff,exampleplots,otherplots,savefolder);
    
    %make two matricies that are pfc cell firing triggered off of the
    %non-local theta sweeps identified in an easy to plot format
    make_PFC_thetaeventtriggeredmat(thisdir)
        
    %spatial preceition within replays
    make_ArmReplay_triggered_spatialprecision(thisdir,'RP')           
    make_ArmReplay_triggered_spatialprecisiondir(thisdir,'RP')         
            
    t = toc;
    disp(['!!! Done with Day ' num2str(id) ', finished in ' num2str(round(t/60,3,'significant')) ' minutes'])
    tic
     catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(id)
        disp('****************************************Error Occured')
    end
end
make_other_cells_touse(dirs)
make_OpenPlaceFields(dirs)

%% Figures across days

toplot = true; toplotcells = true; tosave = true;
for igroup = 1:3 %only used 1 and 3 in paper
    if ~isfolder(['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}])    
        mkdir(['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}]) 
    end
    savefolder = ['F:\XY_matdata\Figures\ForPaperReviews\' grouplabels{igroup}];
    
    
    %%%% Examples or other things that don't depend on stats of which cells to include
    if igroup==1       
         %Figure 1
        Figure1_PlaceFields(dirs,savefolder)    
        
        %Figure 2
        plot_ripples_vs_theta(dirs,savefolder,false) %toplot theta for individual cells
        
        %Figure 3
        Fig2_Examples(dirs,savefolder)
        
        %Figure 4
        Figure2_ExpectedVsReal(savefolder,dirs,igroup,1)
        
        %Figure 5        
        Figure3_AverageThetaSequenceAcrossDays(dirs,savefolder)  
        NonLocalSweeps_ThetaCycle_SummaryFigs(dirs,cutoff,savefolder) 
        
        %Other, in older versions of the paper    
%         alldays_alldat_table(dirs,igroup,savefolder,toplot)                
    end
    
    % Figure 2              
    plot_ripples_vs_theta(dirs,savefolder,false)        
    plot3PFC_arm_modulation('RP',igroup,savefolder)        
    
    % Figure 3
    proportion_replay_modulated_all('RP','',igroup,savefolder)
    
    % Figure 4
    Figure2_ExpectedVsReal(savefolder,dirs,igroup,1,0) %sig/not is 1/0 and then ReveiwerCommentVersion/Normal is 1/0
    Figure2_ExpectedVsReal(savefolder,dirs,igroup,0,0) %sig/not is 1/0 and then ReveiwerCommentVersion/Normal is 1/0
    Figure2_ExpectedVsReal(savefolder,dirs,igroup,1,1) %sig/not is 1/0 and then ReveiwerCommentVersion/Normal is 1/0
    Figure2_ExpectedVsReal(savefolder,dirs,igroup,0,1) %sig/not is 1/0 and then ReveiwerCommentVersion/Normal is 1/0
    
    % Other (in onlder versions of paper)
%     ArmReplay_triggered_spatialprecision(dirs,igroup,savefolder) 
%     ArmReplay_triggered_spatialprecisiondir(dirs,igroup,savefolder)
        
    if igroup~=2
        %Figure 4 and Results and Methods
        linear_classifier_of_replay_all('RP','Arm',1,'',igroup,savefolder,0)
        linear_classifier_of_replay_all('RP','Arm',1,'',igroup,savefolder,1)
        linear_classifier_of_replay_all('RP','Arm_ds',1,'',igroup,savefolder,0)
        linear_classifier_of_replay_all_controls('RP','ArmControls_time','',igroup,savefolder,0)
        linear_classifier_of_replay_all_controls('RP','ArmControls_armpos','',igroup,savefolder,0)
        linear_classifier_of_replay_all_controls('RP','ArmControls_ripple','',igroup,savefolder,0)
        linear_classifier_of_replay_all_controls('RP','ArmControls_replaylength','',igroup,savefolder,0)
    end
        
    %Figure 5
    if igroup==1        
        crosscov_thetatimes(dirs,toplot,toplotcells,tosave,velcutoff,1,igroup,savefolder)      
        crosscov_thetatimes(dirs,toplot,toplotcells,tosave,velcutoff,2,igroup,savefolder)    
    end    
    ThetaCycleHalf_pfc_SummaryFigs(dirs,'0904',igroup,savefolder,1)        
    ThetaCycleHalf_pfc_SummaryFigs(dirs,'0328',igroup,savefolder,2)
    
    %Figure 6 - Internal Fields plus reviewer comments
        if igroup~=2
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.04); % 8cm-8cm         
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,1.5,20,0.04); % 8cm-12cm
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,2.5,20,0.04); % 8cm-20cm  
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,4.5,20,0.04); % 8cm-32cm 
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,15,0.04); %6cm-6cm
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,1.5,15,0.04); %6cm-9cm
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,12,0.04); %5cm-5cm 
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,1.5,12,0.04); %5cm-7cm
            plotInternalErrorBins(savefolder)                          
            
            %across different delays
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.01); 
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.02);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.03);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.05);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.06);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.07);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.09);
            Figure3_Internal2d(dirs,igroup,savefolder,false,false,.04,5,4,5,1000,0.5,20,0.11);    
            plotPFCdelays(savefolder)
            
            Figure3_Internal2d_byphase(dirs,igroup,savefolder,false,false,5,2,4,1000);
        end
        
    %Figure 6 - Internal Fields cont.
    if igroup == 1
        Figure3_make_4sets_barplots(dirs,igroup,savefolder,toplotcells)
    else
        Figure3_make_4sets_barplots(dirs,igroup,savefolder,false)
    end
    
    %%%% Figure 7 - choice
    if igroup~=2
        Theta_Sweeps_Predict_Future_Arm(dirs,cutoff,savefolder,[0 Inf],igroup,1000) 
        decode_choice_by_thetaphase_all(dirs,savefolder,toplot,igroup,1000) 
    end
            
end

%Reviewer comments
OtherDays
Compare_GLM_PFC
Compare_GLM_HP
plot_PFC_ArmReplay_triggered_expected_ReviewerComment
get_PFC_fwdrev_dir_towards
ReviewerComment_CovarianceSimulations
ReviewerComment_ClassifyDownsample


