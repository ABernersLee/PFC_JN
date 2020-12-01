function ProcessData(dirs)

cd(dirs.homedir)
d2 = dir('*.mat');

for id =  1:size(d2,1) %1:size(d2,1)
    cd(dirs.homedir)
    load([d2(id).name],'params','rawpos','rawspikedata','rundat')  
%     load([d2(id).name],'params','rundat','armpos','linposcat','behavior','behave_change_log','pos') 
    
%     adds arms and armslength to params and makes pos that has had incorrect points deleted and linearly interpolated points replacing them
    [pos,params,vel] = GUI_make_arms_cleanup_posdata_day(rawpos,params,rundat);
    params.armslength = [161;81;81]; % just double checking, this is XW data only
%     params.Run = NaN;
    
    %linearizes position data and assigns arms
    [armpos,linposcat,linposnorm,linposcatnan,dirdat,cm_conv] = GUI_make_linear_position(pos,params);
    
    
    %makes behavior epoch marticies (middle, stem, neck, lick area)
    [behavior, behave_change_log, behave_ind] = GUI_make_behavior_day(linposcat,armpos,linposcatnan,params);
        
        
    %makes laps, a few different ways for track type 2, see script
    [laps_coverspace,laps_twoarms,laps_singlepass,headingarm,error_correct] = GUI_makelaps_day(behavior, behave_change_log, pos,armpos, linposcat,params,rundat);
    
    
    if 0
        %indexes the position bins of each spike, seperates types of neurons
        [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params);
    end
    
    save([dirs.spikedatadir params.ident],'pos','params','vel', ...
        'armpos','dirdat','linposcat','linposnorm','linposcatnan','cm_conv', ...
        'behavior','behave_change_log','behave_ind', ...
        'laps_coverspace','laps_twoarms','laps_singlepass','headingarm','error_correct','-append')
    
%     save([dirs.spikedatadir params.ident],'spikedata','hp_cells','other_cells','hpinterneurons','-append')

%     save([dirs.spikedatadir params.ident],'laps_coverspace','laps_twoarms','laps_singlepass','-append')
    
    clearvars -except dirs d2 id
    disp(['Done with processing for day ' num2str(id)])
end