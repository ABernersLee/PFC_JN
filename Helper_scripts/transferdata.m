
cd('D:\XY_matdata\TransferSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)
    dirname = d2(idir).name(1:2);
    dirrun = str2double(d2(idir).name(end-4));
    load(d2(idir).name,...
        'Rat_Name','Date','Track_Type','rawpos','Run','rawspikedata','params','pos','vel','armpos','cm_conv',...
        'dirdat','linposcat','linposcatnan','linposnorm','behave_change_log','behave_ind','behavior','error_correct','headingarm',...
        'laps_coverspace','laps_singlepass','laps_twoarms','hp_cells','hpinterneurons','other_cells','spikedata','Cell_Number')
    if dirrun~=Run
        error('Not the right run')
    end
    if exist(['D:\XY_matdata\AllSessions\' Rat_Name '_' num2str(Date) '_Run' num2str(Run) '.mat'],'file')==0
        save(['D:\XY_matdata\AllSessions\' Rat_Name '_' num2str(Date) '_Run' num2str(Run) '.mat'],'dirname',...
            'Rat_Name','Date','Track_Type','rawpos','Run','rawspikedata','params','pos','vel','armpos','cm_conv',...
            'dirdat','linposcat','linposcatnan','linposnorm','behave_change_log','behave_ind','behavior','error_correct','headingarm',...
            'laps_coverspace','laps_singlepass','laps_twoarms','hp_cells','hpinterneurons','other_cells','spikedata','Cell_Number')
    end
end
    
    