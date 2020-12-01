function transferdata
cd('D:\XY_matdata\TransferSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)    
    dirname = d2(idir).name(1:2);
    dirrun = str2double(d2(idir).name(end-4));
    cd('D:\XY_matdata\TransferSessions\')
    load(d2(idir).name,...
        'Rat_Name','Date','Track_Type','Run','params',...
        'behave_change_log','behave_ind','behavior','error_correct','headingarm',...
        'laps_coverspace','laps_singlepass','laps_twoarms')
    if dirrun~=Run
        error('Not the right run')
    end
    
    
    cd(['D:\XY_matdata\' dirname])
    load(['Run' num2str(Run) '_HP_SlowGamma.mat'])
    load(['Run' num2str(Run) '_HP_Ripple.mat'])
    load(['Run' num2str(Run) '_HP_Theta.mat'])
    load('Spike_Data','Tetrode_Cell_IDs','Spike_Data','HP_Cells','Inhibitory_Neurons','NAc_Cells')
    load('Run_Index')
    load('Velocity_Data')
    load('Field_Data_Directional_byrun.mat')
    load('Experiment_Information','LFP_Electrodes')
    
    InFR = Field_Data(:,:,2,Run)';
    OutFR = Field_Data(:,:,1,Run)';
    
    load('Position_Data.mat','Position_Data')
    pos = Position_Data(Run_Index(:,Run),:);
    vel = Velocity_Data(Run_Index(:,Run),:);

    [~,~,i] = histcounts(Spike_Data(:,1),Position_Data(:,1));
    i(i==0) = NaN;
    i = fillmissing(i,'nearest');
    
    ind = Run_Index(i,Run);
    Spike_Data1 = Spike_Data(ind,:);
    
    spikedata = [Spike_Data1(:,1:2) NaN(size(Spike_Data1,1),1)];
    for icell = 1:max(Spike_Data1(:,2))
        spikedata(spikedata(:,2) ==icell,3) = Tetrode_Cell_IDs(Tetrode_Cell_IDs(:,2)==icell,1);
    end
    hp_cells = HP_Cells;
    other_cells = NAc_Cells;
    hpinterneurons = hp_cells(ismember(hp_cells,Inhibitory_Neurons));

%     if sum(sum(Tetrode_Cell_IDs~=TT))>0
%         disp('Tetrode numbers dont match')
%     end

    clear Tetrode_Cell_IDs HP_Cells
    
    cd('D:\XY_matdata\AllSessions\')
    
%     if exist(['D:\XY_matdata\AllSessions\' Rat_Name '_' num2str(Date) '_Run' num2str(Run) '.mat'],'file')==0
    save(['D:\XY_matdata\AllSessions\' Rat_Name '_' num2str(Date) '_Run' num2str(Run) '.mat'],'dirname',...
            'Rat_Name','Date','pos','vel','Track_Type','Run','params',...
            'behave_change_log','behave_ind','behavior','error_correct','headingarm',...
            'laps_coverspace','laps_singlepass','laps_twoarms','LFP_Electrodes','HP_SlowGamma'...
            ,'HP_Ripple','HP_Theta','armpos','spikedata'...
            ,'hp_cells','hpinterneurons','other_cells','InFR','OutFR')        
%     end
    disp(['done with day ' num2str(idir)])
end
    
    