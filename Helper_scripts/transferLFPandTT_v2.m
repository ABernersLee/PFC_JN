function transferLFPandTT
cd('D:\XY_matdata\AllSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)
    cd('D:\XY_matdata\AllSessions\')
    load(d2(idir).name,'dirname','Run')
    
    cd(['D:\XY_matdata\' dirname])
    load(['Run' num2str(Run) '_HP_SlowGamma.mat'])
    load(['Run' num2str(Run) '_HP_Ripple.mat'])
    load(['Run' num2str(Run) '_HP_Theta.mat'])
    load('Spike_Data','Tetrode_Cell_IDs','Spike_Data','HP_Cells','Inhibitory_Neurons','NAc_Cells')
    load('Run_Index')
    load('Field_Data_Directional_byrun.mat')
    load('Experiment_Information','LFP_Electrodes')
    
    InFR = Field_Data(:,:,2,Run)';
    OutFR = Field_Data(:,:,1,Run)';
    
    Spike_Data = Spike_Data(Run_Index(:,Run),:);
    spikedata = [Spike_Data(:,2) NaN(size(Spike_Data,1),1)];
    for icell = 1:max(Spike_Data(:,2))
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
    save(d2(idir).name,'LFP_Electrodes','HP_SlowGamma','HP_Ripple','HP_Theta','spikedata'...
        ,'hp_cells','hpinterneurons','other_cells','InFR','OutFR','-append')
    disp(['done with day ' num2str(idir)])
end


