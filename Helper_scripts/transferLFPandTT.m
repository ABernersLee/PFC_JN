cd('D:\XY_matdata\AllSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)
    cd('D:\XY_matdata\AllSessions\')
    load(d2(idir).name,'dirname','Run','hp_cells','rawspikedata')
    
    cd(['D:\XY_matdata\' dirname])
    load(['Run' num2str(Run) '_HP_SlowGamma.mat'])
    load(['Run' num2str(Run) '_HP_Ripple.mat'])
    load(['Run' num2str(Run) '_HP_Theta.mat'])
    load('Spike_Data','Tetrode_Cell_IDs','HP_Cells')
    load('Experiment_Information','LFP_Electrodes')
    
%     if sum(HP_Cells~=[hp_cells'])>0
%         disp('Cell numbers dont match')
%     end

    allcells = unique(rawspikedata(:,2));
    TT = NaN(length(allcells),2);
    for icell = 1:length(allcells)
        TT(icell,1) = unique(rawspikedata(rawspikedata(:,2)==allcells(icell),3));
        TT(icell,2) = allcells(icell);
    end

%     if sum(sum(Tetrode_Cell_IDs~=TT))>0
%         disp('Tetrode numbers dont match')
%     end

    clear Tetrode_Cell_IDs HP_Cells
    
    cd('D:\XY_matdata\AllSessions\')
    save(d2(idir).name,'TT','LFP_Electrodes','HP_SlowGamma','HP_Ripple','HP_Theta','-append')
    disp(['done with day ' num2str(idir)])
end


