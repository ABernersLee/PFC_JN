function alldays_alldat_table(dirs,igroup,savefolder,toplot)
%not set up with igroup, tables are not a good format to save out

cd('E:\XY_matdata');
load('dirs.mat')

cd(dirs.homedir)
d2 = dir('*.mat');
alldat = NaN(size(d2,1),17);
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'RP_CandEventTimes','hp_cells','other_cells','other_cells_touse','hpinterneurons','RP_replay_shuffle_p','pos','laps_singlepass','laps_coverspace','rundat'...
        ,'RP_replay_singlebothjoint') %,'SD_replay_singlebothjoint','SD_replay_shuffle_p')
%     error2 = NaN;
    error2 = Behavior_Decoding_Simple(thisdir,toplot,savefolder);       
    other_cells = other_cells(other_cells_touse(:,igroup));
%     SD_Cand_sig_modu_include = SD_Cand_sig_modu_include(other_cells_touse(:,igroup));
    
    tm = 0;
    for irun = 1:size(rundat,2)
        tm = tm+sum(diff(pos(rundat(:,irun)==1,1)));
    end
    
    armacc = get_behavior_accuracy(thisdir);                     
    sz = floor(size(armacc,1)/2);
    h1 = mean(armacc(1:sz,3));
    h2 = mean(armacc(size(armacc,1)-sz+1:end,3));
    learnrate = ((h2-h1)/h1)*100;
    ne = armacc(~isnan(armacc(:,4)),4);
    sz = floor(size(ne,1)/2);
    h1 = mean(ne(1:sz));
    h2 = mean(ne(size(ne,1)-sz+1:end));
    learnrate2 = ((h2-h1)/h1)*100;
    
    load(thisdir,'RP_Cand_sig_modu_include','RP_pSSDarm') %,'SD_Cand_sig_modu_include','SD_pSSDarm')
    RP_Cand_sig_modu_include = RP_Cand_sig_modu_include(other_cells_touse(:,igroup),:);
    RP_pSSDarm = RP_pSSDarm(other_cells_touse(:,igroup),:);
    
    rps = [sum(RP_Cand_sig_modu_include(:,[3 1])) sum(RP_pSSDarm<.05)];
%     sds = [sum(SD_Cand_sig_modu_include(:,[3 1])) sum(SD_pSSDarm<.05)];
    
    alldat(id,:) = [length(setdiff(hp_cells,hpinterneurons)) length(other_cells) ...
        tm/60 size(rundat,2) max(laps_coverspace) max(laps_singlepass) ...
        mean(armacc(:,3)) nanmean(armacc(:,4)) learnrate learnrate2...
        size(RP_CandEventTimes,1) sum(RP_replay_shuffle_p(RP_replay_singlebothjoint~=3)<.05) ...
        mean(error2) median(error2) rps];
end
% T = table(alldat(:,1),alldat(:,2),alldat(:,3),alldat(:,4),alldat(:,5),alldat(:,6)...
%     ,alldat(:,7),alldat(:,8),alldat(:,9),alldat(:,10),alldat(:,11),alldat(:,12),...
%     alldat(:,13),alldat(:,14),alldat(:,15),alldat(:,16),alldat(:,17),alldat(:,18)...
%     ,alldat(:,19),alldat(:,20));
% T.Properties.VariableNames = {'hp_cells';'pfc_cells';'totaltime_min';'runs';'laps_coverspace';'laps_singlepass';...
%     'rewardrate';'alternation_accuracy';'learning';'learning_alternation';...    
%     'SigRippleReplay';'SigSDReplay';'BehaveDecodingMeanError';'BehaveDecodingMedianError'...
%     ;'RP_included';'RP_SigMod';'RP_SigArmMod';'SD_included';'SD_SigMod';'SD_SigArmMod'};
% 
% alldat = round(alldat,2,'significant');

C2 = alldat(1,:)'; C4 = alldat(2,:)'; C5 = alldat(3,:)'; D1 = alldat(4,:)'; D2 = alldat(5,:)'; D3 = alldat(6,:)'; D6 = alldat(7,:)';
E1 = alldat(8,:)'; E3 = alldat(9,:)'; E4 = alldat(10,:)'; I1 = alldat(11,:)';
% Types = {'hp_cells';'pfc_cells';'totaltime_min';'runs';'laps_coverspace';'laps_singlepass';...
%     'rewardrate';'alternation_accuracy';'learning';'learning_alternation';...    
%     'SigRippleReplay';'SigSDReplay';'BehaveDecodingMeanError';'BehaveDecodingMedianError'...
%     ;'RP_included';'RP_SigMod';'RP_SigArmMod';'SD_included';'SD_SigMod';'SD_SigArmMod'};
Types = {'hp_cells';'pfc_cells';'totaltime_min';'runs';'laps_coverspace';'laps_singlepass';...
    'rewardrate';'alternation_accuracy';'learning';'learning_alternation';...    
    'CandEvents';'SigRippleReplay';'BehaveDecodingMeanError';'BehaveDecodingMedianError'...
    ;'RP_included';'RP_SigMod';'RP_SigArmMod'};
T = table(C2,C4,C5,D1,D2,D3,D6,E1,E3,E4,I1,'RowNames',Types);
% T.Properties.ColumnNames = {'C2';'C4';'C5';'D1';'D2';'D3';'D6';'E1';'E3';'E4';'I1'};
% save('E:\XY_matdata\Figures\ForPaper\Basics\alldays_alldat_20190126.mat','alldat','T','d2')

if ~isfolder([savefolder '\Basics\'])
    mkdir([savefolder '\Basics\'])
end

save([savefolder '\Basics\alldays_alldat.mat'],'alldat','T','d2')
f = uifigure; uit = uitable(f,'Data',T);
set(uit,'OuterPosition',[20          20        1112         435])
set(uit,'InnerPosition',[20          20        1112         435])
set(f,'Position',[535         238        1166         481])
uit.BackgroundColor = [1 1 1;1 1 1];

% saveas(f,['E:\XY_matdata\Figures\ForPaper\Basics\Table.fig'],'fig')
% saveas(f,['E:\XY_matdata\Figures\ForPaper\Basics\Table.tiff'],'tiff')


