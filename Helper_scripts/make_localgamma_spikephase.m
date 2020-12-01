function make_localgamma_spikephase(type)
cd('D:\XY_matdata\AllSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)
    cd('D:\XY_matdata\AllSessions\')
    load(d2(idir).name,'dirname','TT','Run','spikedata','LFP_Electrodes')
    
    gamma_localspikephase = NaN(size(spikedata,1),3);
    
    cd(['D:\XY_matdata\' dirname])
    tts = unique(TT(:,1));
    for itt = 1:length(tts)
        cellsT = TT(tts(itt)==TT(:,1),2);
        index = ismember(spikedata(:,2),cellsT);
        s = spikedata(index,1);
        posT = (tts(itt)*4)-3:tts(itt)*4;
        icsc = posT(ismember(posT,LFP_Electrodes));
        icsc = icsc(1); % if there are more than one channel for the tetrode, just pick the first, if one is bad go back and delete it later
        
        if type == 1
            load(['Run' num2str(Run) '_LowGamma_CSC' num2str(icsc) '.mat'])
            eval(['p = Run' num2str(Run) '_LowGamma_CSC' num2str(icsc) ';'])
            clear(['Run' num2str(Run) '_LowGamma_CSC' num2str(icsc)])      
        elseif type == 2
            load(['Run' num2str(Run) '_HighGamma_CSC' num2str(icsc) '.mat'])
            eval(['p = Run' num2str(Run) '_HighGamma_CSC' num2str(icsc) ';'])
            clear(['Run' num2str(Run) '_HighGamma_CSC' num2str(icsc)])      
        end
        
        [~,~,bin] = histcounts(s(:,1),p(:,1));
        if sum(bin==0)>0
            indd = [true(ceil(length(bin)/2),1); false(length(bin)-ceil(length(bin)/2),1)];
            bin(bin==0 & indd)= 1;
            bin(bin==0 & ~indd) = bin(find(bin~=0,1,'last'));
        end
        gamma_localspikephase(index,:) = p(bin,2:4); % power, z-scored power, phase
    end
    
    gamma_localspikephase = [spikedata(:,1:2) gamma_localspikephase]; %spiketime, cell #, % local gamma power, z-scored gamma power, gamma phase
    
    cd('D:\XY_matdata\AllSessions\')
    if type == 1
        slowgamma_localspikephase = gamma_localspikephase;
        save(d2(idir).name,'slowgamma_localspikephase','-append')
    elseif type == 2
        fastgamma_localspikephase = gamma_localspikephase;
        save(d2(idir).name,'fastgamma_localspikephase','-append')
    end
end