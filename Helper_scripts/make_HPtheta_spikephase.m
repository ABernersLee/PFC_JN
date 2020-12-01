function make_HPtheta_spikephase

cd('D:\XY_matdata\AllSessions\')
d2 = dir('*.mat');

for idir = 1:size(d2,1)
    cd('D:\XY_matdata\AllSessions\')
    load(d2(idir).name,'dirname','TT','spikedata','HP_Theta','LFP_Electrodes','pos')        
    s = spikedata(:,1:2); clear spikedata
    p = HP_Theta; clear HP_Theta
    s(s(:,1)<pos(1,1) | s(:,1)>pos(end,1),:) = [];   
    p(p(:,1)<pos(1,1) | p(:,1)>pos(end,1),:) = [];
    clear pos
    
    
    [~,~,bin] = histcounts(s(:,1),p(:,1));
    
    if sum(bin==0)>0
        indd = [true(ceil(length(bin)/2),1); false(length(bin)-ceil(length(bin)/2),1)];
        bin(bin==0 & indd)= 1;
        bin(bin==0 & ~indd) = bin(find(bin~=0,1,'last'));
    end
    HPTheta_spikephase = [s(:,1:2) p(bin,2:4)]; % spike time, spike cell % theta power, theta z-scored power, theta phase
    HPTheta_spikephase(s(:,1)<p(1,1),:) = [];
    HPTheta_spikephase(s(:,1)>p(end,1),:) = [];
    
    cd('D:\XY_matdata\AllSessions\')
    save(d2(idir).name,'HPTheta_spikephase','-append')    
end