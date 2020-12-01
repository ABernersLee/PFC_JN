function make_HPtheta_spikephase2(thisdir)
disp('Start make_HPtheta_spikephase2')
load(thisdir,'spikedata','pos','HP_Theta')     

s = spikedata(:,1:3); clear spikedata
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
HPTheta_spikephase = [s(:,1:3) p(bin,2:4)]; % spike time, spike cell % theta power, theta z-scored power, theta phase
HPTheta_spikephase(s(:,1)<p(1,1),:) = [];
HPTheta_spikephase(s(:,1)>p(end,1),:) = [];

save(thisdir,'HPTheta_spikephase','-append')
disp('Done with make_HPtheta_spikephase2')