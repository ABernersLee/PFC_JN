function make_prospective_theta_triggered(thisdir)

load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','spikedata','other_cells');
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
clear times_armon_thetaof_headingarm_lap_thetahalf_all
%same as for replay
modind1 = [0 .2];
% baseind1 = [-.5 -.1];


spikes = spikedata(ismember(spikedata(:,2),other_cells),:);

tha = NaN(length(other_cells),3,2);
for iarm = 1:3
    otherarms = setdiff(1:3,iarm);
    for ioth = 1:2
        ind = Th(:,3)==iarm & Th(:,4)==otherarms(ioth); % & Th(:,end)>.3;
        if sum(ind)>0
            modind = [Th(ind,2)+modind1(1) Th(ind,2)+modind1(2)];
%             baseind = [Th(ind,2)+baseind1(1) Th(ind,2)+baseind1(2)];

            thall = NaN(size(modind,1),length(other_cells));
            for ievent = 1:size(modind,1)
                modspikes = spikes(spikes(:,1)>=modind(ievent,1) & spikes(:,1)<modind(ievent,2),2);
%                 basespikes = spikes(spikes(:,1)>=baseind(ievent,1) & spikes(:,1)<baseind(ievent,2),2);
                ms = histc(modspikes,other_cells);
                if size(ms,1)==1; ms = ms'; end
%                 bs = histc(basespikes,other_cells);
%                 if size(bs,1)==1; bs = bs'; end
%                 thall(ievent,:) = (ms./range(modind1))-(bs./range(baseind1));
                thall(ievent,:) = (ms./range(modind1));
            end
            tha(:,iarm,ioth) = nanmean(thall);
        end
        
    end
   
end
prospectiveTheta_armcombs = reshape(tha,[size(tha,1) 6]);
 
%need to make other types for controls too
save(thisdir,'prospectiveTheta_armcombs','-append')






