function [Th] = add_time_since_last_reward(thisdir)
load(thisdir,'laps_singlepass','behave_change_log','times_armon_thetaof_headingarm_lap_thetahalf_all','pos')
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
binsize = .75;
binz = binsize/(1/60);
timerevdat = NaN(size(laps_singlepass));
for ilap = 1:max(laps_singlepass)
   leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'last'); %1 leave middle 5 leave platform
   lapind = leave:find(laps_singlepass==ilap,1,'last');
   
   
   
   lapindb = ceil((1:length(lapind))/binz);
   timerevdat(lapind) = lapindb;
end
[~,i] = histc(Th(:,2),pos(:,1));
Th(:,10) = timerevdat(i);