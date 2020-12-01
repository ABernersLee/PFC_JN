function combine_runs

cd('D:\XY_matdata\AllSessions\')
d2 = dir('*.mat');
ident = NaN(size(d2,1),1);
for id = 1:size(d2,1)
    ident(id) = str2double(d2(id).name(3:10));
end

[~,d] = unique(ident);
d = [d;length(ident)];
clear ident

for id = 4:length(d)-1
    
    
    lapnum = zeros(3,1);
    dirname2 = d2(d(id)).name(1:end-9);
    
    armpos1 = [];
    behave_change_log1 = [];
    behave_ind1 = [];
    behavior1 = [];
    dirdat1 = [];
    error_correct1 = [];
    headingarm1 = [];
    laps_coverspace1 = [];
    laps_singlepass1 = [];
    laps_twoarms1 = [];
    linposcat1 = [];
    linposcatnan1 = [];
    linposnorm1 = [];
    pos1 = [];
    rawspikedata1 = [];
    vel1 = [];
    HP_Ripple1 = [];
    HP_SlowGamma1 = [];
    HP_Theta1 = [];
    slowgamma_localspikephase1 =[];
%     HPTheta_spikephase1 = [];
    
    for irun = 1:(d(id+1)-d(id))
        if irun==1
            load([dirname2 '_Run' num2str(irun) '.mat'], 'Date', 'Rat_Name', 'dirname', 'Track_Type', ...
                'armpos', 'behave_change_log', 'behave_ind', 'behavior', 'cm_conv', 'dirdat', ...
                'error_correct', 'headingarm', 'laps_coverspace', ...
                'laps_singlepass', 'laps_twoarms', 'linposcat', 'linposcatnan', 'linposnorm', ...
                'params', 'pos', 'rawspikedata', 'vel', ...
                'HP_Ripple', 'HP_SlowGamma', 'HP_Theta', 'TT', 'armposindex', ...
                'slowgamma_localspikephase', 'HPTheta_spikephase')                            
        else
            load([dirname2 '_Run' num2str(irun) '.mat'], 'armpos', 'behave_change_log', 'behave_ind', 'behavior', 'dirdat', ...
                'error_correct', 'headingarm', 'laps_coverspace', ...
                'laps_singlepass', 'laps_twoarms', 'linposcat', 'linposcatnan', 'linposnorm', ...
                'pos', 'rawspikedata', 'vel', 'HP_Ripple', 'HP_SlowGamma', 'HP_Theta', ...
                'slowgamma_localspikephase', 'HPTheta_spikephase')
        end
        armpos1 = cat(1,armpos1,armpos);
        behave_change_log1 = cat(1,behave_change_log1,behave_change_log);
        behave_ind1 = cat(1,behave_ind1,behave_ind);
        behavior1 = cat(1,behavior1,behavior);
        dirdat1 = cat(1,dirdat1,dirdat);
        error_correct1 = cat(1,error_correct1,error_correct);
        headingarm1 = cat(1,headingarm1,headingarm);
        linposcat1 = cat(1,linposcat1,linposcat);
        linposcatnan1 = cat(1,linposcatnan1,linposcatnan);
        linposnorm1 = cat(1,linposnorm1,linposnorm);
        pos1 = cat(1,pos1,pos);
        rawspikedata1 = cat(1,rawspikedata1,rawspikedata);
        vel1 = cat(1,vel1,vel);
        HP_Ripple1 = cat(1,HP_Ripple1,HP_Ripple);
        HP_SlowGamma1 = cat(1,HP_SlowGamma1,HP_SlowGamma);
        HP_Theta1 = cat(1,HP_Theta1,HP_Theta);
        slowgamma_localspikephase1 = cat(1,slowgamma_localspikephase1,slowgamma_localspikephase);
%         HPTheta_spikephase1 = cat(1,HPTheta_spikephase1,HPTheta_spikephase);
        
        
        laps_coverspace1 = cat(1,laps_coverspace1,laps_coverspace+lapnum(1));
        lapnum(1) = max(laps_coverspace);
        laps_singlepass1 = cat(1,laps_singlepass1,laps_singlepass+lapnum(2));
        lapnum(2) = max(laps_singlepass);
        laps_twoarms1 = cat(1,laps_twoarms1,laps_twoarms+lapnum(3));
        lapnum(3) = max(laps_twoarms);
        
    end 
    
    armpos = armpos1;
    behave_change_log = behave_change_log1;
    behave_ind = behave_ind1;
    behavior = behavior1;
    dirdat = dirdat1;
    error_correct = error_correct1;
    headingarm = headingarm1;
    laps_coverspace = laps_coverspace1;
    laps_singlepass = laps_singlepass1;
    laps_twoarms = laps_twoarms1;
    linposcat = linposcat1;
    linposcatnan = linposcatnan1;
    linposnorm = linposnorm1;
    pos = pos1;
    rawspikedata = rawspikedata1;
    vel = vel1;
    HP_Ripple = HP_Ripple1;
    HP_SlowGamma = HP_SlowGamma1;
    HP_Theta = HP_Theta1;
    slowgamma_localspikephase = slowgamma_localspikephase1;
%     HPTheta_spikephase = HPTheta_spikephase1;
    
    [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params);

    slowgamma_localspikephase(:,2) = spikedata(:,2);
    HPTheta_spikephase = make_HPtheta_spikephase2(spikedata,HP_Theta,pos);
    
    save(['D:\XY_matdata\AllDays\' dirname2 '.mat'], 'Date', 'Rat_Name', 'dirname', 'Track_Type', ...
                'armpos', 'behave_change_log', 'behave_ind', 'behavior', 'cm_conv', 'dirdat', ...
                'error_correct', 'headingarm', 'laps_coverspace', ...
                'laps_singlepass', 'laps_twoarms', 'linposcat', 'linposcatnan', 'linposnorm', ...
                'params', 'pos', 'rawspikedata', 'vel', ...
                'HP_Ripple', 'HP_SlowGamma', 'HP_Theta', 'TT', 'armposindex', ...
                'slowgamma_localspikephase', 'HPTheta_spikephase','hp_cells','other_cells','spikedata','hpinterneurons')          
    
end