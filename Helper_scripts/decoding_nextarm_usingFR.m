function decoding_nextarm_usingFR(thisdir,num_shuffles,num_hp_shuffles,velcutoff)

%formerly prospective_coding_test
load(thisdir,'other_cells','spikedata','vel','hp_cells','hpinterneurons','behave_change_log','laps_singlepass','armpos')
hpcells = setdiff(hp_cells,hpinterneurons);
% get out lap data for both PFC and HP

spksPFC = spikedata(ismember(spikedata(:,2),other_cells),:);
spksHP = spikedata(ismember(spikedata(:,2),hpcells),:);
traj = []; FRpfc = []; lap_ind = []; FRhp = [];
for ilap = 1:max(laps_singlepass)
    
    %real
   leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
   lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle
   
   
       %testing
%    leave = find(behave_change_log(:,4) & laps_singlepass == ilap,1,'first'); %4 leave lick
%    if isempty(leave)
%        leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %4 leave lick
%    end
%    lapind = find(laps_singlepass==ilap,1,'first'):leave;
   
   lapind(vel(lapind)<velcutoff) = []; %changed 10/22 at 6:40pm
   
   binz = binsize/(1/60);
   lapindb = ceil((1:length(lapind))/binz);
%     lapindb = 1;
   
   nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
   thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
   
   sP = spksPFC(ismember(spksPFC(:,3), lapind),1:2);
   sH = spksHP(ismember(spksHP(:,3), lapind),1:2);
   frP = zeros(max(lapindb),length(other_cells)); 
   frH = zeros(max(lapindb),length(hpcells)); 
   
   [~,~,i] = histcounts(sP(:,1),pos(lapind(1),1):binsize:pos(lapind(end),1));
   for icell = 1:length(other_cells)
       if sum(sP(:,2)==other_cells(icell))>0
          frP(:,icell) = histc(i(sP(:,2)==other_cells(icell)),1:max(lapindb));
       end
   end
   
   [~,~,i] = histcounts(sH(:,1),pos(lapind(1),1):binsize:pos(lapind(end),1));
   for icell = 1:length(hpcells)
       if sum(sH(:,2)==hpcells(icell))>0
          frH(:,icell) = histc(i(sH(:,2)==hpcells(icell)),1:max(lapindb));
       end
   end
   
   FRpfc = cat(1,FRpfc,frP./binsize);
   FRhp = cat(1,FRhp,frH./binsize);

   
%    for icell = 1:length(other_cells)
%        if sum(sP(:,2)==other_cells(icell))>0
%           frP(:,icell) = sum(sP(:,2)==other_cells(icell));
%        end
%    end
%    
%    
%    for icell = 1:length(hpcells)
%        if sum(sH(:,2)==hpcells(icell))>0
%           frH(:,icell) = sum(sH(:,2)==hpcells(icell));
%        end
%    end
%    
%    FRpfc = cat(1,FRpfc,frP./length(lapind));
%    FRhp = cat(1,FRhp,frH./length(lapind));
   
   traj = cat(1,traj,ones(max(lapindb),1)*[nextarm thisarm]);
   lap_ind = cat(1,lap_ind,ones(max(lapindb),1)*ilap);
end

            

% n = floor(max(lap_ind)/2);
% m = max(lap_ind)+1-n;
exclude = sum(~isnan(FRpfc),2)==0 | sum(~isnan(FRhp),2)==0;
FRpfc(exclude,:) = [];
FRhp(exclude,:) = [];
traj(exclude,:) = [];
lap_ind(exclude,:) = [];

PFCreal =[];
PFCshuff = [];
HPreal =[];
HPshuff = [];
Armind = [];
CoeffsPFC = NaN(size(FRpfc,2),3); 
CoeffsHP = NaN(size(FRhp,2),3);
% CoeffsPFC12 = NaN(size(FRpfc,2),3,2); 
% CoeffsHP12 = NaN(size(FRhp,2),3,2);


for iarm = 1:3     
        
    touse = traj(:,2)==iarm;
    
    % touse = traj(:,2)==iarm & lap_ind<=n; %     touse = traj(:,2)==iarm & lap_ind>=m;
    Labels = traj(touse,1);    
    PFCarm = FRpfc(touse,:);
%     PFCind = sum(PFCarm~=0)>5; 
    PFCind = sum(PFCarm)~=0;
    PFCarm(:,~PFCind) = [];
    PFCarm = nanzscore(PFCarm);
    
    HParm = FRhp(touse,:);
    HPind = find(sum(HParm)~=0);     
    HParm(:,sum(HParm)==0) = []; 
    HParm = nanzscore(HParm);
%     HPind = find(sum(HParm~=0)>5);
%     HParm(:,sum(HParm~=0)<=5) = []; 
    
    num_folds = length(Labels);
    
    iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles*round(num_shuffles/num_hp_shuffles));
    iterXH = NaN(num_folds,num_hp_shuffles); iterSXH = NaN(num_folds,round(num_shuffles/num_hp_shuffles),num_shuffles);
    
    hs(iarm,1) = length(Labels);
    h = hist(Labels,setdiff(1:3,iarm));
    hs(iarm,2) = min(h);    
    randss = randi(num_folds-1,num_folds-1,num_shuffles*round(num_shuffles/num_hp_shuffles));
    CoeffsPFC1 = NaN(size(FRpfc,2),num_folds);
    CoeffsHP1 = NaN(size(FRhp,2),num_folds,num_hp_shuffles);
    if size(PFCarm,2)<2 || hs(iarm,2)<2 %%|| length(Labels)<20 % length(Labels)<size(PFCarm,2) ||
            continue
    else
         indices = 1:length(Labels);            
        for i = 1:num_folds
                test = (indices == i); train = ~test;

                %PFC real data and shuffle
                traindata = PFCarm(train,:);
                testdata = PFCarm(test,:);
                            
                
                Mdl = fitcdiscr(traindata,Labels(train),'discrimType','pseudoLinear');
                CoeffsPFC1(PFCind,i) = Mdl.Coeffs(1,2).Linear;
                guess = predict(Mdl,testdata);
                                
%                 [guess,~,~,~,d] = classify(testdata,traindata,Labels(train),'linear');     
%                 CoeffsPFC1(PFCind,i) = d(2,1).linear;
                                
                real = Labels(test);
                iterX(i) = real==guess;
                ind2 = find(train);   
                 for j = 1:num_shuffles*round(num_shuffles/num_hp_shuffles)
                    ss = randss(:,j);            
%                     guessS = classify(testdata,traindata,Labels(ind2(ss)),'linear');   
                    Mdl = fitcdiscr(traindata,Labels(ind2(ss)),'discrimType','pseudoLinear');
                    guessS = predict(Mdl,testdata);
                    iterSX(i,j) = real==guessS;
                 end

                 %HP downsampled distribution and shuffle
                for ihp = 1:num_hp_shuffles
                    randhp = randperm(size(HParm,2),size(PFCarm,2));
                    thisHP = HParm(:,randhp);
                    
%                     thisHP = zscore(thisHP);
                    
                    traindataHP = thisHP(train,:);
                    testdataHP = thisHP(test,:);
%                     [guess,~,~,~,d] = classify(testdataHP,traindataHP,Labels(train),'linear');     
%                     CoeffsHP1(HPind(randhp),i,ihp) = d(2,1).linear;
                    Mdl = fitcdiscr(traindataHP,Labels(train),'discrimType','pseudoLinear');
                    CoeffsHP1(HPind(randhp),i,ihp) = Mdl.Coeffs(1,2).Linear;
                    guess = predict(Mdl,testdataHP);                    
                    iterXH(i,ihp) = real==guess;
                    for j = 1:round(num_shuffles/num_hp_shuffles)
                        ss = randss(:,j);            
%                         guessS = classify(testdataHP,traindataHP,Labels(ind2(ss)),'linear');   
                        Mdl = fitcdiscr(traindataHP,Labels(ind2(ss)),'discrimType','pseudoLinear');
                        guessS = predict(Mdl,testdataHP);
                        iterSXH(i,j,ihp) = real==guessS;                                                
                    end
                end  
                
        end        
    end
    CoeffsPFC(:,iarm) = nanmean(CoeffsPFC1,2);
    CoeffsHP(:,iarm) = nanmean(nanmean(CoeffsHP1,2),3);
    Armind = [Armind; iarm*ones(size(iterX))];
    PFCreal = [PFCreal; iterX];
    PFCshuff = [PFCshuff;iterSX];
    HPreal = [HPreal; iterXH];
    HPshuff = [HPshuff;reshape(iterSXH,size(iterSXH,1),size(iterSXH,2)*size(iterSXH,3))];
    disp(['Done with Arm ' num2str(iarm)])
end
pPFCvShuff = (sum(nanmean(PFCshuff,1)>nanmean(PFCreal))+1)/(sum(~isnan(nanmean(PFCshuff,1)))+1);
pPFCvHP = (sum(nanmean(HPreal,1)>nanmean(PFCreal))+1)/(sum(~isnan(nanmean(HPreal,1)))+1);
pHPvShuff = (sum(nanmean(HPshuff,1)>nanmean(nanmean(HPreal)))+1)/(sum(~isnan(nanmean(HPshuff,1))+1));


if sum(~isnan(nanmean(PFCshuff,1)))==0
    pPFCvShuff = NaN;
    pPFCvHP = NaN;
    pHPvShuff = NaN;
end



ps = [pPFCvShuff,pPFCvHP,pHPvShuff];
decodenextarmFR_PFCreal = PFCreal;
decodenextarmFR_PFCshuff = PFCshuff;
decodenextarmFR_HPreal = HPreal;
decodenextarmFR_HPshuff = HPshuff;
decodenextarmFR_CoeffsHP = CoeffsHP;
decodenextarmFR_CoeffsPFC = CoeffsPFC;
decodenextarmFR_ps = ps;

save(thisdir,'decodenextarmFR_PFCreal','decodenextarmFR_ps','decodenextarmFR_PFCshuff','decodenextarmFR_HPreal'...
    ,'decodenextarmFR_HPshuff','decodenextarmFR_CoeffsHP','decodenextarmFR_CoeffsPFC','-append')

disp('Done with decoding_nextarm_usingFR')