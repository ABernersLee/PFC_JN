%this is a problem script because I break up each arm into bins but then I
%can use any of the bins to predict other bins so theres a chance the
%decoder is getting information about the arm/trial/choice the animal is on in the test
%data in the training data, but real predictiveness should be more than
%predicting what arm you are going to given you know that a second ago the
%FR meant that you were going to a given arm. That is too weak. 
%Changed this to a stricter version where I still use binsize but only test on the one lap.
function [pPFCvShuff,pPFCvHP,pHPvShuff,Armind,PFCreal,PFCshuff,HPreal,HPshuff,CoeffsPFC,CoeffsHP] = prospective_coding_test(thisdir,toplot,label,num_shuffles,num_hp_shuffles)
binsize = 6;
load(thisdir,'other_cells','spikedata','vel','hp_cells','hpinterneurons','behave_change_log','laps_singlepass','linposcat','armpos')
hpcells = setdiff(hp_cells,hpinterneurons);
velcutoff = 5;
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
   
   lappos = linposcat(lapind);
   subtracter = min(lappos);
   lappos2 = lappos-subtracter+.0001;
   [h,c,~] = histcounts(lappos2,binsize);
   if sum(h(4:end)==0)>0
       exclude = lappos2==max(lappos2);
       lapind(exclude) = [];
       lappos2(exclude) = [];
       [h,c,~] = histcounts(lappos2,binsize);
   end
   
   
   nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
   thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
   
   sP = spksPFC(ismember(spksPFC(:,3), lapind),1:3);
   sH = spksHP(ismember(spksHP(:,3), lapind),1:3);
   frP = zeros(binsize,length(other_cells)); 
   frH = zeros(binsize,length(hpcells)); 
     
   for icell = 1:length(other_cells)
       if sum(sP(:,2)==other_cells(icell))>0
          frP(:,icell) = histcounts(linposcat(sP(sP(:,2)==other_cells(icell),3))-subtracter,c)./h;
       end
   end
   
   
   for icell = 1:length(hpcells)
       if sum(sH(:,2)==hpcells(icell))>0          
          frH(:,icell) = histcounts(linposcat(sH(sH(:,2)==hpcells(icell),3))-subtracter,c)./h;
       end
   end
   frP = frP(:)'; frH = frH(:)';
   FRpfc = cat(1,FRpfc,frP./binsize);
   FRhp = cat(1,FRhp,frH./binsize);

%   frP = zeros(1,length(other_cells)); 
%    frH = zeros(1,length(hpcells)); 
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
%    FRpfc = cat(1,FRpfc,frP);
%    FRhp = cat(1,FRhp,frH);
   
   traj = cat(1,traj,[nextarm thisarm ilap]);
%    traj = cat(1,traj,ones(max(lapindb),1)*[nextarm thisarm]);
%    lap_ind = cat(1,lap_ind,ones(max(lapindb),1)*ilap);
end

            

% n = floor(max(lap_ind)/2);
% m = max(lap_ind)+1-n;
exclude = sum(~isnan(FRpfc),2)==0 | sum(~isnan(FRhp),2)==0;
FRpfc(exclude,:) = [];
FRhp(exclude,:) = [];
traj(exclude,:) = [];
% lap_ind(exclude,:) = [];

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
    LabelLap = traj(touse,3);
    PFCarm = FRpfc(touse,:);
%     lap_ind2 = lap_ind(touse,1);
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
    
%     num_folds = length(Labels);
    num_folds = length(unique(LabelLap));
    
    iterX = NaN(num_folds,1); 
%     iterSX = NaN(num_folds,num_shuffles*round(num_shuffles/num_hp_shuffles));
    iterSX = NaN(num_folds,num_shuffles);
    iterXH = NaN(num_folds,num_hp_shuffles); 
    iterSXH = NaN(num_folds,round(num_shuffles/num_hp_shuffles),num_shuffles);
    
    hs(iarm,1) = length(Labels);
    h = hist(Labels,setdiff(1:3,iarm));
    hs(iarm,2) = min(h);    
    randss = randi(num_folds-1,num_folds-1,num_shuffles*round(num_shuffles/num_hp_shuffles));
    CoeffsPFC1 = NaN(size(FRpfc,2),num_folds);
    CoeffsHP1 = NaN(size(FRhp,2),num_folds,num_hp_shuffles);
    if size(PFCarm,2)<2 || hs(iarm,2)<2 %%|| length(Labels)<20 % length(Labels)<size(PFCarm,2) ||
            continue
    else
%          indices = 1:length(Labels);
         indices = unique(LabelLap);
        for i = 1:num_folds
%                 test = (indices == i); train = ~test;
                test = (LabelLap == indices(i)); train = ~test;


%                 ind2 = find(train);  
                trainlabels = Labels(train);
                 ns = histc(trainlabels,unique(trainlabels));
                        
                 if sum(ns)<171                    
                   numu =  factorial(sum(ns))./(factorial(ns(1))*factorial(ns(2)));
                 else numu = numshuff+1;
                 end

                 if numu<num_shuffles
                    disp([numu iarm i])
                    continue                
                 end
                    
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
                    
                 n = numel(trainlabels);
                 if n >= 33
%                             disp('Maximum variable size allowed by the program is exceeded.')
                    perm_list = NaN(num_shuffles,n);
                    for ishuff = 1:num_shuffles
                       perm_list(ishuff,:) = trainlabels(randperm(n));
                    end
                    numwhile = 0;
                    while (size(unique(perm_list,'rows'),1)<num_shuffles) && numwhile<5000
                         for ishuff = 1:num_shuffles
                           perm_list(ishuff,:) = trainlabels(randperm(n));
                         end
                         numwhile = numwhile+1;
                    end
                    if numwhile>=5000                                
                        perm_list = unique_perms(trainlabels');
                        perm_list = perm_list(randperm(num_shuffles),:);
                    end
                else
                    perm_list = unique_perms(trainlabels');
                    perm_list = perm_list(randperm(num_shuffles),:);

                end
                 for j = 1:num_shuffles %*round(num_shuffles/num_hp_shuffles)
%                     ss = randss(:,j);                        
%                     guessS = classify(testdata,traindata,Labels(ind2(ss)),'linear');   
%                     Mdl = fitcdiscr(traindata,Labels(ind2(ss)),'discrimType','pseudoLinear');
                    Mdl = fitcdiscr(traindata,perm_list(j,:),'discrimType','pseudoLinear');
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
%                         Mdl = fitcdiscr(traindataHP,Labels(ind2(ss)),'discrimType','pseudoLinear');
                        Mdl = fitcdiscr(traindataHP,perm_list(j,:),'discrimType','pseudoLinear');                        
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



dirname = thisdir(1:end-4);
if toplot
    figure; hold on;
%     histogram(allaccS(:),20); histogram(allacc(:),20);  legend('Shuffle','Data')
    histogram(nanmean(PFCshuff),'FaceColor','k');
    yl= get(gca,'ylim');
    plot([nanmean(PFCreal) nanmean(PFCreal)],yl,'r-','LineWidth',3)
    title([dirname ' p = ' num2str(round(pPFCvShuff,3,'significant'))])
    helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\ProspectiveFRDecoding\New_AllClassifyTest_' dirname])
end