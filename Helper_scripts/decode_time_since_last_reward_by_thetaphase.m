function [accdiff1,accdiff2,accall,accSall] = decode_time_since_last_reward_by_thetaphase(thisdir,numshuff,indinds,indinds2)
load(thisdir,'PFCthetaspikes_binned');
Th = add_time_since_last_reward(thisdir);
% ThInd = ~isnan(Th(:,10));
% Th(ThInd,10) = Th(ThInd,10)>nanmedian(Th(ThInd,10));

binsize = .02;
window = [-2 2];
ind = [window(1):binsize:window(2)];

% modind = [ind>=-.03 & ind<.03]; %this too, STARK contrast in the difference in decoding ability
% modind2 = [ind>=-.09 & ind<-.03];

% modind = [ind>=-.02 & ind<.04]; %these do better for session 1 at least, sig vs nonsig
% modind2 = [ind>=-.08 & ind<-.02];

% 0 .06
modind = [ind>=indinds(1) & ind<(indinds(2))];
% modind2 = [ind>=-.06 & ind<0];
modind2 = [ind>=indinds2(1) & ind<(indinds2(2))];

% modind = [ind>=0 & ind<.05];
% modind2 = [ind>=-.05 & ind<0];

% modind2 = true(size(ind));
% numshuff = 50;


accall = [];
accSall = [];
for iarm = 1:3
%     for ilap = 1:max(Th(:,6))
%         touse = Th(:,3)==iarm & Th(:,6)==ilap;
    touse = Th(:,3)==iarm;
    if sum(touse)==0; continue; end
    labels = Th(touse,10);
    ThL = Th(touse,6);
    laplaps = unique(ThL(~isnan(ThL)));
    numlaps = length(unique(ThL(~isnan(ThL))));
    
%     dat = squeeze(nansum(PFCthetaspikes_binned(:,modind,touse),2)); %./(sum(modind)*binsize));
    dat2 = squeeze(nansum(PFCthetaspikes_binned(:,modind2,touse),2)); %./(sum(modind2)*binsize));
%     dat(sum(dat,2)<5,:) = [];
    dat2(mean(dat2,2)<.015,:) = [];
    dat = nanzscore(dat2,[],2);
    if length(unique(labels))==1; continue; end

    
    for ii = 1:numlaps
       testdat = dat(:,ThL==laplaps(ii));
       train = find(ThL~=laplaps(ii) & ~isnan(ThL));
       traindat = dat(:,train);       
       guess = classify(testdat',traindat',labels(train),'diaglinear');  
       real = labels(ThL==laplaps(ii));
       accall = [accall;guess==real];
       accuracy_shuff = NaN(length(real),numshuff);
       for ishuff = 1:numshuff
          train2 = train(randperm(length(train)));
          while sum(train2~=train)==0
              train2 = train(randperm(length(train)));
          end
          guess2 = classify(testdat',traindat',labels(train2),'diaglinear');         
          accuracy_shuff(:,ishuff) = guess2==real;
       end
        accSall = [accSall;accuracy_shuff];       
    end
%     accall = [accall;accuracy];
%     accSall = [accSall;accuracy_shuff];
%     end
end
p1 = (sum(nanmean(accSall)>=mean(accall))+1)/(numshuff+1);
disp(['local: p = ' num2str(p1)])
figure; 
histogram(nanmean(accSall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(accall) mean(accall)],yl,'r-','LineWidth',3)
title([thisdir ' local p = ' num2str(p1)])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Classify\Sessions\' thisdir(1) thisdir(3:end-4) '_Time_Since_Last_Reward__Local_' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'ms_numshuff_' num2str(numshuff)])
    
accdiff2 = mean(accall)-nanmean(accSall);

accall = [];
accSall = [];
for iarm = 1:3
%     for ilap = 1:max(Th(:,6))
%         touse = Th(:,3)==iarm & Th(:,6)==ilap;
    touse = Th(:,3)==iarm;
    if sum(touse)==0; continue; end
    labels = Th(touse,10);
    ThL = Th(touse,6);
    laplaps = unique(ThL(~isnan(ThL)));
    numlaps = length(unique(ThL(~isnan(ThL))));
    
    dat = squeeze(nansum(PFCthetaspikes_binned(:,modind,touse),2)); %./(sum(modind)*binsize));
%     dat2 = squeeze(nansum(PFCthetaspikes_binned(:,modind2,touse),2)); %./(sum(modind2)*binsize));
    dat(mean(dat,2)<.015,:) = [];
%     dat2(sum(dat2,2)<5,:) = [];
    dat = nanzscore(dat,[],2);
    if length(unique(labels))==1; continue; end

    
    for ii = 1:numlaps
       testdat = dat(:,ThL==laplaps(ii));
       train = find(ThL~=laplaps(ii) & ~isnan(ThL));
       traindat = dat(:,train);       
       guess = classify(testdat',traindat',labels(train),'diaglinear');  
       real = labels(ThL==laplaps(ii));
       accall = [accall;guess==real];
       accuracy_shuff = NaN(length(real),numshuff);
       for ishuff = 1:numshuff
          train2 = train(randperm(length(train)));
          while sum(train2~=train)==0
              train2 = train(randperm(length(train)));
          end
          guess2 = classify(testdat',traindat',labels(train2),'diaglinear');         
          accuracy_shuff(:,ishuff) = guess2==real;
       end
        accSall = [accSall;accuracy_shuff];       
    end
%     accall = [accall;accuracy];
%     accSall = [accSall;accuracy_shuff];
%     end
end
p1 = (sum(nanmean(accSall)>=mean(accall))+1)/(numshuff+1);
disp(['sweep: p = ' num2str(p1)])
figure; histogram(nanmean(accSall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(accall) mean(accall)],yl,'r-','LineWidth',3)
title([thisdir ' sweep foward, p = ' num2str(p1)])
helper_saveandclosefig(['E:\XY_matdata\Figures\ForPaper\Classify\Sessions\' thisdir(1) thisdir(3:end-4) '_Time_Since_Last_Reward_NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff)])
accdiff1 = mean(accall)-nanmean(accSall);


