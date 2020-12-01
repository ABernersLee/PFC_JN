function [accdiff1,accdiff2,accdiff1v2,accdiff2v2,accall,accSall,accallL,accSallL,p1,lapaccL,lapaccA] = decode_choice_by_thetaphase(thisdir,numshuff,indinds,indinds2,toplot,savefolder)
load(thisdir,'times_armon_thetaof_headingarm_lap_thetahalf_all','PFCthetaspikes_binned');
Th = times_armon_thetaof_headingarm_lap_thetahalf_all;
clear times_armon_thetaof_headingarm_lap_thetahalf_all

% indinds = indinds-(theta_seqeuence_zero/360)*.125;
% indinds2 = indinds2-(theta_seqeuence_zero/360)*.125;
if ~isfolder([savefolder '\Classify\Sessions\'])
    mkdir([savefolder '\Classify\Sessions\'])
end
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

lapacc = [];
accall = [];
accSall = [];
tousearm1 = true(3,1);
for iarm = 1:3
%     for ilap = 1:max(Th(:,6))
%         touse = Th(:,3)==iarm & Th(:,6)==ilap;
    touse = Th(:,3)==iarm;
%     touse = Th(:,3)==iarm & Th(:,end)>.34;
    if sum(touse)==0; continue; end
    labels = Th(touse,5);
    ThL = Th(touse,6);
    laplaps = unique(ThL(~isnan(ThL)));
    numlaps = length(unique(ThL(~isnan(ThL))));
    
    if numlaps<6; continue; end % testing
    
%     dat = squeeze(nansum(PFCthetaspikes_binned(:,modind,touse),2)); %./(sum(modind)*binsize));
    dat2 = squeeze(nansum(PFCthetaspikes_binned(:,modind2,touse),2)); %./(sum(modind2)*binsize));
%     dat2(sum(dat2,2)<5,:) = [];
    dat2(mean(dat2,2)<.1,:) = [];
%     dat2(mean(dat2,2)<.015,:) = [];
    dat = nanzscore(dat2,[],2);
    if length(unique(labels))==1 || isempty(dat)
        tousearm1(iarm,1) = false;
        continue; 
    end

    
    for ii = 1:numlaps
       testdat = dat(:,ThL==laplaps(ii));
       train = find(ThL~=laplaps(ii) & ~isnan(ThL));
       real = labels(ThL==laplaps(ii));
       h = histc(labels(train),unique(labels(train))); 
       if nchoosek(sum(h),min(h))<numshuff
           disp(['Skipping, only ' num2str(nchoosek(sum(h),min(h))) ' perms'])
           continue
       end
       traindat = dat(:,train);       
       excl = std(traindat,[],2)<10^-10;
       traindat(excl,:) = []; testdat(excl,:) = [];
       if isempty(traindat)
           continue
       end
       guess = classify(testdat',traindat',labels(train),'diaglinear');  
       
%        indall = [indall; 
       accall = [accall;guess==real];
       accuracy_shuff = NaN(length(real),numshuff);
       for ishuff = 1:numshuff
          train2 = train(randperm(length(train)));
          while sum(labels(train2)~=labels(train))==0
              train2 = train(randperm(length(train)));
          end
          if ishuff == 1; shuffsave = labels(train2); 
          else
              while any(sum(shuffsave~=labels(train2))==0)
                  train2 = train(randperm(length(train))); 
              end
              shuffsave = [shuffsave labels(train2)]; 
          end
          guess2 = classify(testdat',traindat',labels(train2),'diaglinear');         
          accuracy_shuff(:,ishuff) = guess2==real;
       end
       clear shuffsave
        accSall = [accSall;accuracy_shuff];  
        lapacc = [lapacc; mean(guess==real)-mean(mean(accuracy_shuff))];        
    end
%     accall = [accall;accuracy];
%     accSall = [accSall;accuracy_shuff];
%     end
end
if ~isempty(accall) && toplot
p1 = (sum(nanmean(accSall)>=mean(accall))+1)/(numshuff+1);
disp(['local: p = ' num2str(p1)])
figure; 
histogram(nanmean(accSall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(accall) mean(accall)],yl,'r-','LineWidth',3)
title([thisdir ' local p = ' num2str(p1)])
helper_saveandclosefig([savefolder '\Classify\Sessions\' thisdir(1) thisdir(3:end-4) '_Local_' num2str(indinds2(1)) 'ms_to_' num2str(indinds2(2)) 'ms_numshuff_' num2str(numshuff)])    
end
accdiff2 = mean(accall)-nanmean(accSall);
accdiff2v2 = accall-nanmean(accSall,2);
accallL = accall; accSallL = accSall;
lapaccL =lapacc;

accall = [];
accSall = [];
lapacc = [];
tousearm2 = true(3,1);
for iarm = 1:3
%     for ilap = 1:max(Th(:,6))
%         touse = Th(:,3)==iarm & Th(:,6)==ilap;
    touse = Th(:,3)==iarm;
    if sum(touse)==0; continue; end
    labels = Th(touse,5);
    ThL = Th(touse,6);
    laplaps = unique(ThL(~isnan(ThL)));
    numlaps = length(unique(ThL(~isnan(ThL))));
    
    if numlaps<6; continue; end %testing
    
    dat2 = squeeze(nansum(PFCthetaspikes_binned(:,modind,touse),2)); %./(sum(modind)*binsize));
%     dat2 = squeeze(nansum(PFCthetaspikes_binned(:,modind2,touse),2)); %./(sum(modind2)*binsize));
    dat2(mean(dat2,2)<.1,:) = [];
%     dat2(sum(dat2,2)<5,:) = [];
    dat = nanzscore(dat2,[],2);
     if length(unique(labels))==1 || isempty(dat) || tousearm1(iarm,1)==false
        tousearm2(iarm,1) = false;
        continue; 
    end

    
    for ii = 1:numlaps
       testdat = dat(:,ThL==laplaps(ii));
       train = find(ThL~=laplaps(ii) & ~isnan(ThL));
       real = labels(ThL==laplaps(ii));
       h = histc(labels(train),unique(labels(train))); 
       if nchoosek(sum(h),min(h))<numshuff
           disp(['Skipping, only ' num2str(nchoosek(sum(h),min(h))) ' perms'])
           continue
       end
       traindat = dat(:,train);    
       excl = std(traindat,[],2)<10^-10;
       traindat(excl,:) = []; testdat(excl,:) = [];
       if isempty(traindat)
           continue
       end
       guess = classify(testdat',traindat',labels(train),'diaglinear');  
       
       accall = [accall;guess==real];
       accuracy_shuff = NaN(length(real),numshuff);
       for ishuff = 1:numshuff
          train2 = train(randperm(length(train)));
          while sum(labels(train2)~=labels(train))==0
              train2 = train(randperm(length(train)));
          end
          if ishuff == 1; shuffsave = labels(train2); 
          else
              while any(sum(shuffsave~=labels(train2))==0)
                  train2 = train(randperm(length(train))); 
              end
              shuffsave = [shuffsave labels(train2)]; 
          end
          guess2 = classify(testdat',traindat',labels(train2),'diaglinear');         
          accuracy_shuff(:,ishuff) = guess2==real;
       end
        accSall = [accSall;accuracy_shuff]; 
        lapacc = [lapacc; mean(guess==real)-mean(mean(accuracy_shuff))];
    end
%     accall = [accall;accuracy];
%     accSall = [accSall;accuracy_shuff];
%     end
end
lapaccA =lapacc;

if sum(tousearm1==tousearm2)~=3
    error('Mismatch of arm index')
end

if ~isempty(accall) && toplot
p1 = (sum(nanmean(accSall)>=mean(accall))+1)/(numshuff+1);
disp(['sweep: p = ' num2str(p1)])
figure; histogram(nanmean(accSall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([mean(accall) mean(accall)],yl,'r-','LineWidth',3)
title([thisdir ' sweep foward, p = ' num2str(p1)])
helper_saveandclosefig([savefolder '\Classify\Sessions\' thisdir(1) thisdir(3:end-4) '_NonLocal_' num2str(indinds(1)) 'ms_to_' num2str(indinds(2)) 'ms_numshuff_' num2str(numshuff)])
end
if ~isempty(accall)
    p1 = (sum(nanmean(accSall)>=mean(accall))+1)/(numshuff+1);
else
    p1 = NaN;
end
accdiff1 = mean(accall)-nanmean(accSall);
accdiff1v2 = accall-nanmean(accSall,2);


