function linear_classifier_of_replay_withPFC(thisdir,label,type,Labels,savelab,downsample)
if downsample
    savelab = [savelab '_ds'];
end
%type 1 = linear discriminant analysis
%type 2 = Multinomial logistic regression
%type 3 = psuedolinear discriminant analysis
load(thisdir,[label '_PFCreplayspikes_binned'],...
    [label  '_replay_replayarm'],[label  '_replay_singlebothjoint'],'other_cells',[label '_Cand_sig_modu_include'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['PFC_binned = ' label '_PFCreplayspikes_binned;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
clear([label '_Cand_sig_modu_include'],[label '_PFCreplayspikes_binned'],[label '_replay_singlebothjoint'],[label '_replay_replayarm'])

pfc = other_cells(Cand_sig_modu_include(:,3)==1); clear other_cells;
PFC_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];


Labels= Labels(singlebothjoint~=3);
if downsample
    h = hist(Labels,1:3);
    gto = min(h);
    loarm = find(gto==h);
    armsout = setdiff(1:3,loarm);
    for iarm = 1:2
       indlab = find(Labels==armsout(iarm));
       randx = randperm(length(indlab));
       randxx = indlab(randx);
       Labels(randxx(1:length(indlab)-gto)) = NaN;
    end
    exevent = isnan(Labels);
    Labels(exevent) = [];
else
    exevent = false(size(Labels));
end


binsize = .02;
window = [-2 2];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
num_shuffles = 500;
mind = [0 .2]; bind = [-.5 -.1];
modind = ind>=mind(1) & ind<mind(2);
baseind = ind>=bind(1) & ind<bind(2);

PFC_binned = PFC_binned(:,:,singlebothjoint~=3);
bindat = PFC_binned; 


eacharm = zeros(3,3);

if ~isempty(pfc) 
    tic
    indices = 1:length(Labels);
    num_folds = length(Labels);
    
    randss = randi(num_folds-1,num_folds-1,num_shuffles);
    while any(sum(randss~=([1:length(Labels)-1]'))==0)
        randss = randi(num_folds-1,num_folds-1,num_shuffles);
    end
    iterX = NaN(num_folds,1); iterSX = NaN(3,num_folds,num_shuffles); iterXarm = NaN(3,num_folds);
    bindat1 = squeeze((sum(bindat(:,modind,:),2)/sum(modind))-(sum(bindat(:,baseind,:),2))/sum(baseind));
    bindat1(:,exevent) = []; %added 9/8/20
    
    excl = sum(unique(bindat1','rows')~=0)<4; %added 5/20/19
    bindat2 = zscore(bindat1(~excl,:)');
    
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        traindata = bindat2(train,:);
        testdat = bindat2(test,:);
        excl = std(traindata,[],1)<10^-10;
        traindata(:,excl) = []; testdat(:,excl) = [];
        if type == 1
            guess = classify(testdat,traindata,Labels(train),'diaglinear');  %changed from linear only for day 11, on 10/25, revisit
        elseif type == 2
            B = mnrfit(traindata',Labels(train));
            f = mnrval(B,testdat');
            [~,guess] = max(f,[],2);
        elseif type == 3
            Mdl = fitcdiscr(traindata',Labels(train),'discrimType','pseudoLinear');
            guess = predict(Mdl,testdat');
        end
        real = Labels(test);
        iterX(i) = real==guess;
        eacharm(real,guess) = eacharm(real,guess)+1; 
        iterXarm(real,i) = real==guess;
        
        ind2 = find(train);
        for j = 1:num_shuffles  
            randoth = randss(:,~ismember(1:size(randss,2),j));
            while any(sum(randoth~=randss(:,j))==0)
                randss(:,j) = randi(num_folds-1,num_folds-1,1);
            end
            ss = randss(:,j);
            if type == 1
                guessS = classify(testdat,traindata,Labels(ind2(ss)),'diaglinear');    %changed from linear only for day 11, on 10/25, revisit
            elseif type == 2
                B = mnrfit(traindata',Labels(ind2(ss)));
                f = mnrval(B,testdat');
                [~,guessS] = max(f,[],2);
            elseif type == 3
                Mdl = fitcdiscr(traindata',Labels(ind2(ss)),'discrimType','pseudoLinear');
                guessS = predict(Mdl,testdat');
            end
            iterSX(real,i,j) = real==guessS;            
        end
        
        if rem(i,500)==0
            t = toc;
           disp([num2str(i) ' Out of ' num2str(num_folds) ' Folds Finished in ' num2str(t/60) ' Minutes']) 
           tic
        end
    end
    
    r = sum(iterX)/num_folds;    
    s = sum(squeeze(nansum(iterSX,1)),1)./num_folds;
    p = (sum(s>=r)+1)/(length(s)+1);
    a = nansum(iterXarm,2)./sum(~isnan(iterSX(:,:,1)),2);
    pa = (sum((squeeze(nansum(iterSX,2))./sum(~isnan(iterSX(:,:,1)),2))>=a,2)+1)/(length(s)+1);
elseif isempty(pfc) 
   p = NaN; r = NaN; s = NaN; pa = NaN; a = NaN;
end


figure; imagesc((eacharm-min(eacharm,[],2))./(max(eacharm,[],2)-min(eacharm,[],2))); hold on; axis xy; ylabel('Real'); xlabel('Guess'); title(['p = ' num2str(round(p,2,'significant')) ', armPs = ' num2str(round(pa',2,'significant'))])
helper_saveandclosefig(['F:\XY_matdata\Figures\ForPaperJN\AllCells_AllArmEvents\Classify\Sessions\ByArm_' thisdir(1:end-4) '_' savelab])

figure; imagesc(eacharm./sum(eacharm,2)); hold on; axis xy; ylabel('Real'); xlabel('Guess'); title(['p = ' num2str(round(p,2,'significant')) ', armPs = ' num2str(round(pa',2,'significant'))])
helper_saveandclosefig(['F:\XY_matdata\Figures\ForPaperJN\AllCells_AllArmEvents\Classify\Sessions\ByArm_prop_' thisdir(1:end-4) '_' savelab])

figure; imagesc(eacharm); hold on; axis xy; ylabel('Real'); xlabel('Guess'); title(['p = ' num2str(round(p,2,'significant')) ', armPs = ' num2str(round(pa',2,'significant'))])
helper_saveandclosefig(['F:\XY_matdata\Figures\ForPaperJN\AllCells_AllArmEvents\Classify\Sessions\ByArm_Raw_' thisdir(1:end-4) '_' savelab])

eval([[label '_' savelab '_Classifier_type' num2str(type) '_p'] '= p;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_eacharm'] '= eacharm;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_r'] '= r;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_s'] '= s;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_a'] '= a;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_pa'] '= pa;'])

save(thisdir,[label '_' savelab '_Classifier_type' num2str(type) '_p'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_r'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_s'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_a'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_pa'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_eacharm'],'-append')


