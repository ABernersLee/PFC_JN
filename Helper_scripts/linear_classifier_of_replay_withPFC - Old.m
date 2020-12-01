function linear_classifier_of_replay_withPFC(thisdir,label,type,Labels,savelab)

%type 1 = linear discriminant analysis
%type 2 = Multinomial logistic regression
load(thisdir,[label '_PFCreplayspikes_binned'],...
    [label  '_replay_replayarm'],[label  '_replay_singlebothjoint'],'other_cells',[label '_Cand_sig_modu_include'])
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['PFC_binned = ' label '_PFCreplayspikes_binned;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
clear([label '_Cand_sig_modu_include'],[label '_PFCreplayspikes_binned'],[label '_replay_singlebothjoint'],[label '_replay_replayarm'])
pfc = other_cells(Cand_sig_modu_include(:,3)==1); clear other_cells;
PFC_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];


binsize = .02;
window = [-2 2];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
num_shuffles = 250;
numz = 1;
num_folds = 2;
mind = [0 .2]; bind = [-.5 -.1];
modind = ind>=mind(1) & ind<mind(2);
baseind = ind>=bind(1) & ind<bind(2);
bindat = PFC_binned; 
if type==1 && ~isempty(pfc)     
   
    
    iter = NaN(num_shuffles,num_folds);
    iterS = NaN(num_shuffles,num_folds);
    for j = 1:num_shuffles
        indices = crossvalind('Kfold',Labels,num_folds); 

        for i = 1:num_folds
            test = (indices == i); train = ~test;
            traindata = (squeeze(sum(bindat(:,modind,train),2))./range(mind) - ...
                squeeze(sum(bindat(:,baseind,train),2))./range(bind))./ ...
                (squeeze(sum(bindat(:,modind,train),2))./range(mind) + ...
                squeeze(sum(bindat(:,baseind,train),2))./range(bind));
            traindata(isnan(traindata)) = 0;
            testdat = squeeze(sum(bindat(:,modind,test),2));
            guess = classify(testdat',traindata',Labels(train),'linear');        
            real = Labels(test);
            iter(j,i) = sum(real==guess)./length(guess);

%             for iz = 1:numz
                ind2 = find(train);
                ss = randperm(sum(train));
                guessS = classify(testdat',traindata',Labels(ind2(ss)),'linear');                
                iterS(j,i) = sum(real==guessS)./length(guess);
%             end
        end
    end
%     p3 = (sum(nanmean(iter(:))<iterS(:))+1)/(size(iterS,1)*size(iterS,2)+1);
%     p3 = (sum(nanmean(iter(:))<nanmean(iterS,2))+1)/(size(iterS,1)+1);
%     p3 = ranksum(iter(:),iterS(:),'tail','right');
    if lillietest([iterS(:);iter(:)])==1
        p3 = ranksum(iter(:),iterS(:));
    elseif lillietest([iterS(:);iter(:)])==0
        [~,p3] = ttest(iter(:),iterS(:));
    end
    if nanmedian(iter(:))<nanmedian(iter(:))
        p3 = 1;
    end
    dist = iter;
    sdist = iterS;
elseif type==1 && isempty(pfc)
   p3 = NaN; dist = NaN; sdist = NaN;
end
% figure; hold on; 
% errorbar(1,nanmean(iter(:)),std(mean(iter,2))./sqrt(size(iter,1)),'r'); 
% errorbar(1,nanmean(iterS(:)),std(mean(iterS,2))./sqrt(size(iterS,1)),'k')

if type == 2 && ~isempty(pfc)

    
%     X= NaN(length(pfc),sum(singlebothjoint~=3));
%     for icell = 1:length(pfc)
        X =  squeeze(sum(bindat(:,modind,:),2))./range(mind) - ...
                squeeze(sum(bindat(:,baseind,:),2))./range(bind);
%     end
    X = X';
    
    
    Fit = NaN(num_folds,num_shuffles);
    FitS = NaN(num_folds,num_shuffles);
    for j = 1:num_shuffles
        indices = crossvalind('Kfold',Labels,num_folds); 

        for i = 1:num_folds
            test = (indices == i); train = ~test;
            LL = Labels(train);
            B = mnrfit(X(train,:),LL);
            f = mnrval(B,X(test,:));
            [~,guess] = max(f,[],2);
            Fit(i,j) = sum(guess==Labels(test))./length(f);
            

            B = mnrfit(X(train,:),LL(randperm(length(LL))));
            f = mnrval(B,X(test,:));
            [~,guessS] = max(f,[],2);
            FitS(i,j) = sum(guessS==Labels(test))./length(f);
        end
    end
    p3 = signrank(Fit(:),FitS(:),'tail','right');
    dist = Fit;
    sdist = FitS;

elseif type ==2 && isempty(pfc) 
       p3 = NaN; dist = NaN; sdist = NaN;
end

eval([[label '_' savelab '_Classifier_type' num2str(type) '_p'] '= p3;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_dist'] '= dist;'])
eval([[label '_' savelab '_Classifier_type' num2str(type) '_sdist'] '= sdist;'])
save(thisdir,[label '_' savelab '_Classifier_type' num2str(type) '_p'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_dist'],...
    [label '_' savelab '_Classifier_type' num2str(type) '_sdist'],'-append')

end

