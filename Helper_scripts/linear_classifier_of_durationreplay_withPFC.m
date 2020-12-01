function p2 =  linear_classifier_of_durationreplay_withPFC(thisdir,label)

load(thisdir,[label '_PFCreplayspikes_binned'],[label '_PFCreplayspikes_list'],...
    [label  '_replay_stnd'],[label  '_replay_singlebothjoint'],'other_cells',[label '_Cand_sig_modu_include'])
eval(['replaystnd = ' label '_replay_stnd;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['PFC_binned = ' label '_PFCreplayspikes_binned;'])
eval(['PFC_list = ' label '_PFCreplayspikes_list;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
clear([label '_Cand_sig_modu_include'],[label '_PFCreplayspikes_list'],[label '_PFCreplayspikes_binned'],[label '_replay_singlebothjoint'],[label '_replay_replayarm'])
touse = Cand_sig_modu_include(:,3)==1 & Cand_sig_modu_include(:,1)==1;
pfc = other_cells(touse); clear other_cells;
PFC_binned(~touse,:,:) = [];
dur = abs(diff(replaystnd'));

binsize = .002;
window = [-.5 .5];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];

% Multinomial logistic regression
modind = ind>.0 & ind<=.3;
X= NaN(length(pfc),sum(singlebothjoint==1));
for icell = 1:length(pfc)
    X(icell,:) =  squeeze(sum(PFC_binned(icell,modind,singlebothjoint==1),2));
end
X = X';
Labels = dur(singlebothjoint==1);
num_shuffles = 10;
num_folds = 10;
Fit = NaN(num_shuffles,num_folds);
Fit2 = Fit;
for j = 1:num_shuffles
    indices = crossvalind('Kfold',Labels,num_folds);
    Labels2 = Labels(randperm(length(Labels)));
    for i = 1:num_folds
        test = (indices == i); train = ~test;
        B = glmfit(X(train,:),Labels(train),'poisson');
        guess = glmval(B,X(test,:),'log');        
        Fit(j,i) = sum((guess'-dur(test)).^2);
        
        B = glmfit(X(train,:),Labels2(train),'poisson');
        guess = glmval(B,X(test,:),'log');        
        Fit2(j,i) = sum((guess'-dur(test)).^2);
    end
end

if mean(Fit2(:))>mean(Fit(:))
    p2 = signrank(Fit(:),Fit2(:));
else
    p2 = 1;
end

% if length(pfc)<4
%     p2 = NaN;
% end