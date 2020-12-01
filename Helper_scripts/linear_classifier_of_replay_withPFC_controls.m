function linear_classifier_of_replay_withPFC_controls(thisdir,label,Labels,replaypos,savelab,id)
try
%type 1 = linear discriminant analysis

%control 1 only taking ones that happen within the same median split of
%ripple power
%control 2 only taking ones that happen at the same arm at a time
%control 3 only taking ones that happen in the same half of the session at a time
load(thisdir,[label '_PFCreplayspikes_binned'],[label  '_replay_stnd'],...
    [label  '_replay_replayarm'],[label  '_replay_singlebothjoint'],'other_cells',[label '_Cand_sig_modu_include'],'HP_Ripple')
eval(['replayarm = ' label '_replay_replayarm;'])
eval(['stnd = ' label '_replay_stnd;'])
eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['PFC_binned = ' label '_PFCreplayspikes_binned;'])
eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
clear([label '_Cand_sig_modu_include'],[label '_PFCreplayspikes_binned'],[label '_replay_singlebothjoint'],[label '_replay_replayarm'])

pfc = other_cells(Cand_sig_modu_include(:,3)==1); clear other_cells;
PFC_binned(Cand_sig_modu_include(:,3)==0,:,:) = [];


replaypos= replaypos(singlebothjoint~=3);
stnd = stnd(singlebothjoint~=3,:);
Labels = Labels(singlebothjoint~=3);
PFC_binned = PFC_binned(:,:,singlebothjoint~=3);

binsize = .02;
window = [-2 2];
ind = [window(1)+(binsize/2):binsize:window(2)-(binsize/2)];
num_shuffles = 500;
% num_shuffles = 10;
mind = [0 .2]; bind = [-.5 -.1];
modind = ind>=mind(1) & ind<mind(2);
baseind = ind>=bind(1) & ind<bind(2);
% bindat = PFC_binned; 
Labelssave = Labels;


% seperately for smaller ripples and larger ripples
ripple_max = NaN(size(stnd,1),1);
for ievent = 1:size(stnd,1)
    ripple_max(ievent,1) = max(HP_Ripple(HP_Ripple(:,1)>=stnd(ievent,1) & HP_Ripple(:,1)<=stnd(ievent,2),2));
end
med = nanmedian(ripple_max);
ind1 = ripple_max<med;
ind2 = ripple_max>med;
if ~isempty(pfc) 
    tic
    itX = []; itSX = [];    
    for ihalf = 1:2
        if ihalf==1
            Labels = Labelssave(ind1);
            bindat = PFC_binned(:,:,ind1);
        elseif ihalf==2 
            Labels = Labelssave(ind2);
            bindat = PFC_binned(:,:,ind2);
        end
            
        indices = 1:length(Labels);
        num_folds = length(Labels);
        randss = randi(num_folds-1,num_folds-1,num_shuffles);    
        while any(sum(randss~=([1:length(Labels)-1]'))==0)
            randss = randi(num_folds-1,num_folds-1,num_shuffles);
        end
        iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles);
        bindat1 = squeeze((sum(bindat(:,modind,:),2)/sum(modind))-(sum(bindat(:,baseind,:),2))/sum(baseind));
        excl = sum(unique(bindat1','rows')~=0)<4; %added 5/20/19
        bindat2 = zscore(bindat1(~excl,:)');
        for i = 1:num_folds            
            test = (indices == i); train = ~test;
            traindata = bindat2(train,:);
            testdat = bindat2(test,:);   
            excl = std(traindata,[],1)<10^-10;
            traindata(:,excl) = []; testdat(:,excl) = [];
            guess = classify(testdat,traindata,Labels(train),'diaglinear');  %changed from linear only for day 11, on 10/25, revisit        
            real = Labels(test);
            iterX(i) = real==guess;
            ind2 = find(train);
            for j = 1:num_shuffles            
                randoth = randss(:,~ismember(1:size(randss,2),j));
                while any(sum(randoth~=randss(:,j))==0)
                    randss(:,j) = randi(num_folds-1,num_folds-1,1);
                end
                ss = randss(:,j);                
                guessS = classify(testdat,traindata,Labels(ind2(ss)),'diaglinear');    %changed from linear only for day 11, on 10/25, revisit                
                iterSX(i,j) = real==guessS;
            end

            if rem(i,500)==0
                t = toc;
               disp([num2str(i) ' Out of ' num2str(num_folds) ' Folds Finished in ' num2str(t/60) ' Minutes']) 
               tic
            end
        end
        itX = cat(1,itX,iterX);
        itSX = cat(1,itSX,iterSX);
    end
    r = sum(itX)/length(Labelssave);    
    s = sum(itSX,1)./length(Labelssave);
    p = (sum(s>=r)+1)/(length(s)+1);
    
elseif isempty(pfc) 
   p = NaN; r = NaN; s = NaN;

end
eval([[label '_' savelab '_ripple_Classifier_type1_p'] '= p;'])
eval([[label '_' savelab '_ripple_Classifier_type1_r'] '= r;'])
eval([[label '_' savelab '_ripple_Classifier_type1_s'] '= s;'])
save(thisdir,[label '_' savelab '_ripple_Classifier_type1_p'],...
    [label '_' savelab '_ripple_Classifier_type1_r'],...
    [label '_' savelab '_ripple_Classifier_type1_s'],'-append')

disp(['Ripple control p = ' num2str(p)])



% seperately for first half and second half (session time)
if ~isempty(pfc) 
    tic
    itX = []; itSX = [];
    midp = round(length(Labelssave)/2);
    for ihalf = 1:2
        if ihalf==1
            Labels = Labelssave(1:midp);
            bindat = PFC_binned(:,:,1:midp);
        elseif ihalf==2 
            Labels = Labelssave(midp+1:end);
            bindat = PFC_binned(:,:,midp+1:end);
        end
            
        indices = 1:length(Labels);
        num_folds = length(Labels);
        randss = randi(num_folds-1,num_folds-1,num_shuffles);    
        while any(sum(randss~=([1:length(Labels)-1]'))==0)
            randss = randi(num_folds-1,num_folds-1,num_shuffles);
        end
        iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles);
        bindat1 = squeeze((sum(bindat(:,modind,:),2)/sum(modind))-(sum(bindat(:,baseind,:),2))/sum(baseind));
        excl = sum(unique(bindat1','rows')~=0)<4;
        bindat2 = zscore(bindat1(~excl,:)');
        for i = 1:num_folds            
            test = (indices == i); train = ~test;
            traindata = bindat2(train,:);
            testdat = bindat2(test,:);   
            excl = std(traindata,[],1)<10^-10;
            traindata(:,excl) = []; testdat(:,excl) = [];
            guess = classify(testdat,traindata,Labels(train),'diaglinear');  %changed from linear only for day 11, on 10/25, revisit        
            real = Labels(test);
            iterX(i) = real==guess;
            ind2 = find(train);
            for j = 1:num_shuffles     
                randoth = randss(:,~ismember(1:size(randss,2),j));
                while any(sum(randoth~=randss(:,j))==0)
                    randss(:,j) = randi(num_folds-1,num_folds-1,1);
                end
                ss = randss(:,j);                
                guessS = classify(testdat,traindata,Labels(ind2(ss)),'diaglinear');    %changed from linear only for day 11, on 10/25, revisit                
                iterSX(i,j) = real==guessS;
            end

            if rem(i,500)==0
                t = toc;
               disp([num2str(i) ' Out of ' num2str(num_folds) ' Folds Finished in ' num2str(t/60) ' Minutes']) 
               tic
            end
        end
        itX = cat(1,itX,iterX);
        itSX = cat(1,itSX,iterSX);
    end
    r = sum(itX)/length(Labelssave);    
    s = sum(itSX,1)./length(Labelssave);
    p = (sum(s>=r)+1)/(length(s)+1);
    
elseif isempty(pfc) 
   p = NaN; r = NaN; s = NaN;

end
eval([[label '_' savelab '_time_Classifier_type1_p'] '= p;'])
eval([[label '_' savelab '_time_Classifier_type1_r'] '= r;'])
eval([[label '_' savelab '_time_Classifier_type1_s'] '= s;'])
save(thisdir,[label '_' savelab '_time_Classifier_type1_p'],...
    [label '_' savelab '_time_Classifier_type1_r'],...
    [label '_' savelab '_time_Classifier_type1_s'],'-append')

disp(['Time control p = ' num2str(p)])

% Seperately depending on which arm the animal is in space (armpos)
if ~isempty(pfc) 
    tic
    itX = []; itSX = [];    nf = 0;
    for iarm = 1:3
        
        Labels = Labelssave(replaypos==iarm);
        if length(unique(Labels))<2 || length(Labels)<60 % figure out how much data is too little data (20 too little)
            disp(['Not enough data for arm ' num2str(iarm)])
            continue
        end
        bindat = PFC_binned(:,:,replaypos==iarm);        
            
        indices = 1:length(Labels);
        num_folds = length(Labels);
        nf = nf+num_folds;
        randss = randi(num_folds-1,num_folds-1,num_shuffles);  
        while any(sum(randss~=([1:length(Labels)-1]'))==0)
            randss = randi(num_folds-1,num_folds-1,num_shuffles);
        end
        iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles);
%         bindat2 = squeeze((sum(bindat(:,modind,:),2)-sum(bindat(:,baseind,:),2))./(sum(bindat(:,modind,:),2)+sum(bindat(:,baseind,:),2)));
        bindat1 = squeeze((sum(bindat(:,modind,:),2)/sum(modind))-(sum(bindat(:,baseind,:),2))/sum(baseind));
        excl = sum(unique(bindat1','rows')~=0)<4;
        bindat2 = zscore(bindat1(~excl,:)');
        
        for i = 1:num_folds
            test = (indices == i); train = ~test;
            traindata = bindat2(train,:);
            testdat = bindat2(test,:);   
            excl = std(traindata,[],1)<10^-10;
            traindata(:,excl) = []; testdat(:,excl) = [];
            guess = classify(testdat,traindata,Labels(train),'diaglinear');  
            real = Labels(test);
            iterX(i) = real==guess;
            ind2 = find(train);
            for j = 1:num_shuffles     
                randoth = randss(:,~ismember(1:size(randss,2),j));
                while any(sum(randoth~=randss(:,j))==0)
                    randss(:,j) = randi(num_folds-1,num_folds-1,1);
                end
                ss = randss(:,j);                                
                guessS = classify(testdat,traindata,Labels(ind2(ss)),'diaglinear');   
                iterSX(i,j) = real==guessS;
            end

            if rem(i,500)==0
                t = toc;
               disp([num2str(i) ' Out of ' num2str(num_folds) ' Folds Finished in ' num2str(t/60) ' Minutes']) 
               tic
            end
        end
        itX = cat(1,itX,iterX);
        itSX = cat(1,itSX,iterSX);
    end
    r = sum(itX)/nf;    
    s = sum(itSX,1)./nf;
    p = (sum(s>=r)+1)/(length(s)+1);
    
elseif isempty(pfc) 
   p = NaN; r = NaN; s = NaN;
end

disp(['Arm control p = ' num2str(p)])
eval([[label '_' savelab '_armpos_Classifier_type1_p'] '= p;'])
eval([[label '_' savelab '_armpos_Classifier_type1_r'] '= r;'])
eval([[label '_' savelab '_armpos_Classifier_type1_s'] '= s;'])
save(thisdir,[label '_' savelab '_armpos_Classifier_type1_p'],...
    [label '_' savelab '_armpos_Classifier_type1_r'],...
    [label '_' savelab '_armpos_Classifier_type1_s'],'-append')

    


% replay length
ddiff = diff(stnd,[],2);
med = nanmedian(ddiff);
ind1 = ddiff<med;
ind2 = ddiff>med;
if ~isempty(pfc) 
    tic
    itX = []; itSX = [];    
    for ihalf = 1:2
        if ihalf==1
            Labels = Labelssave(ind1);
            bindat = PFC_binned(:,:,ind1);
        elseif ihalf==2 
            Labels = Labelssave(ind2);
            bindat = PFC_binned(:,:,ind2);
        end
            
        indices = 1:length(Labels);
        num_folds = length(Labels);
        randss = randi(num_folds-1,num_folds-1,num_shuffles);    
        iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles);
        bindat1 = squeeze((sum(bindat(:,modind,:),2)/sum(modind))-(sum(bindat(:,baseind,:),2))/sum(baseind));
        excl = sum(unique(bindat1','rows')~=0)<4; %added 5/20/19
        bindat2 = zscore(bindat1(~excl,:)');
        for i = 1:num_folds            
            test = (indices == i); train = ~test;
            traindata = bindat2(train,:);
            testdat = bindat2(test,:);   
            excl = std(traindata,[],1)<10^-10;
            traindata(:,excl) = []; testdat(:,excl) = [];
            guess = classify(testdat,traindata,Labels(train),'diaglinear');  %changed from linear only for day 11, on 10/25, revisit        
            real = Labels(test);
            iterX(i) = real==guess;
            ind2 = find(train);
            for j = 1:num_shuffles            
                ss = randss(:,j);                
                guessS = classify(testdat,traindata,Labels(ind2(ss)),'diaglinear');    %changed from linear only for day 11, on 10/25, revisit                
                iterSX(i,j) = real==guessS;
            end

            if rem(i,500)==0
                t = toc;
               disp([num2str(i) ' Out of ' num2str(num_folds) ' Folds Finished in ' num2str(t/60) ' Minutes']) 
               tic
            end
        end
        itX = cat(1,itX,iterX);
        itSX = cat(1,itSX,iterSX);
    end
    r = sum(itX)/length(Labelssave);    
    s = sum(itSX,1)./length(Labelssave);
    p = (sum(s>=r)+1)/(length(s)+1);
    
elseif isempty(pfc) 
   p = NaN; r = NaN; s = NaN;

end
eval([[label '_' savelab '_replaylength_Classifier_type1_p'] '= p;'])
eval([[label '_' savelab '_replaylength_Classifier_type1_r'] '= r;'])
eval([[label '_' savelab '_replaylength_Classifier_type1_s'] '= s;'])
save(thisdir,[label '_' savelab '_replaylength_Classifier_type1_p'],...
    [label '_' savelab '_replaylength_Classifier_type1_r'],...
    [label '_' savelab '_replaylength_Classifier_type1_s'],'-append')

disp(['Replay length control p = ' num2str(p)])
 catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(id)
        disp('****************************************Error Occured')
end