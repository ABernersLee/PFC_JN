function decoding_jointreplays(thisdir,label,num_shuffles,num_hp_shuffles)

%formerly make_jointreplayeventtriggeredmat
load(thisdir,[label '_replay_stnd'],'other_cells','hp_cells','hpinterneurons','spikedata',[label '_replay_singlebothjoint'],...
    [label  '_replay_jointreplayarm'],[label  '_replay_jointjumptime'],[label  '_replay_post'],[label  '_replay_postseqindex'])

eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
eval(['Event = ' label '_replay_stnd(:,1);'])
eval(['JointArm = ' label  '_replay_jointreplayarm;'])
eval(['JointTime = ' label  '_replay_jointjumptime;'])
clear([label '_replay_stnd'])

pfc = other_cells; 

hp = setdiff(hp_cells,hpinterneurons);

clear other_cells hp_cells hpinterneurons
Event(singlebothjoint~=3) = [];
JointArm(singlebothjoint~=3,:) = [];
JointTime(singlebothjoint~=3,:) = [];

PFCspk = spikedata(ismember(spikedata(:,2),pfc),1:2); 
HPspk = spikedata(ismember(spikedata(:,2),hp),1:2); 

FRpfc = NaN(size(Event,1),length(pfc));
FRhp = NaN(size(Event,1),length(hp));
for iE = 1:size(Event,1)
    for icell = 1:length(pfc)
        FRpfc(iE,icell) = sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)>=Event(iE) & PFCspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE));      %changed 10/22 to change in FR 
%         FRpfc(iE,icell) = (sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)>=Event(iE) & PFCspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)))-...
%             (sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)<Event(iE) & PFCspk(:,1)>=(Event(iE)-.2))./.2); 
    end
    for icell = 1:length(hp)
        FRhp(iE,icell) = sum(HPspk(:,2)==hp(icell) & HPspk(:,1)>=Event(iE) & HPspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)); %changed 10/22 to change in FR  
%         FRhp(iE,icell) = (sum(HPspk(:,2)==hp(icell) & HPspk(:,1)>=Event(iE) & HPspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)))- ...
%         (sum(HPspk(:,2)==hp(icell) & HPspk(:,1)<Event(iE) & HPspk(:,1)>=(Event(iE)-.2))./.2);
    end
end


PFCreal =[];
PFCshuff = [];
HPreal =[];
HPshuff = [];
Armind = [];
CoeffsPFC = NaN(size(FRpfc,2),3); 
CoeffsHP = NaN(size(FRhp,2),3); 
hs = NaN(3,3);
for iarm = 1:3
    
    touse = JointArm(:,1)==iarm; % & sum(JFR,2)~=0;
    Labels = JointArm(touse,2);
    PFCarm = FRpfc(touse,:);
    PFCind = sum(PFCarm)~=0;
    PFCarm(:,sum(PFCarm)==0) = [];
    PFCarm = zscore(PFCarm);
    
    HParm = FRhp(touse,:);
    HPind = find(sum(HParm)~=0);
    HParm(:,sum(HParm)==0) = [];
    
    num_folds = length(Labels);
    iterX = NaN(num_folds,1); iterSX = NaN(num_folds,num_shuffles*round(num_shuffles/num_hp_shuffles));
    iterXH = NaN(num_folds,num_hp_shuffles); iterSXH = NaN(num_folds,round(num_shuffles/num_hp_shuffles),num_shuffles);
    
    hs(iarm,1) = length(Labels);
    h = hist(Labels,setdiff(1:3,iarm));
    hs(iarm,2) = min(h);    
    randss = randi(num_folds-1,num_folds-1,num_shuffles*round(num_shuffles/num_hp_shuffles));
    CoeffsPFC1 = NaN(size(FRpfc,2),num_folds);
    CoeffsHP1 = NaN(size(FRhp,2),num_folds,num_hp_shuffles);
        if size(PFCarm,2)<2  || hs(iarm,2)<2 %%|| length(Labels)<20 % length(Labels)<size(PFCarm,2) ||
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
                    
                    thisHP = zscore(thisHP);
                    
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
    disp(['Done with Arm ' num2str(iarm) ' in decoding_jointreplays'])

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
eval([label  '_decodejointreplays_PFCreal = PFCreal;'])
eval([label  '_decodejointreplays_PFCshuff = PFCshuff;'])
eval([label  '_decodejointreplays_HPreal = HPreal;'])
eval([label  '_decodejointreplays_HPshuff = HPshuff;'])
eval([label  '_decodejointreplays_CoeffsHP = CoeffsHP;'])
eval([label  '_decodejointreplays_CoeffsPFC = CoeffsPFC;'])
eval([label  '_decodejointreplays_ps = ps;'])

save(thisdir,[label  '_decodejointreplays_PFCreal'],[label  '_decodejointreplays_ps'],[label  '_decodejointreplays_PFCshuff'],[label  '_decodejointreplays_HPreal']...
    ,[label  '_decodejointreplays_HPshuff'],[label  '_decodejointreplays_CoeffsHP'],[label  '_decodejointreplays_CoeffsPFC'],'-append')

disp(['Done with decoding_jointreplays ' label])