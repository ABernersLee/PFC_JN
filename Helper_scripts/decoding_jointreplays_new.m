% num_shuffles = 10;
% num_hp_shuffles = 10;

label = 'RP';


cd(dirs.homedir)
d2 = dir('*.mat');

realall = []; shuffall = [];
numshuff= 500;
realacc = [];
shuffacc = [];
realacc2 = [];
shuffacc2 = [];
pid = NaN(size(d2,1),1);

for id = 1:size(d2,1)
        
        
    shuff_acc = [];
    real_acc = [];
    thisdir = d2(id).name;
    load(thisdir,[label '_replay_stnd'],'other_cells','hp_cells','hpinterneurons','spikedata',[label '_replay_singlebothjoint'],...
        [label  '_replay_jointreplayarm'],[label '_Cand_sig_modu_include'],[label  '_replay_jointjumptime'],[label  '_replay_post'],[label  '_replay_postseqindex'])

    eval(['singlebothjoint = ' label '_replay_singlebothjoint;'])
    eval(['Event = ' label '_replay_stnd(:,1);'])
    eval(['JointArm = ' label  '_replay_jointreplayarm;'])
    eval(['JointTime = ' label  '_replay_jointjumptime;'])
    eval(['Cand_sig_modu_include = ' label '_Cand_sig_modu_include;'])
    clear([label '_replay_stnd'])

%     pfc = other_cells(Cand_sig_modu_include(:,1)==1); 
    pfc = other_cells; 

    hp = setdiff(hp_cells,hpinterneurons);
    
    cellnums = 1:length(hp); cellnums = cellnums(randperm(length(cellnums)));
    hp = hp(cellnums(1:length(pfc)));

    clear other_cells hp_cells hpinterneurons
    Event(singlebothjoint~=3) = [];
    JointArm(singlebothjoint~=3,:) = [];
    JointTime(singlebothjoint~=3,:) = [];

    PFCspk = spikedata(ismember(spikedata(:,2),pfc),1:2); 
    HPspk = spikedata(ismember(spikedata(:,2),hp),1:2); 

    FRpfc = NaN(size(Event,1),length(pfc));
    FRhp = NaN(size(Event,1),length(hp));
    spikewindow = .2;
    for iE = 1:size(Event,1)
        for icell = 1:length(pfc)
            FRpfc(iE,icell) = sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)>=(Event(iE)-spikewindow) & PFCspk(:,1)<JointTime(iE))./((JointTime(iE)-Event(iE)+spikewindow));      %!! changed 10/22 to change in FR 
%             FRpfc(iE,icell) = (sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)>=Event(iE) & PFCspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)))-...
%                 (sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)<Event(iE) & PFCspk(:,1)>=(Event(iE)-.2))./.2); 
    
%             FRpfc(iE,icell) = sum(PFCspk(:,2)==pfc(icell) & PFCspk(:,1)>=(Event(iE)-spikewindow) & PFCspk(:,1)<Event(iE))./(spikewindow);     
        end
%         for icell = 1:length(hp)
%             FRhp(iE,icell) = sum(HPspk(:,2)==hp(icell) & HPspk(:,1)>=(Event(iE)-spikewindow) & HPspk(:,1)<JointTime(iE))./((JointTime(iE)-Event(iE))+spikewindow); % !!changed 10/22 to change in FR  
% %     %         FRhp(iE,icell) = (sum(HPspk(:,2)==hp(icell) & HPspk(:,1)>=Event(iE) & HPspk(:,1)<JointTime(iE))./(JointTime(iE)-Event(iE)))- ...
% %     %         (sum(HPspk(:,2)==hp(icell) & HPspk(:,1)<Event(iE) & HPspk(:,1)>=(Event(iE)-.2))./.2);
% %                 FRhp(iE,icell) = sum(HPspk(:,2)==hp(icell) & HPspk(:,1)>=(Event(iE)-spikewindow) & HPspk(:,1)<Event(iE))./(spikewindow);
%         end
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
        HParm = zscore(HParm);
        

        pfc = PFCarm;
        
        if all(histc(Labels,unique(Labels))>1) && length(unique(Labels))>1
            folds = 10;
%             folds = length(Labels);
            jnk = ones(ceil(length(Labels)/folds),1)*[1:folds];
            inx = jnk(:); inx = inx(randperm(length(inx))); inx = inx(1:length(Labels));

            for ifold = 1:folds
                train = inx~=ifold; test = inx==ifold;
                traindata = pfc(train,:); testdata = pfc(test,:);
                trainlabels = Labels(train);
                Mdl = fitcdiscr(traindata,trainlabels,'discrimType','pseudoLinear');                        
                guess = predict(Mdl,testdata);
                corrReal = guess==Labels(test);
                real_acc = [real_acc; corrReal]; 
                corrShuff = NaN(length(corrReal),numshuff);

                for ishuff = 1:numshuff
                    shufflabels = trainlabels(randperm(length(trainlabels)));                                
                    Mdl = fitcdiscr(traindata,shufflabels,'discrimType','pseudoLinear');                        
                    guess = predict(Mdl,testdata);
                    corrShuff(:,ishuff) = guess==Labels(test);                                
                end
                shuff_acc = [shuff_acc; corrShuff];

            end
        else
            disp(['skip arm ' num2str(iarm) ' day ' num2str(id)])
        end
    end
    pid(id) = (sum(nanmean(shuff_acc)>=nanmean(real_acc))+1)/(size(shuff_acc,2)+1);
    realacc = [realacc;nanmean(real_acc)];
    realacc2 = [realacc2;real_acc];
    shuffacc = [shuffacc;nanmean(shuff_acc)];
    shuffacc2 = [shuffacc2;shuff_acc];
    disp(['Done with day ' num2str(id) ' p ' num2str(pid(id)) ' # cells: ' num2str(size(pfc,2))])
end

p = (sum(nanmean(shuffacc)>=nanmean(realacc))+1)/(size(shuffacc,2)+1);
figure; histogram(nanmean(shuffacc),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realacc) nanmean(realacc)],yl,'r-','LineWidth',3)
text(nanmean(realacc),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title('Session Averages')

p = (sum(nanmean(shuffacc2)>=nanmean(realacc2))+1)/(size(shuffacc2,2)+1);
figure; histogram(nanmean(shuffacc2),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realacc2) nanmean(realacc2)],yl,'r-','LineWidth',3)
text(nanmean(realacc2),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
title('All Replays')