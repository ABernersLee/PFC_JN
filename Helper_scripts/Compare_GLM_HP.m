% function Compare_GLM_HP(dirs)

cd(dirs.spikedatadir)
d2 = dir('*.mat');
torun = 1:size(d2,1);

bn = [.02; .04; .06; .1; .2; .4; .6];
% bn = [.2; .8];
bn = bn(end:-1:1);
numpredict = 9; %should be same as size(X,2)




num_folds = 10;
Rall = []; Rpopall = []; Rallid = NaN(length(torun),numpredict,length(bn));
for id = 1:length(torun)
   load(d2(torun(id)).name,'hpinterneurons','hp_cells','spikedata','vel','pos','dirdat','linposcat','armpos','laps_singlepass','behave_change_log')
    
    other_cells = hp_cells(~ismember(hp_cells,hpinterneurons));
    
    dc = linposcat;
    for iarm = 1:max(armpos)
        dc(armpos==iarm) = dc(armpos==iarm)-min(dc(armpos==iarm));
    end
    
    
    spks = spikedata(ismember(spikedata(:,2),other_cells),:);
    
    armacc = get_behavior_accuracy(d2(torun(id)).name);    
    iscorrectlap = false(size(laps_singlepass));
    iscorrectlap(ismember(laps_singlepass,find(armacc(:,3)==1))) = true;
    
    R = NaN(length(other_cells),numpredict,length(bn),num_folds); %R(icell,itype,ib,ifold)
    RPop = NaN(1,numpredict,length(bn),num_folds); 
    for ib = 1:length(bn)

%         pos1 = []; 
        a1 = []; fa = []; v = []; d = []; dc1 = []; iscor = []; tt = []; ts = []; lap_fr = []; dhr = []; dt = [];
        for ilap = 1:max(laps_singlepass)
             leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'last'); %1 leave middle 5 leave platform
            if isempty(leave)
               leave = find(laps_singlepass==ilap,1,'first');
            end
           lapind = leave:find(laps_singlepass==ilap,1,'last');

            posdat1 = dc(lapind); %distance to heading reward
            posdat3 = posdat1; %distance from any reward
            
            laparm = armpos(lapind);        
            lastarm = laparm(end);
            otharms = setdiff(1:3,lastarm);
            for ioth = 1:2
                posdat1(laparm==otharms(ioth)) = max(posdat1(laparm==lastarm))+(posdat1(laparm==otharms(ioth))-min(posdat1(laparm==otharms(ioth))));
            end        
            posdat1(laparm==lastarm) = max(posdat1(laparm==lastarm))-posdat1(laparm==lastarm);
            posdat2 = [0;cumsum(-diff(posdat1))];  %distance traveled
            if sum(posdat2<0)>0
%                 if ib==1
%                     disp(['problem for ' num2str(sum(posdat2<0))])
%                 end
                posdat2 = posdat2-min(posdat2);
            end
            timeto = pos(lapind(end),1)-pos(lapind,1);
            timesince = pos(lapind,1)-pos(lapind(1),1);

            futurearm = ones(size(lapind))*lastarm;
%             futurearm(laparm~=lastarm) = lastarm;
%             futurearm(laparm==lastarm) = armpos(find(armpos(lapind(end)+1:end)~=lastarm,1,'first')+lapind(end));


           timebins = pos(lapind(1),1):bn(ib):pos(lapind(end),1)+bn(ib);
            
           s = spks(spks(:,1)>=pos(lapind(1),1) & spks(:,1)<=pos(lapind(end),1),1:2);
           [~,~,i] = histcounts(s(:,1),timebins);
           fr = zeros(length(timebins)-1,length(other_cells));
           for icell = 1:length(other_cells)
               if sum(s(:,2)==other_cells(icell))>0
                  fr(:,icell) = histc(i(s(:,2)==other_cells(icell)),1:length(timebins)-1);
               end
           end
                
%             timebins = pos(lapind(1),1):bn(ib):pos(lapind(end),1);
            [~,~,i] = histcounts(pos(lapind,1),timebins);
            laparmX = NaN(max(i),1); 
            futurearmX = laparmX; velX = laparmX; dirdatX = laparmX; posdat3X= laparmX;
            posdat1X = laparmX; posdat2X = laparmX; iscorrectlapX = laparmX; timetoX = laparmX; timesinceX = laparmX;
            for ii = 1:max(i)
                if sum(i==ii)>0
                    laparmX(ii) = mode(laparm(i==ii));
                    futurearmX(ii) = mode(futurearm(i==ii));
                    velX(ii) = mean(vel(i==ii));
                    dirdatX(ii) = mode(dirdat(i==ii));
                    posdat3X(ii) = mean(posdat3(i==ii));
                    posdat1X(ii) = mean(posdat1(i==ii));
                    posdat2X(ii) = mean(posdat2(i==ii));
                    iscorrectlapX(ii) = mode(iscorrectlap(i==ii));
                    timetoX(ii) = mean(timeto(i==ii));
                    timesinceX(ii) = mean(timesince(i==ii));
                end
            end
            exl = find(isnan(posdat3X));
%             disp(['Exclude ' num2str(length(exl)) ', ' num2str(bn(ib))])
            
            laparmX(exl) = [];
            futurearmX(exl) = [];
            velX(exl) = [];
            dirdatX(exl) = [];
            posdat3X(exl) = [];
            posdat1X(exl) = [];
            posdat2X(exl) = [];
            iscorrectlapX(exl) = [];
            timetoX(exl) = [];
            timesinceX(exl) = [];
            fr(exl,:) = [];
            
%             pos1 = [pos1;pos(lapind(i),1)]; 
            a1 = [a1;laparmX]; fa = [fa;futurearmX]; v = [v;velX]; 
            d = [d;dirdatX]; dc1 = [dc1;posdat3X]; dhr = [dhr;posdat1X]; dt = [dt; posdat2X];
            iscor = [iscor;iscorrectlapX]; tt = [tt;timetoX]; ts = [ts;timesinceX]; 
            lap_fr = cat(1,lap_fr,fr);
            if size(fr,1)~=length(laparmX)
                disp('error')
            end
        end


%             X = [a1 fa iscor d v dc1 dhr dt tt ts];
        X = [a1 iscor d v dc1 dhr dt tt ts];

        indices = crossvalind('Kfold',size(X,1),num_folds)';             
        for ifold = 1:num_folds
            test = (indices == ifold); train = ~test;

            
            for itype = 1:size(X,2)
                if itype< 4
                    mdl = fitglm(lap_fr(train,:),X(train,itype),'linear','CategoricalVars',itype);
                else
                    mdl = fitglm(lap_fr(train,:),X(train,itype),'linear');
                end
                [ypred,~] = predict(mdl,lap_fr(test,:));
                RPop(1,itype,ib,ifold) = corr(X(test,itype),double(ypred)).^2;
                
                for icell = 1:size(lap_fr,2)
                    if sum(lap_fr(train,icell))==0
                        continue
                    end
                    if itype < 4 %if arm, future arm, correct/incorrect, direction
                            lm = fitcdiscr(lap_fr(train,icell),categorical(X(train,itype))); %predictor, response                            
                            [ypred,~] = predict(lm,lap_fr(test,icell));
                            
                    else % if velocity or distance from center, time to, time since
%                             lm = fitlm(X(train,itype),lap_fr(train,icell));
%                             [ypred,~] = predict(lm,X(test,itype));
                        lm = fitlm(lap_fr(train,icell),X(train,itype));
                        [ypred,~] = predict(lm,lap_fr(test,icell));
                    end


                    R(icell,itype,ib,ifold) = corr(X(test,itype),double(ypred)).^2;
                end
            end

        end


        disp(['Done id ' num2str(id) ' ib ' num2str(ib)])
    end
    
    Rpopall = cat(1,Rpopall,RPop);
    Rall = cat(1,Rall,nanmean(R,4));
    Rallid(id,:,:) = nanmean(nanmean(R,4),1);
end


figlab = 'F:\XY_matdata\Figures\ForPaperReviews\AllCells_AllArmEvents\GLM\';

clear a
figure; hold on
cc = varycolor(size(X,2));
for itype = [1:size(X,2)]
%     for ib = 1:length(bn)
        m = squeeze(nanmean(Rall(:,itype,:),1));
        sem = squeeze(nanstd(Rall(:,itype,:)))./squeeze(sqrt(sum(~isnan(Rall(:,itype,:)))));
        a(itype) = errorbar(1:size(Rall,3),m,sem,'LineWidth',2,'Color',cc(itype,:));
%     end    
end
 legend(a,{'Arm';'Corr/Incorr';'Direction';'Speed';'Distance from center';'Distance to heading reward';'Distance traveled';'Time to Reward';'Time since left'},'Location','Best')
set(gca,'xtick',1:length(bn),'xticklabel',bn)
tit = ['10-fold crossval, all HP neruons (N = ' num2str(size(Rall,1)) ')'];
title(tit)
xlabel('Time Window/Bin (Seconds)')
ylabel('R-squared')
set(gca,'FontSize',18)
set(gcf,'Position',[    2006          68         906         696])
helper_saveandclosefig([figlab '\CompLinearModels_' num2str(num_folds) 'fold_Neurons_N' num2str(size(Rall,1)) '_HP'])



clear a
figure; hold on
cc = varycolor(size(X,2));
for itype = [1:size(X,2)]
%     for ib = 1:length(bn)
        m = squeeze(nanmean(Rallid(:,itype,:),1));
        sem = squeeze(nanstd(Rallid(:,itype,:)))./squeeze(sqrt(sum(~isnan(Rallid(:,itype,:)))));
        a(itype) = errorbar(1:size(Rallid,3),m,sem,'LineWidth',2,'Color',cc(itype,:));
%     end    
end
 legend(a,{'Arm';'Corr/Incorr';'Direction';'Speed';'Distance from center';'Distance to heading reward';'Distance traveled';'Time to Reward';'Time since left'},'Location','Best')
set(gca,'xtick',1:length(bn),'xticklabel',bn)
tit = ['10-fold crossval, all days (N = ' num2str(size(Rallid,1)) ')'];
title(tit)
xlabel('Time Window/Bin (Seconds)')
ylabel('R-squared')
set(gca,'FontSize',18)
set(gcf,'Position',[    2006          68         906         696])
helper_saveandclosefig([figlab '\CompLinearModels_' num2str(num_folds) 'fold_Dayd_N' num2str(size(Rallid,1)) '_HP'])



Rpopall2 = nanmean(Rpopall,4);
clear a
figure; hold on
cc = varycolor(size(X,2));
for itype = [1:size(X,2)]
%     for ib = 1:length(bn)
        m = squeeze(nanmean(Rpopall2(:,itype,:),1));
        sem = squeeze(nanstd(Rpopall2(:,itype,:)))./squeeze(sqrt(sum(~isnan(Rpopall2(:,itype,:)))));
        a(itype) = errorbar(1:size(Rpopall2,3),m,sem,'LineWidth',2,'Color',cc(itype,:));
%     end    
end
 legend(a,{'Arm';'Corr/Incorr';'Direction';'Speed';'Distance from center';'Distance to heading reward';'Distance traveled';'Time to Reward';'Time since left'},'Location','Best')
set(gca,'xtick',1:length(bn),'xticklabel',bn)
tit = ['10-fold crossval, GLM all HP neruons (N = ' num2str(size(Rall,1)) ')'];
title(tit)
xlabel('Time Window/Bin (Seconds)')
ylabel('R-squared')
set(gca,'FontSize',18)
set(gcf,'Position',[    2006          68         906         696])
helper_saveandclosefig([figlab '\CompLinearModels_GLM_' num2str(num_folds) 'fold_Neurons_N' num2str(size(Rpopall2,1)) '_HP'])



Rpopall2 = nanmean(Rpopall,4);
a = [];
figure; hold on
cc = varycolor(size(X,2));
for itype = [1:5 7 9]
%     for ib = 1:length(bn)
        m = squeeze(nanmean(Rpopall2(:,itype,:),1));
        sem = squeeze(nanstd(Rpopall2(:,itype,:)))./squeeze(sqrt(sum(~isnan(Rpopall2(:,itype,:)))));
        aa = errorbar(1:size(Rpopall2,3),m,sem,'LineWidth',2,'Color',cc(itype,:));
        a = cat(1,a,aa); 
%     end    
end
 legend(a,{'Arm';'Corr/Incorr';'Direction';'Speed';'Distance from center';'Distance traveled';'Time since left'},'Location','Best')
set(gca,'xtick',1:length(bn),'xticklabel',bn)
tit = ['10-fold crossval, GLM all HP neruons (N = ' num2str(size(Rall,1)) ')'];
title(tit)
xlabel('Time Window/Bin (Seconds)')
ylabel('R-squared')
set(gca,'FontSize',18)
set(gcf,'Position',[    2006          68         906         696])
set(gcf,'renderer','Painters')
helper_saveandclosefig([figlab '\CompLinearModels_GLM_fewer_' num2str(num_folds) 'fold_Neurons_N' num2str(size(Rpopall2,1)) '_HP'])
