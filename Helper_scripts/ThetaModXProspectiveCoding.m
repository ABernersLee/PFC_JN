
cd(dirs.homedir)
d2 = dir('*.mat');
lapinfo = cell(size(d2,1),1);
for id = 1:size(d2,1)
    thisdir = d2(id).name;

    load(thisdir,'vel','behave_change_log','laps_singlepass','armpos','pos','armposindex','dirdat','binpos')
    [~,armp] = max(armposindex,[],2); clear armposindex
    lappro = NaN(size(laps_singlepass,1),3);
    lapout = NaN(size(laps_singlepass,1),2);
    % lapdir = NaN(size(laps_singlepass,1),1);
    lapi = NaN(max(laps_singlepass),2);
    for ilap = 1:max(laps_singlepass)
       leave = find(behave_change_log(:,5) & laps_singlepass == ilap,1,'first'); %1 leave middle 5 leave platform
       lapind = leave:find(behave_change_log(:,6) & laps_singlepass == ilap,1,'first'); % 6 is arrive middle      
       lapout1 = find(behave_change_log(:,1) & laps_singlepass == ilap,1,'first'):find(behave_change_log(:,2) & laps_singlepass == ilap,1,'first'); % 2 is arrive platform

       nextarm = armpos(find(laps_singlepass == ilap,1,'last'));
       thisarm = armpos(find(laps_singlepass == ilap,1,'first'));
       lapi(ilap,:) = [thisarm nextarm];
    end
    lapinfo{id} = lapi;
end

%% Does HP tend to have theta sweeps of where the animal is about to go? yes
cd(dirs.homedir)
d2 = dir('*.mat');
pid = NaN(size(d2,1),1);
realall = []; shuffall = [];
numshuff= 1000;

for id = 1:size(d2,1)
    lapi = lapinfo{id};
    thisdir = d2(id).name;
    load(thisdir,'times_armon_thetaof_headingarm_lap')
    aa = times_armon_thetaof_headingarm_lap;
    realprop = NaN(max(aa(:,end)),1);
    for ilap = 1:max(aa(:,end))
        dat = aa(aa(:,6)==ilap,1:5);
        a = sum(dat(:,4)==unique(dat(:,5))); b = sum(dat(:,4)~=unique(dat(:,5)));
        realprop(ilap,1) = (a-b)./(a+b);
    end
    shufflap = [];
    for iarm = 1:3
        shufflaps = lapi(lapi(:,1)==iarm,2); 
        currlaps = unique(aa(aa(:,3)==iarm,6));
        shufflap1 = NaN(length(shufflaps),numshuff);
        for is = 1:numshuff
            sslaps = shufflaps(randperm(length(shufflaps)));
            for ilap = 1:length(currlaps)
                a = sum(aa(aa(:,6)==currlaps(ilap),4)==sslaps(ilap)); b = sum(aa(aa(:,6)==currlaps(ilap),4)~=sslaps(ilap));
                shufflap1(ilap,is) =  (a-b)./(a+b);
            end
        end
        shufflap = [shufflap;shufflap1];
    end
    pid(id) = (sum(nanmean(shufflap)>=nanmean(realprop))+1)/(size(shufflap,2)+1);
    realall = [realall; nanmean(realprop)];
    shuffall = [shuffall;nanmean(shufflap)];
end

p = (sum(nanmean(shuffall)>=nanmean(realall))+1)/(size(shuffall,2)+1);
figure; histogram(nanmean(shuffall),'FaceColor','k'); 
hold on; 
yl = get(gca,'ylim');
plot([nanmean(realall) nanmean(realall)],yl,'r-','LineWidth',3)
text(nanmean(realall),yl(2)*.8,['      p = ' num2str(round(p,2,'significant'))])
%% what can you decode - heading arm better than theta sweep - except tested arm2 and it is opposite
cd(dirs.homedir)
d2 = dir('*.mat');

realall = []; shuffall = [];
numshuff= 200;
realacc = [];
shuffacc = [];
pid = NaN(size(d2,1),1);

for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap')
    aa = times_armon_thetaof_headingarm_lap;
    pfcfr = zscore(prsopectiveTheta_pfcfr);
    shuff_acc = [];
    real_acc = [];
    for iarm = 1:3
        
        Labels = aa(aa(:,3)==iarm,5); % 4 is theta of, 5 is heading arm
        pfc = pfcfr(aa(:,3)==iarm,:);
        %leave one out
        thetacycles = 1:size(pfc,1);
        corrReal = NaN(max(thetacycles),1);
        corrShuff = NaN(max(thetacycles),numshuff);
%         for itheta = 1:max(thetacycles)
%             train = thetacycles~=itheta; test = thetacycles==itheta;
%         folds = max(thetacycles); 
        folds = 10;
        
        jnk = ones(ceil(max(thetacycles)/folds),1)*[1:folds];
        inx = jnk(:); inx = inx(randperm(length(inx))); inx = inx(1:max(thetacycles));
        for itheta = 1:folds
            train = inx~=itheta; test = inx==itheta;
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
        
        
        
        disp(['Done with arm ' num2str(iarm)])
    end
    pid(id) = (sum(nanmean(shuff_acc)>nanmean(real_acc))+1)/(size(shuff_acc,2)+1);
    realacc = [realacc;real_acc];
    shuffacc = [shuffacc;shuff_acc];
    disp(['Done with day ' num2str(id)])
end
p = (sum(nanmean(shuffacc)>nanmean(realacc))+1)/(size(shuffacc,2)+1)


%% ridge regression, theta is better explanation on center arm, and sig more so than other arms
% label = 'RP'; (just using sig mod cells doesnt change it much)
cd(dirs.homedir)
d2 = dir('*.mat');
k = 1e-5; %0:1e-5:5e-3;
bsall = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap')
    pfcfr = prsopectiveTheta_pfcfr;    
%     load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap',[label '_Cand_sig_modu_include'])
%     eval(['Cand = ' label '_Cand_sig_modu_include;'])
%     pfcfr = prsopectiveTheta_pfcfr(:,Cand(:,1)==1);    
    
    
    aa = times_armon_thetaof_headingarm_lap;
%     bs = NaN(size(pfc,2),3,length(k),3);
%     bs = NaN(size(pfc,2),3,3); 
    bs = NaN(size(pfc,2),2,3); 
    for iarm = 1:3
        
        Y = aa(aa(:,3)==iarm,4:5); % 4 is theta of, 5 is heading arm
        D = x2fx(Y,'interaction');
        pfc = pfcfr(aa(:,3)==iarm,:); pfc2 = zscore(pfc);
        for icell = 1:size(pfc,2)
            if sum(pfc(:,icell))>0                
%                 [b,~,~] = glmfit(Y,pfc2(:,icell));        
%                 b = ridge(pfc(:,icell),D(:,2:end),k);                        
                b = ridge(pfc(:,icell),Y,k);        
%                 bs(icell,:,:,iarm) = b;
                bs(icell,:,iarm) = b;
            end
        end
        disp(['Done with arm ' num2str(iarm)])
    end    
    bsall = cat(1,bsall,bs);
    disp(['Done with day ' num2str(id)])
end
%%
bsall = abs(bsall);
%% results from single factor test
 [p,tbl,stats]= kruskalwallis(bsall(:,1,:)-bsall(:,2,:)); % sig difference across arms
 p = multcompare(stats) %center arm is sig more than other two arms
 
 p1 = ranksum(bsall(:,1,1),bsall(:,2,1))
 nanmedian(bsall(:,1,1))-nanmedian(bsall(:,2,1)) % left arm is non-significantly more heading than theta
 
 p3 = ranksum(bsall(:,1,3),bsall(:,2,3))
 nanmedian(bsall(:,1,3))-nanmedian(bsall(:,2,3)) % right arm is non-significantly more heading than theta (switches but not sig when abs(beta))
 
 p1 = ranksum(bsall(:,1,2),bsall(:,2,2))
 nanmedian(bsall(:,1,2))-nanmedian(bsall(:,2,2)) % center arm is significantly more about theta than heading

 %% results from interaction test

 [p,tbl,stats]= kruskalwallis(nanmean(bsall,3)); % no difference when average over arms
 
[p,tbl,stats]= kruskalwallis(bsall(:,:,1)); % no difference
 
 [p,tbl,stats]= kruskalwallis(bsall(:,:,2)); % on center arm, significantly more theta than heading or interaction
 p = multcompare(stats) 
 
 [p,tbl,stats]= kruskalwallis(bsall(:,:,3)); % no difference
 
 [p,tbl,stats]= kruskalwallis(bsall(:,1,:)); % significantly more theta on center arm than other two arms (not when I take abs value)
 p = multcompare(stats) 
%% glm - averaged abs(betas) no significant differences within or across arms
%Warning: X is ill conditioned, or the model is overparameterized, and
% some coefficients are not identifiable.  You should use caution
% in making predictions. 
cd(dirs.homedir)
d2 = dir('*.mat');
k = 0:1e-5:5e-3;
bsall = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prsopectiveTheta_pfcfr','times_armon_thetaof_headingarm_lap')
    aa = times_armon_thetaof_headingarm_lap;
    pfcfr = prsopectiveTheta_pfcfr;    
    bs = NaN(size(pfc,2),2,3);
    for iarm = 1:3
        
        Y = aa(aa(:,3)==iarm,4:5); % 4 is theta of, 5 is heading arm
        pfc = pfcfr(aa(:,3)==iarm,:); pfc2 = zscore(pfc);
        for icell = 1:size(pfc,2)
            if sum(pfc(:,icell))>0                
                [b,~,~] = glmfit(Y,pfc2(:,icell));                   
%                 b = ridge(pfc(:,icell),Y,k);        
                bs(icell,:,iarm) = abs(b(2:3));
            end
        end
        disp(['Done with arm ' num2str(iarm)])
    end    
    bsall = cat(1,bsall,bs);
    disp(['Done with day ' num2str(id)])
end
%%

[p,tbl,stats]= kruskalwallis(bsall(:,1,:)-bsall(:,2,:)); % no difference across arms
 
 p1 = ranksum(bsall(:,1,1),bsall(:,2,1))
 nanmedian(bsall(:,1,1))-nanmedian(bsall(:,2,1)) % left arm is no different
 
 p3 = ranksum(bsall(:,1,3),bsall(:,2,3))
 nanmedian(bsall(:,1,3))-nanmedian(bsall(:,2,3)) % right arm is no different
 
 p1 = ranksum(bsall(:,1,2),bsall(:,2,2))
 nanmedian(bsall(:,1,2))-nanmedian(bsall(:,2,2)) % center arm is no difference

%% looking at diffferences in theta sweeps
cd(dirs.homedir)
d2 = dir('*.mat');
thdiff = NaN(size(d2,1),3,2);
thdiffall = NaN(size(d2,1),3,2,103);
aaa = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'times_armon_thetaof_headingarm_lap')
    aa = times_armon_thetaof_headingarm_lap;
    for iarm = 1:3
        otharms = setdiff(1:3,iarm);
        for ihead = 1:2        
            a = sum(aa(:,3)==iarm & aa(:,4)==otharms(ihead) & aa(:,5)==otharms(ihead));
            b = sum(aa(:,3)==iarm & aa(:,4)==otharms(setdiff(1:2,ihead)) & aa(:,5)==otharms(ihead));
            thdiff(id,iarm,ihead) = (a-b); %./(a+b);
            for ilap = 1:max(aa(:,end))
                a = sum(aa(:,3)==iarm & aa(:,4)==otharms(ihead) & aa(:,5)==otharms(ihead) & aa(:,6)==ilap);
                b = sum(aa(:,3)==iarm & aa(:,4)==otharms(setdiff(1:2,ihead)) & aa(:,5)==otharms(ihead) & aa(:,6)==ilap);
                if (a+b)>0
                    thdiffall(id,iarm,ihead,ilap) = a-b;
                end
            end
            
        end
    end
    aaa = [aaa;aa(unique(aa(:,6)),[3 5])];
end


signrank(thdiff(:),0,'tail','right')
thdiff2 = nanmean(thdiff,3); signrank(thdiff2(:),0,'tail','right')
thdiff2 = nanmean(thdiff,2); signrank(thdiff2(:),0,'tail','right')
signrank(nanmean(nanmean(thdiff,3),2),0,'tail','right')
%% two way anova on single cell level - not great, confusing results
label = 'RP';
Psa = []; Fsa = []; pSS = []; pSS2 = [];
for id = 1:size(d2,1)
    thisdir = d2(id).name;
    load(thisdir,'prospectiveTheta_toplot','prospectiveTheta_toplot_head',[label '_Cand_sig_modu_include'])
    eval(['Cand = ' label '_Cand_sig_modu_include;'])
    %ask, if it is different between theta sweeops (3rd dim) is it because
    %of heading direction (head last dim). also need to take into account 3 arms
    phead = prospectiveTheta_toplot_head(Cand(:,1)==1,:,:,:,:);
    Fs = NaN(size(phead,1),3,3); Ps = Fs;
    for iarm = 1:3
        for icell = 1:size(phead,1)
            b = squeeze(phead(icell,iarm,:,:,:));
%             bb = reshape(b,[size(b,1) 2*size(b,2)])';
%             [p,tbl] = anova2(bb,size(bb,1)/2,'off'); 
            dat = [];
            for i = 1:size(b,2)
                [x,y] =  find(~isnan(b(:,i,:)));
                if ~isempty(x)
                    for ix = 1:size(x,1)
                        dat = [dat; b(x(ix),i,y(ix)) x(ix) y(ix)];
                    end
                end
            end
            dat(isnan(dat(:,1)),:) = [];
            [p,tbl] = anovan(zscore(dat(:,1)),dat(:,2:3),'model','interaction','display','off');
            Fs(icell,1,iarm) = tbl{2,5};
            Fs(icell,2,iarm) = tbl{3,5};
            Fs(icell,3,iarm) = tbl{4,5};
            Ps(icell,:,iarm) = p;
    %         c = multcompare(stats);
            %column is theta, heading direction is row, interaction
        end
    end
    Fsa = [Fsa; Fs];
    Psa = [Psa; Ps];
    
    pS = []; pS2 = NaN(size(prospectiveTheta_toplot,1),3);
    for icell = 1:size(prospectiveTheta_toplot,1)
        b = permute(squeeze(prospectiveTheta_toplot(icell,:,:,:)),[3 1 2]);
        bb = reshape(b,[size(b,1) 2*size(b,2)]);
%         pS(icell) = kruskalwallis(bb,[],'off');
        pS(icell) = anova1('kruskalwallis',bb,[],'off'); % 38% of cells have significant difference in theta across arms and direction
        
        for iarm = 1:3
            if sum(~isnan(b(:,iarm,1)))>1 && sum(~isnan(b(:,iarm,2)))>1
                pS2(icell,iarm) = ranksum(b(:,iarm,1),b(:,iarm,2)); % 19% have sig difference between R and L when on center
            end
        end
    end
    pSS = [pSS; pS'];
    pSS2 = [pSS2; pS2];
end

%% trying to interpret the resulst of all those two-way anovas
touse = SMI(:,1)==1 & SMI(:,3)==1;
Fsa(isnan(Psa)) = NaN;
touse = true(size(Fsa,1),1);
[p,tbl,stats] = kruskalwallis(nanmean(Fsa,3));
c = multcompare(stats)

p = ranksum(Fsa(touse,1,1),Fsa(touse,2,1))
p = ranksum(Fsa(touse,1,2),Fsa(touse,2,2))
p = ranksum(Fsa(touse,1,3),Fsa(touse,2,3))
ab = Fsa(touse,1,:);
ba = Fsa(touse,2,:);
p = ranksum(ab(:),ba(:));

[p,tbl,stats] = kruskalwallis(nanmean(Fsa(touse,:,:),3)); %interaction sig lower than theta or heading (heading non sig higher)
c = multcompare(stats)
[p,tbl,stats] = kruskalwallis(Fsa(touse,:,1)); % left arm, interaction sig lower than theta or heading (heading non sig higher)
c = multcompare(stats)
[p,tbl,stats] = kruskalwallis(Fsa(touse,:,2)); % center arm, interaction sig lower than theta or heading (theta trending higher)
c = multcompare(stats);
[p,tbl,stats] = kruskalwallis(Fsa(touse,:,3)); % right arm, interaction sig lower than theta or heading (heading non sig higher)
c = multcompare(stats);

[p1,tbl1,stats1] = kruskalwallis(nanmean(Psa,3));
c1 = multcompare(stats1);